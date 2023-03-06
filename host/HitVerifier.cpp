#include "HitVerifier.hpp"
#include "device/PublicDefines.h"
#include "../PhmmReprojection/PhmmReprojection.h"
#include <exception>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <cstring>
//private


HitVerifier::HitVerifier(shared_ptr<FastaVector> fastaVector, shared_ptr<P7HmmList> phmmList)
  :fastaVector(fastaVector),
  phmmList(phmmList){
  this->ssvCellScores = std::make_shared<vector<float>>(SSV_PHMM_VERIFICATION_RANGE);
}



shared_ptr<vector<VerifiedHit>> HitVerifier::verify(shared_ptr<vector<HardwareHitReport>> hits) {
  shared_ptr<vector<VerifiedHit>> verifiedHitList = std::make_shared<vector<VerifiedHit>>();
  for (const auto& hit : *hits) {
    verifyHit(hit, verifiedHitList);
  }

  return verifiedHitList;
}


std::tuple<uint32_t, uint32_t> HitVerifier::getPhmmIndexFromPosition(const HardwareHitReport& hardwareHitReport) {
  const uint32_t numPhmmsInList = this->phmmList->count;

  //set initially to -1, if it ends as -1, the position was outside the possible range of values for the phmmList
  int32_t hitLocatedInPhmmIndex = -1;
  uint32_t globalPhmmStartPosition = 0;

  //find which hit the report is in. The hit gives an exact phmm location, so it can't be in multiple phmms
  for (uint32_t phmmIndex = 0; phmmIndex < numPhmmsInList; phmmIndex++) {
    uint32_t phmmModelLength = this->phmmList->phmms[phmmIndex].header.modelLength;
    uint32_t globalPhmmEndPosition = globalPhmmStartPosition + phmmModelLength;
    if (globalPhmmEndPosition > hardwareHitReport.phmmPosition) {
      hitLocatedInPhmmIndex = phmmIndex;
      break;
    }
    else {
      globalPhmmStartPosition = globalPhmmEndPosition;
    }
  }

  //now we should have the index for the phmm, as well as it's global position we can use to get a local position,
  //but first we need to check for errors. if the phmmIndexForHit is still -1, the given location wasn't in the phmmList's range,
  //so we're going to throw an exception
  if (hitLocatedInPhmmIndex == -1) {
    throw std::domain_error("phmm position was outside the possible range of loaded phmm vectors.");
  }

  if (globalPhmmStartPosition > hardwareHitReport.phmmPosition) {
    throw std::logic_error("global phmm start position was greater than the hit phmm index, which should never happen.");
  }

  uint32_t localPhmmHitPosition = hardwareHitReport.phmmPosition - globalPhmmStartPosition;

  return std::make_pair((uint32_t)hitLocatedInPhmmIndex, localPhmmHitPosition);
}

//can throw std::domain_error
void HitVerifier::verifyHit(const HardwareHitReport& hardwareHitReport, shared_ptr<vector<VerifiedHit>> verifiedHitList) {

  //find which phmm the hit was found in, and the local position of the hit.
  uint32_t phmmIndex;
  uint32_t localPhmmHitPosition;
  std::tie(phmmIndex, localPhmmHitPosition) = getPhmmIndexFromPosition(hardwareHitReport);

  uint32_t thisPhmmLength = this->phmmList->phmms[phmmIndex].header.modelLength;
  uint32_t localSsvPhmmStartPosition = std::max(0L, (int64_t)localPhmmHitPosition - (int64_t)SSV_VERIFICATION_PRIOR_FLANK_RANGE);
  uint32_t localSsvPhmmEndPosition = std::min(thisPhmmLength, localPhmmHitPosition + SSV_VERIFICATION_POST_FLANK_RANGE);
  uint32_t localPhmmSsvLength = localSsvPhmmEndPosition - localSsvPhmmStartPosition;

  //now we have enough info to get the phmm data needed for the local ssv verification. Now, 
  //figure out which sequence segments might be in the range. 
  //this represents the leftmost cell that could've caused the hit.
  uint32_t sequenceHitRangePosition = (hardwareHitReport.sequencePassIndex * NUM_CELL_PROCESSORS) +
    (hardwareHitReport.sequenceGroupIndex * CELLS_PER_GROUP);

  //total length of all the sequences, concatenated together.
  uint32_t totalSequenceAggregatedLength = this->fastaVector->metadata.data[this->fastaVector->metadata.count - 1].sequenceEndPosition;
  //leftmost cell in the DP-matrix to perform reference SSV on, clamped to zero (aka, a ReLU)
  uint32_t sequencePossibleStartPosition = std::max(0L, (int64_t)sequenceHitRangePosition - (int64_t)SSV_VERIFICATION_PRIOR_FLANK_RANGE);
  //this end position adds CELLS_PER_GROUP, since any cell in the group could've caused the hit.
  uint32_t sequencePossibleEndPosition = std::min(totalSequenceAggregatedLength,
    sequencePossibleStartPosition + CELLS_PER_GROUP + SSV_VERIFICATION_POST_FLANK_RANGE);


  //iterate through the sequences, and determine which sequences overlap this range
  uint32_t thisSequenceStartPosition = 0;
  for (uint32_t sequenceIndex = 0; sequenceIndex < this->fastaVector->metadata.count; sequenceIndex++) {
    //end position is the position 1 after index of the last character in the sequence
    uint32_t thisSequenceEndPosition = this->fastaVector->metadata.data[sequenceIndex].sequenceEndPosition;

    //if there is overlap, it's worth checking for an actual hit.
    bool sequenceOverlapsWithPossibleHit = 
    (thisSequenceStartPosition <= sequencePossibleEndPosition) &&
      (thisSequenceEndPosition > sequencePossibleStartPosition);
    
    if (sequenceOverlapsWithPossibleHit) {

      //these sequences overlap, the hit might be in this sequence. 
      //the phmm will need local coordinates, but the sequence will need global coordinants.
      //keep the ssv range in the range of the actual sequence.
      const uint32_t ssvSequenceRangeBegin = std::max(thisSequenceStartPosition, sequencePossibleStartPosition);
      //ssvSequenceRangeEnd is exclusive (this value may be out of bounds)
      const uint32_t ssvSequenceRangeEnd = std::min(thisSequenceEndPosition, sequencePossibleEndPosition);

      this->verifyWithReferenceSsv(phmmIndex, sequenceIndex, localSsvPhmmStartPosition, localPhmmSsvLength, ssvSequenceRangeBegin, ssvSequenceRangeEnd);
    }
  }
}

//TODO: debug-check this, and make sure the phmm scaling is right, we might have to rescale
void HitVerifier::verifyWithReferenceSsv(const uint32_t hitLocatedInPhmmNumber, const uint32_t sequenceNumber, const uint32_t localSsvPhmmStartPosition,
  const uint32_t localPhmmSsvLength, const uint32_t ssvSequenceRangeBegin, const uint32_t ssvSequenceRangeEnd) {

  const float ssvThresholdScore = 256.0f; //the scores in the PHMM struct are prescaled for this 

  const uint32_t ssvSequenceRangeLength = ssvSequenceRangeEnd - ssvSequenceRangeBegin;
  //reinitialize the ssvCellScores array
  std::memset(this->ssvCellScores.get()->data(), 0, sizeof(this->ssvCellScores));

  uint8_t phmmCardinality = p7HmmGetAlphabetCardinality(&this->phmmList->phmms[hitLocatedInPhmmNumber]);
  const float* phmmRange = &this->phmmList->phmms[hitLocatedInPhmmNumber].model.matchEmissionScores[localSsvPhmmStartPosition * phmmCardinality];

  for (uint32_t sequencePosition = ssvSequenceRangeBegin; sequencePosition < ssvSequenceRangeEnd; sequencePosition++) {
    char sequenceSymbol = this->fastaVector->sequence.charData[sequencePosition];
    uint8_t sequenceSymbolAsIndex;
    switch (sequenceSymbol) {
    case 'a': case 'A': sequenceSymbolAsIndex = 0; break;
    case 'c': case 'C': sequenceSymbolAsIndex = 1; break;
    case 'g': case 'G': sequenceSymbolAsIndex = 2; break;
    default:            sequenceSymbolAsIndex = 3;
    }

    auto cellScores = *this->ssvCellScores;
    for (uint32_t phmmPosition = localPhmmSsvLength - 1; phmmPosition > 0; phmmPosition--) {
      float matchScore = phmmRange[(phmmPosition * phmmCardinality) + sequenceSymbolAsIndex];
      // const uint32_t positionInSsvCellScores = sequencePosition - ssvSequenceRangeBegin;
      cellScores[phmmPosition] = cellScores[phmmPosition - 1] + matchScore;
    }
    //now compute the first cell
    float matchScore = phmmRange[sequenceSymbolAsIndex];
    cellScores[0] = matchScore;


    for (uint32_t ssvCell = 0; ssvCell < localPhmmSsvLength; ssvCell++) {
      if (cellScores[ssvCell] >= ssvThresholdScore) {

        std::vector<VerifiedHit> verifiedHitLocations;
        uint64_t sequencePosition = ssvSequenceRangeBegin + ssvCell;
        uint32_t sequenceIndex = sequenceNumber;
        uint32_t phmmPosition = localPhmmSsvLength + phmmPosition;
        uint32_t phmmIndex = hitLocatedInPhmmNumber;
        VerifiedHit verifiedHit(sequencePosition, sequenceIndex, phmmPosition, phmmIndex);
        verifiedHitLocations.push_back(verifiedHit);
      }
    }
  }
}

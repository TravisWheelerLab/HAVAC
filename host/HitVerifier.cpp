#include "HitVerifier.hpp"
#include "device/PublicDefines.h"
#include <exception>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include "../../PhmmReprojection/PhmmReprojection.h"

//private


HitVerifier::HitVerifier(shared_ptr<FastaVector>& fastaVector, shared_ptr<P7HmmList>& phmmList)
  :fastaVector(fastaVector),
  phmmList(phmmList) {
  this->ssvCellScores = std::make_unique<float[]>(SSV_VERIFICATION_RANGE);
}

void HitVerifier::setFastaVector(shared_ptr<FastaVector>& fastaVector) {
  this->fastaVector = fastaVector;
}

void HitVerifier::setPhmm(shared_ptr<P7HmmList>& phmmList) {
  this->phmmList = phmmList;
}


void HitVerifier::clear() {
  this->verifiedHitLocations.clear();
}



std::tuple<std::shared_ptr<std::vector<SsvVerifiedHit>>,
  std::shared_ptr<std::vector<HavacHardwareHit>>>
  HitVerifier::verify(std::shared_ptr<std::vector<HavacHardwareHitReport>> hits) {
  for (const auto& hit : hits) {
    verifyHit(hit);
  }
}


std::tuple<int32_t, uint32_t>HitVerifier::getPhmmIndexFromPosition(const struct HavacHardwareHit& hardwareHitReport) {
  const uint32_t numPhmmsInList  this->phmmList->count;
  const uint32_t phmmModelLength = this->phmmList->phmms[phmmIndex].header.modelLength;

  //set initially to -1, if it ends as -1, the position was outside the possible range of values for the phmmList
  int32_t hitLocatedInPhmmIndex = -1;
  uint32_t globalPhmmStartPosition = 0;

  //find which hit the report is in. The hit gives an exact phmm location, so it can't be in multiple phmms
  for (uint32_t phmmIndex = 0; phmmIndex < numPhmmsInList; phmmIndex++) {
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
  if (hitLocatedInPhmmNumber == -1) {
    throw std::domain_error("phmm position was outside the possible range of loaded phmm vectors.");
  }

  if (globalPhmmStartPosition > hardwareHitReport.phmmIndex) {
    throw std::logic_error("global phmm start position was greater than the hit phmm index, which should never happen.");
  }

  uint32_t localPhmmHitIndex = hardwareHitReport.phmmIndex - globalPhmmStartPosition;

  return std::make_pair(hitLocatedInPhmmIndex, localPhmmHitIndex);
}

//can throw std::domain_error
bool HitVerifier::verifyHit(const struct HavacHardwareHit& hardwareHitReport) {

  //find which phmm the hit was found in, and the local position of the hit.
  int32_t phmmIndex;
  uint32_t localPhmmHitPosition;
  std::tie(phmmIndex, localPhmmHitPosition) = getPhmmIndexFromPosition(hardwareHitReport);

  uint32_t thisPhmmLength = this->phmmList->phmms[phmmIndex].header.modelLength;
  uint32_t localSsvPhmmStartPosition = std::max(0, (int64_t)localPhmmHitPosition - (int64_t)SSV_VERIFICATION_PRIOR_FLANK_RANGE);
  uint32_t localSsvPhmmEndPosition = std::min(thisPhmmLength, localPhmmHitPosition + SSV_VERIFICATION_POST_FLANK_RANGE);
  uint32_t localPhmmSsvLength = localSsvPhmmEndPosition - localSsvPhmmStartPosition;

  //now we have enough info to get the phmm data needed for the local ssv verification. Now, 
  //figure out which sequence segments might be in the range.
  uint32_t sequenceHitRangePosition = (hardwareHitReport.sequencePassIndex * NUM_CELL_PROCESSORS) +
    (hardwareHitReport.sequenceGroupIndex * CELLS_PER_GROUP);
  uint32_t totalSequenceAggregatedLength = this->fastaVector->metadata.data[this->fastaVector->metadata.count - 1].sequenceEndPosition;
  uint32_t sequencePossibleStartPosition = std::max(0, (int64_t)sequenceHitRangePosition - (int64_t)SSV_VERIFICATION_PRIOR_FLANK_RANGE);
  uint32_t sequencePossibleEndPosition = std::min(totalSequenceAggregatedLength,
    sequencePossibleStartPosition + (CELLS_PER_GROUP - 1) + SSV_VERIFICATION_POST_FLANK_RANGE);


  //iterate through the sequences, and determine which sequences overlap this range
  uint32_t thisSequenceStartPosition = 0;
  for (uint32_t sequenceIndex = 0; sequenceIndex < this->fastaVector->metadata.count; sequenceIndex++) {
    uint32_t thisSequenceEndPosition = this->fastaVector->metadata.data[sequenceIndex].sequenceEndPosition;

    //if there is overlap, it's worth checking for an actual hit.
    bool sequenceOverlapsWithPossibleHit = (thisSequenceStartPosition <= sequencePossibleEndPosition) &&
      (thisSequenceEndPosition >= sequencePossibleStartPosition);
    if (sequenceOverlapsWithPossibleHit) {

      //these sequences overlap, the hit might be in this sequence. 
      //the phmm will need local coordinates, but the sequence will need global coordinants.
      //keep the ssv range in the range of the actual sequence.
      const uint32_t ssvSequenceRangeBegin = std::max(thisSequenceStartPosition, sequencePossibleStartPosition);
      const uint32_t ssvSequenceRangeEnd = std::min(thisSequenceEndPosition, sequencePossibleEndPosition);

      this->verifySsv(phmmIndex, localSsvPhmmStartPosition, localPhmmSsvLength, ssvSequenceRangeBegin, ssvSequenceRangeEnd);
    }
  }
}

void HitVerifier::verifyWithReferenceSsv(const uint32_t hitLocatedInPhmmNumber, const uint32_t sequenceNumber, const uint32_t localSsvPhmmStartPosition,
  const uint32_t localPhmmSsvLength, const uint32_t ssvSequenceRangeBegin, const uint32_t ssvSequenceRangeEnd) {

  const float ssvThresholdScore = 256.0f; //the scores in the PHMM struct are prescaled for this 

  const uint32_t ssvSequenceRangeLength = ssvSequenceRangeEnd - ssvSequenceRangeBegin;
  //reinitialize the ssvCellScores array
  std::memset(this->ssvCellScores, 0, sizeof(this->ssvCellScores));

  uint8_t phmmCardinality = p7HmmGetAlphabetCardinality(&this->phmmList->phmms[hitLocatedInPhmmNumber]);
  const float* phmmRange = &this->phmmList->phmms[hitLocatedInPhmmNumber].model.matchEmissionScores[localSsvPhmmStartPosition * phmmCardinality];

  for (uint32_t phmmPosition = 0; phmmPosition < localPhmmSsvLength; phmmPosition++) {

    for (uint32_t sequencePosition = ssvSequenceRangeEnd - 1; sequencePosition > ssvSequenceRangeBegin; sequencePosition--) {
      char sequenceSymbol = this->fastaVector->sequence.charData[sequencePosition];
      uint8_t sequenceSymbolAsIndex;
      switch (sequenceSymbol) {
      case 'a': case 'A': sequenceSymbolAsIndex = 0; break;
      case 'c': case 'C': sequenceSymbolAsIndex = 1; break;
      case 'g': case 'G': sequenceSymbolAsIndex = 2; break;
      default:            sequenceSymbolAsIndex = 3;
      }

      float phmmValue = phmmRange[(phmmPosition * phmmCardinality) + sequenceSymbolAsIndex];
      const uint32_t positionInSsvCellScores = sequencePosition - ssvSequenceRangeBegin;
      this->ssvCellScores[positionInSsvCellScores] = this->ssvCellScores[positionInSsvCellScores - 1] + phmmValue;


    }
    //now compute the first cell
    char sequenceSymbol = this->fastaVector->sequence.charData[ssvSequenceRangeBegin];
    uint8_t sequenceSymbolAsIndex;
    switch (sequenceSymbol) {
    case 'a': case 'A': sequenceSymbolAsIndex = 0; break;
    case 'c': case 'C': sequenceSymbolAsIndex = 1; break;
    case 'g': case 'G': sequenceSymbolAsIndex = 2; break;
    default:            sequenceSymbolAsIndex = 3;
    }
    float phmmValue = phmmRange[(phmmPosition * phmmCardinality) + sequenceSymbolAsIndex];
    this->ssvCellScores[0] = phmmValue;


    for (uint32_t ssvCell = 0; ssvCell < ssvSequenceRangeLength; ssvCell++) {
      if (this->ssvCellScores[ssvCell] >= ssvThresholdScore) {

        std::vector<struct SsvVerifiedHitLocation> verifiedHitLocations;
        struct SsvVerifiedHitLocation verifiedHit;
        verifiedHit.phmmIndex = hitLocatedInPhmmNumber;
        verifiedHit.phmmPosition = localPhmmSsvLength + phmmPosition;
        verifiedHit.sequenceIndex = sequenceNumber;
        verifiedHit.sequencePosition = ssvSequenceRangeBegin + ssvCell;
        verifiedHitLocations.push_back(verifiedHit);
      }
    }
  }
}

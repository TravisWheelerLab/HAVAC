#include "HitVerifier.hpp"
#include "device/PublicDefines.h"
#include "../PhmmReprojection/PhmmReprojection.h"
#include <exception>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <cstring>
//private


HitVerifier::HitVerifier(FastaVector* fastaVector, P7HmmList* phmmList)
  :fastaVector(fastaVector),
  phmmList(phmmList) {
  this->ssvCellScores = std::make_shared<vector<float>>(SSV_PHMM_VERIFICATION_RANGE);
}



shared_ptr<vector<VerifiedHit>> HitVerifier::verify(shared_ptr<vector<HardwareHitReport>> hits, const float desiredPValue) {
  shared_ptr<vector<VerifiedHit>> verifiedHitList = std::make_shared<vector<VerifiedHit>>();
  for (const auto& hit : *hits) {
    verifyHit(hit, verifiedHitList, desiredPValue);
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


void HitVerifier::verifyHitForGroup(uint32_t localPhmmPosition, uint32_t phmmIndex, uint32_t sequenceSegmentIndex, uint32_t groupIndex,
  shared_ptr<vector<VerifiedHit>> verifiedHitList, const float desiredPValue) {
  std::cout << "phmm: "<< localPhmmPosition<<"."<<phmmIndex<<", ssi:"<< sequenceSegmentIndex<< "group: " << groupIndex << "." << std::endl;
  //find which phmm the hit was found in, and the local position of the hit.

  uint32_t sequenceHitRangePosition = (sequenceSegmentIndex * NUM_CELL_PROCESSORS) +
    (groupIndex * CELLS_PER_GROUP);

  //total length of all the sequences, concatenated together.
  uint32_t totalSequenceAggregatedLength = this->fastaVector->metadata.data[this->fastaVector->metadata.count - 1].sequenceEndPosition;
  //leftmost cell in the DP-matrix to perform reference SSV on, clamped to zero
  uint32_t sequencePossibleStartPosition = std::max(0L, (int64_t)sequenceHitRangePosition - (int64_t)SSV_VERIFICATION_PRIOR_FLANK_RANGE);
  //this end position adds CELLS_PER_GROUP, since any cell in the group could've caused the hit.
  uint32_t sequencePossibleEndPosition = std::min(totalSequenceAggregatedLength,
    sequenceHitRangePosition + CELLS_PER_GROUP + SSV_VERIFICATION_POST_FLANK_RANGE);


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
      std::cout << "checking for hit in sequence #" << sequenceIndex << std::endl;
      //the sequence overlaps with the possible hit area, the hit might be in this sequence. 
      //the phmm will need local coordinates, but the sequence will need global coordinants.
      //keep the ssv range in the range of the actual sequence.
      const uint32_t ssvSequenceRangeBegin = std::max(thisSequenceStartPosition, sequencePossibleStartPosition);
      //ssvSequenceRangeEnd is exclusive (this value may be out of bounds)
      const uint32_t ssvSequenceRangeEnd = std::min(thisSequenceEndPosition, sequencePossibleEndPosition);
      std::cout << "[" << ssvSequenceRangeBegin << "," << ssvSequenceRangeEnd << "]" << std::endl;
      this->verifyWithReferenceSsv(phmmIndex, sequenceIndex, localPhmmPosition,
        ssvSequenceRangeBegin, ssvSequenceRangeEnd, desiredPValue, verifiedHitList);
    }
  }
}

//can throw std::domain_error
void HitVerifier::verifyHit(const HardwareHitReport& hardwareHitReport, shared_ptr<vector<VerifiedHit>> verifiedHitList, const float desiredPValue) {
  uint32_t groupsReportingHit = hardwareHitReport.sequenceGroupIndex;
  while (groupsReportingHit) {
    //find the next bit that's set in the sequence group index. while this will almost always only 
    //have 1 bit set, we need to handle the possibility of multiple groups set in the same report.
    uint32_t groupIndex = _bit_scan_forward(hardwareHitReport.sequenceGroupIndex);
    //remove the group we're about to validaten from the sequenceGroupIndex data
    groupsReportingHit ^= 1 << groupIndex;

    uint32_t phmmIndex;
    uint32_t localPhmmHitPosition;
    std::tie(phmmIndex, localPhmmHitPosition) = getPhmmIndexFromPosition(hardwareHitReport);
    verifyHitForGroup(localPhmmHitPosition, phmmIndex, hardwareHitReport.sequencePassIndex, groupIndex,
      verifiedHitList, desiredPValue);
  }
}



void HitVerifier::verifyWithReferenceSsv(const uint32_t hitLocatedInPhmmNumber, const uint32_t sequenceNumber,
  const uint32_t hitOccurredAtPhmmIndex, const uint32_t possibleSequenceIndexStart, const uint32_t possibleSequenceIndexEnd,
  const float desiredPValue, shared_ptr<vector<VerifiedHit>> verifiedHitList) {
  int16_t debugBestScore = -1;
  uint_fast16_t SSV_THRESHOLD = 256;
  P7Hmm* phmmPointer = &this->phmmList->phmms[hitLocatedInPhmmNumber];
  float ssvScoreMultiplier = generateScoreMultiplierForPhmmScore(phmmPointer, desiredPValue);

  uint8_t phmmCardinality = p7HmmGetAlphabetCardinality(phmmPointer);

  //iterate over each cell that could have caused the hit
  for (uint32_t sequenceConsiderationIndex = possibleSequenceIndexStart; sequenceConsiderationIndex < possibleSequenceIndexEnd;
    sequenceConsiderationIndex++) {
    //the length of the diagonal to accumulate score on. this could be shorter if that length would go past the start of the phmm or sequence  
    const uint32_t MAX_SSV_WALKBACK_DISTANCE = 64;
    uint32_t diagonalWalkbackLength = std::min({ MAX_SSV_WALKBACK_DISTANCE, hitOccurredAtPhmmIndex, sequenceConsiderationIndex });
    int16_t accumulatedScore = 0;
    bool foundThresholdHitOnDiagonal = false;
    for (uint32_t ssvWalkbackIndex = 0; ssvWalkbackIndex < diagonalWalkbackLength; ssvWalkbackIndex++) {
      const uint32_t walkbackPhmmIndex = hitOccurredAtPhmmIndex - ssvWalkbackIndex;
      float* phmmVectorAsFloat = &phmmPointer->model.matchEmissionScores[walkbackPhmmIndex * phmmCardinality];

      const uint32_t walkbackSequenceIndex = sequenceConsiderationIndex - ssvWalkbackIndex;
      char sequenceSymbol = this->fastaVector->sequence.charData[walkbackSequenceIndex];
      uint8_t sequenceSymbolAsIndex;
      switch (sequenceSymbol) {
      case 'a': case 'A': sequenceSymbolAsIndex = 0; break;
      case 'c': case 'C': sequenceSymbolAsIndex = 1; break;
      case 'g': case 'G': sequenceSymbolAsIndex = 2; break;
      default:            sequenceSymbolAsIndex = 3;
      }

      float unprojectedMatchScore = phmmVectorAsFloat[sequenceSymbolAsIndex];
      float projectedPhmmMatchScore = projectPhmmScoreWithMultiplier(unprojectedMatchScore, ssvScoreMultiplier);
      int8_t matchScoreAsInt = projectedPhmmMatchScore;
      accumulatedScore += matchScoreAsInt;


      if (accumulatedScore > debugBestScore) {
        debugBestScore = accumulatedScore;
      }
      //the following method of saturated addition and threshold checking should cause
      //no branching instructions on x86 GCC with -03 according to CompilerExplorer
      // first, we check for a threshold hit, and OR it into the boolean. Then,
      // we reset on either underflow or overflow
      foundThresholdHitOnDiagonal |= accumulatedScore >= 256;
      if (accumulatedScore < 0 || accumulatedScore >= 256) {
        accumulatedScore = 0;
      }
    }
    if (foundThresholdHitOnDiagonal) {
      VerifiedHit hit(sequenceConsiderationIndex, sequenceNumber, hitOccurredAtPhmmIndex, hitLocatedInPhmmNumber);
      verifiedHitList->push_back(hit);
    }
  }
  std::cout << "best acc: " << debugBestScore << std::endl;
}
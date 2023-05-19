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
  for (auto& hit : *hits) {
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
  //find which phmm the hit was found in, and the local position of the hit.

  uint64_t hitGroupStartCellIndex = (sequenceSegmentIndex * NUM_CELL_PROCESSORS) +
    (groupIndex * CELLS_PER_GROUP);


  //total length of all the sequences, concatenated together.
  uint64_t totalSequenceAggregatedLength = this->fastaVector->metadata.data[this->fastaVector->metadata.count - 1].sequenceEndPosition;

  //this end position adds CELLS_PER_GROUP, since any cell in the group could've caused the hit.
  uint64_t hitGroupEndCellIndex = std::min(totalSequenceAggregatedLength,
    hitGroupStartCellIndex + CELLS_PER_GROUP);


  //iterate through the sequences, and determine which sequences overlap this range
  uint64_t localSequenceStartPosition = 0;
  for (uint32_t sequenceIndex = 0; sequenceIndex < this->fastaVector->metadata.count; sequenceIndex++) {
    //end position is the position 1 after index of the last character in the sequence
    uint64_t localSequenceEndPosition = this->fastaVector->metadata.data[sequenceIndex].sequenceEndPosition;

    //if there is overlap, it's worth checking for an actual hit.
    bool sequenceOverlapsWithPossibleHit =
      (localSequenceStartPosition <= hitGroupEndCellIndex) &&
      (localSequenceEndPosition > hitGroupStartCellIndex);

    if (sequenceOverlapsWithPossibleHit) {

      //the sequence overlaps with the possible hit area, the hit might be in this sequence. 
      //the phmm will need local coordinates, but the sequence will need global coordinants.
      //keep the ssv range in the range of the actual sequence.
      const uint64_t ssvSequenceRangeBegin = std::max(localSequenceStartPosition, hitGroupStartCellIndex);
      //ssvSequenceRangeEnd is exclusive (this value may be out of bounds)
      const uint64_t ssvSequenceRangeEnd = std::min(localSequenceEndPosition, hitGroupEndCellIndex);
      this->verifyWithReferenceSsv(phmmIndex, sequenceIndex, localPhmmPosition,
        ssvSequenceRangeBegin, ssvSequenceRangeEnd, desiredPValue, verifiedHitList);
    }
    //set the start position using the end of the sequence we just looked at
    localSequenceStartPosition = localSequenceEndPosition;
  }
}

//can throw std::domain_error
void HitVerifier::verifyHit(const HardwareHitReport& hardwareHitReport, shared_ptr<vector<VerifiedHit>> verifiedHitList, const float desiredPValue) {
  uint8_t groupsReportingHits[NUM_CELL_GROUPS / 8];
  std::memcpy(groupsReportingHits, hardwareHitReport.groupHitBits, GROUP_HIT_BITS_NUM_BYTES);

  for (uint32_t thisGroupHitByte = 0; thisGroupHitByte < GROUP_HIT_BITS_NUM_BYTES; thisGroupHitByte++) {

    while (groupsReportingHits[thisGroupHitByte]) {
      //find the next bit that's set in the sequence group index. while this will almost always only 
      //have 1 bit set, we need to handle the possibility of multiple groups set in the same report.
      uint32_t groupIndex = _bit_scan_forward(groupsReportingHits[thisGroupHitByte]);
      //remove the group we're about to validaten from the sequenceGroupIndex data
      groupsReportingHits[thisGroupHitByte] ^= 1 << groupIndex;

      groupIndex += (thisGroupHitByte * 8);

      uint32_t phmmIndex;
      uint32_t localPhmmHitPosition;
      std::tie(phmmIndex, localPhmmHitPosition) = getPhmmIndexFromPosition(hardwareHitReport);
      verifyHitForGroup(localPhmmHitPosition, phmmIndex, hardwareHitReport.sequencePassIndex, groupIndex,
        verifiedHitList, desiredPValue);
    }
  }
}

void HitVerifier::verifyWithReferenceSsv(const uint32_t hitLocatedInPhmmNumber, const uint32_t sequenceNumber,
  const uint32_t hitOccurredAtPhmmIndex, const uint64_t possibleSequenceIndexStart, const uint64_t possibleSequenceIndexEnd,
  const float desiredPValue, shared_ptr<vector<VerifiedHit>> verifiedHitList) {
  P7Hmm* phmmPointer = &this->phmmList->phmms[hitLocatedInPhmmNumber];
  const uint_fast16_t SSV_THRESHOLD = 256;
  const float ssvScoreMultiplier = generateScoreMultiplierForPhmmScore(phmmPointer, desiredPValue);
  const uint8_t phmmCardinality = p7HmmGetAlphabetCardinality(phmmPointer);
  const uint64_t phmmLength = phmmPointer->header.modelLength;


  uint64_t localSequenceStartIndex = sequenceNumber == 0 ?
    0 : this->fastaVector->metadata.data[sequenceNumber - 1].sequenceEndPosition;
  uint64_t localSequenceEndIndex = this->fastaVector->metadata.data[sequenceNumber].sequenceEndPosition;
  uint64_t localSequenceLength = localSequenceEndIndex - localSequenceStartIndex;


  //iterate over each cell that could have caused the hit
  for (uint64_t sequenceConsiderationIndex = possibleSequenceIndexStart; sequenceConsiderationIndex < possibleSequenceIndexEnd;
    sequenceConsiderationIndex++) {

    const uint64_t localSequenceConsiderationIndex = sequenceConsiderationIndex - localSequenceStartIndex;
    const uint64_t remainingVectorsInPhmm = phmmLength - hitOccurredAtPhmmIndex;
    const uint64_t remainingSymbolsInSequence = localSequenceEndIndex - localSequenceConsiderationIndex;

    const uint32_t maximumSsvWalkback = std::min({ (uint64_t)hitOccurredAtPhmmIndex, localSequenceConsiderationIndex });
    const uint32_t maximumSsvWalkforward = std::min({ (uint64_t)(phmmPointer->header.modelLength - hitOccurredAtPhmmIndex),
    this->fastaVector->metadata.data[sequenceNumber].sequenceEndPosition - sequenceConsiderationIndex });

    float* phmmVectorAsFloat = phmmPointer->model.matchEmissionScores;
    const char* sequenceData = &this->fastaVector->sequence.charData[localSequenceStartIndex];
    verifySsvWalkback(hitLocatedInPhmmNumber, sequenceNumber, hitOccurredAtPhmmIndex, localSequenceConsiderationIndex, maximumSsvWalkback,
      maximumSsvWalkforward, phmmVectorAsFloat, sequenceData, ssvScoreMultiplier, verifiedHitList);

  }
}

void HitVerifier::verifySsvWalkback(const uint32_t hitLocatedInPhmmNumber, const uint32_t sequenceNumber,
  const uint32_t phmmPosition, const uint64_t sequencePosition, const uint64_t maxWalkbackLength,
  const uint64_t maxWalkforwardLength, const float* phmmVectorAsFloat, const char* sequenceData, const float ssvScoreMultiplier,
  shared_ptr<vector<VerifiedHit>> verifiedHitList) {

  const uint32_t PHMM_CARDINALITY = 4;
  int16_t accumulatedScore = 0;
  uint8_t highestScoreSeen = 0;
  for (uint32_t ssvWalkbackIndex = 0; ssvWalkbackIndex <= maxWalkbackLength; ssvWalkbackIndex++) {

    char sequenceSymbol = sequenceData[sequencePosition - ssvWalkbackIndex];

    uint8_t sequenceSymbolAsIndex;
    switch (sequenceSymbol) {
    case 'a': case 'A': sequenceSymbolAsIndex = 0; break;
    case 'c': case 'C': sequenceSymbolAsIndex = 1; break;
    case 'g': case 'G': sequenceSymbolAsIndex = 2; break;
    default:            sequenceSymbolAsIndex = 3;
    }

    float unprojectedMatchScore = phmmVectorAsFloat[((phmmPosition - ssvWalkbackIndex) * PHMM_CARDINALITY) + sequenceSymbolAsIndex];
    float projectedPhmmMatchScore = projectPhmmScoreWithMultiplier(unprojectedMatchScore, ssvScoreMultiplier);
    int8_t matchScoreAsInt = projectedPhmmMatchScore;
    uint8_t originalScore = accumulatedScore;
    accumulatedScore += (int16_t)matchScoreAsInt;

    if (accumulatedScore > highestScoreSeen) {
      highestScoreSeen = accumulatedScore;
    }


    //if we've accumulated to a score of 0 or lower, then nothing before this cell could've caused a hit.
    if (accumulatedScore <= 0) {
      return;
    }
    else if (accumulatedScore >= 256) {
      VerifiedHit hit(sequencePosition, sequenceNumber, phmmPosition, hitLocatedInPhmmNumber);
      verifiedHitList->push_back(hit);
      std::cout << "\t WALKBACK HIT: seq#" << sequenceNumber << " (" << sequencePosition << "), p#" << hitLocatedInPhmmNumber << " (" << phmmPosition << ")" << std::endl;
      return;
    }
  }

  //the bgest way we could extend this is by using the highest score we've seen so far during the walkback.
  accumulatedScore = highestScoreSeen;
  std::cout << "walkback seems to have hit the start of the phmm or sequence, beginning walk forward with residual score" << (int)accumulatedScore << std::endl;

  //now, if we still haven't hit 256 or 0, walk forward past the place that hit to see if we'd get a hit on the diagonal.
  //the accumulator isn't cleared here, since we need to collect more score to see a hit
  for (size_t ssvWalkforwardIndex = 0; ssvWalkforwardIndex < maxWalkforwardLength; ssvWalkforwardIndex++) {
    char sequenceSymbol = sequenceData[sequencePosition + ssvWalkforwardIndex];
    uint8_t sequenceSymbolAsIndex;
    switch (sequenceSymbol) {
    case 'a': case 'A': sequenceSymbolAsIndex = 0; break;
    case 'c': case 'C': sequenceSymbolAsIndex = 1; break;
    case 'g': case 'G': sequenceSymbolAsIndex = 2; break;
    default:            sequenceSymbolAsIndex = 3;
    }

    float unprojectedMatchScore = phmmVectorAsFloat[((phmmPosition + ssvWalkforwardIndex) * PHMM_CARDINALITY) + sequenceSymbolAsIndex];
    float projectedPhmmMatchScore = projectPhmmScoreWithMultiplier(unprojectedMatchScore, ssvScoreMultiplier);
    int8_t matchScoreAsInt = projectedPhmmMatchScore;
    uint8_t originalScore = accumulatedScore;
    accumulatedScore += (int16_t)matchScoreAsInt;

    //if we've accumulated to a score of 0 or lower, then no use extending the hit further
    if (accumulatedScore <= 0) {
      return;
    }
    else if (accumulatedScore >= 256) {
      VerifiedHit hit(sequencePosition + ssvWalkforwardIndex, sequenceNumber, phmmPosition + ssvWalkforwardIndex, hitLocatedInPhmmNumber);
      verifiedHitList->push_back(hit);
      return;
    }
  }
}

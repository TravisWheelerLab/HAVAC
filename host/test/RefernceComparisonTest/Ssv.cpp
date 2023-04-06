#include "Ssv.hpp"
#include "../../../PhmmReprojection/PhmmReprojection.h"
#include <cstring>


shared_ptr<vector<ReferenceSsvHit>> HitsFromSsv(FastaVector* fastaVector,
  P7HmmList* phmmList, const float desiredPValue) {
  shared_ptr<vector<ReferenceSsvHit>> hitsFromReference = make_shared<vector<ReferenceSsvHit>>();
  const uint8_t cardinality = p7HmmGetAlphabetCardinality(&phmmList->phmms[0]);
  const uint32_t sequenceLength = fastaVector->sequence.count;
  const uint32_t phmmLength = phmmList->phmms[0].header.modelLength;

  float scoreMultiplier = generateScoreMultiplierForPhmmScore(&phmmList->phmms[0], desiredPValue);

  shared_ptr<vector<uint8_t>> cellScores = make_shared<vector<uint8_t>>(phmmLength);
  std::memset(cellScores->data(), 0, phmmLength * sizeof(uint8_t));

  for (uint32_t sequencePosition = 0; sequencePosition < sequenceLength; sequencePosition++) {
    char symbolFromSequence = fastaVector->sequence.charData[sequencePosition];
    uint8_t symbolEncoding;
    switch (symbolFromSequence) {
    case 'a': case 'A': symbolEncoding = 0; break;
    case 'c': case 'C': symbolEncoding = 1; break;
    case 'g': case 'G': symbolEncoding = 2; break;
    default:            symbolEncoding = 3;
    }

    for (uint32_t phmmPosition = phmmLength - 1; phmmPosition > 0; phmmPosition--) {
      float matchScoreFromPhmmFile = phmmList->phmms[0].model.matchEmissionScores[(phmmPosition * cardinality) + symbolEncoding];
      int8_t projectedMatchScore = projectPhmmScoreWithMultiplier(matchScoreFromPhmmFile, scoreMultiplier);
      int16_t cellSumResult = (int16_t)cellScores->data()[phmmPosition] + (int16_t)projectedMatchScore;

      if (cellSumResult > 256) {
        ReferenceSsvHit hit = ReferenceSsvHit();
        hit.sequenceNumber = 0;
        hit.phmmNumber = 0;
        hit.sequencePosition = sequencePosition;
        hit.phmmPosition = phmmPosition;
        hitsFromReference->push_back(hit);

        cellSumResult = 0;
      }
      cellScores->data()[phmmPosition] = cellSumResult;

      //now, handle the first cell in the row. it should be impossible to generate a hit here, 
      // as it can't ever accumulate enough score for a hit yet.
      matchScoreFromPhmmFile = phmmList->phmms[0].model.matchEmissionScores[symbolEncoding];
      projectedMatchScore = projectPhmmScoreWithMultiplier(matchScoreFromPhmmFile, scoreMultiplier);

      cellSumResult = projectedMatchScore >= 0 ? projectedMatchScore : 0;
      cellScores->data()[0] = cellSumResult;
    }
  }

  return hitsFromReference;
}

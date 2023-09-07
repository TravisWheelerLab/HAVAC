#include "Ssv.hpp"
#include "../../PhmmReprojection/PhmmReprojection.h"
#include <cstring>
#include <iostream>
#include <math.h>


shared_ptr<vector<ReferenceSsvHit>> HitsFromSsv(FastaVector* fastaVector,
  P7HmmList* phmmList, const float desiredPValue) {
  shared_ptr<vector<ReferenceSsvHit>> hitsFromReference = make_shared<vector<ReferenceSsvHit>>();
  const uint8_t cardinality = p7HmmGetAlphabetCardinality(&phmmList->phmms[0]);

  for (uint32_t phmmNumber = 0; phmmNumber < phmmList->count; phmmNumber++) {
    float* phmmMatchEmissionScores = phmmList->phmms[phmmNumber].model.matchEmissionScores;
    float scoreMultiplier = findThreshold256ScalingFactor(&phmmList->phmms[phmmNumber], desiredPValue);
    const uint32_t phmmLength = phmmList->phmms[phmmNumber].header.modelLength;
    shared_ptr<vector<uint8_t>> cellScores = make_shared<vector<uint8_t>>(phmmLength);

    for (uint32_t sequenceNumber = 0; sequenceNumber < fastaVector->metadata.count;sequenceNumber++) {
      std::memset(cellScores->data(), 0, phmmLength * sizeof(uint8_t));
      uint32_t sequenceLength;
      uint32_t sequenceStartingPosition = sequenceNumber == 0 ? 0 : fastaVector->metadata.data[sequenceNumber - 1].sequenceEndPosition;
      sequenceLength = (fastaVector->metadata.data[sequenceNumber].sequenceEndPosition - sequenceStartingPosition);


      std::cout << "phmm #"<< phmmNumber <<", seq #"<< sequenceNumber<< ", debug: seq len " << sequenceLength << std::endl;

      for (uint32_t sequencePosition = 0; sequencePosition < sequenceLength; sequencePosition++) {
        char symbolFromSequence = fastaVector->sequence.charData[sequenceStartingPosition + sequencePosition];
        uint8_t symbolEncoding;
        switch (symbolFromSequence) {
        case 'a': case 'A': symbolEncoding = 0; break;
        case 'c': case 'C': symbolEncoding = 1; break;
        case 'g': case 'G': symbolEncoding = 2; break;
        default:            symbolEncoding = 3;
        }

        uint8_t previousCellValue = 0;
        for (uint32_t phmmPosition = 0; phmmPosition < phmmLength; phmmPosition++) {
          float matchScoreFromPhmmFile = phmmMatchEmissionScores[(phmmPosition * cardinality) + symbolEncoding];
          float projectedMatchScore = emissionScoreToProjectedScore(matchScoreFromPhmmFile, scoreMultiplier);
          int8_t projectedMatchScoreAsInt = round(projectedMatchScore);
          uint8_t tempCellValue = cellScores->data()[phmmPosition]; //store what was here before the summation so we can give it to the next cell.
          int16_t cellSumResults = (int16_t)previousCellValue + (int16_t)projectedMatchScoreAsInt;
          previousCellValue = tempCellValue;


          if (cellSumResults >= 256) {
            ReferenceSsvHit hit = ReferenceSsvHit();
            hit.sequenceNumber = sequenceNumber;
            hit.phmmNumber = phmmNumber;
            hit.sequencePosition = sequencePosition;
            hit.phmmPosition = phmmPosition;
            hitsFromReference->push_back(hit);

            cellSumResults = 0;
          }
          if (cellSumResults <= 0) {
            cellSumResults = 0;
          }
          cellScores->data()[phmmPosition] = cellSumResults;
        }
      }
    }
  }

  return hitsFromReference;
}

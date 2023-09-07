#include "byCellComparator.hpp"
#include <map>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <string>
#include <sstream>

// #define PRINT_MATCHING_CELLS

std::map<CellCompareKey, CellCompareValue> hardwareMap;
std::map<CellCompareKey, CellCompareValue> softwareMap;


std::string cellKeyToString(struct CellCompareKey key) {
  std::ostringstream stringStream;
  stringStream << "[s: " << key.globalSequenceIndex << ", pI: " << key.phmmIndex << "]";
  return stringStream.str();
}

std::string cellValueToString(struct CellCompareValue value) {
  std::ostringstream stringStream;
  stringStream << "prev: " << (int)value.prevValue << "\tmatch: " << (int)value.matchScore << "\tscore" <<
    (int)value.cellValue << "\thit?" << (int)value.passesThreshold << "\tsym:" << (int)value.symbol <<" P[" <<
    (int)value.phmmVector[0] << ", " << (int)value.phmmVector[1] << ", " << (int)value.phmmVector[2] << ", " << (int)value.phmmVector[3] << "]";
  return stringStream.str();
}

bool compareCellValues(struct CellCompareValue v1, struct CellCompareValue v2) {
  return (v1.cellValue == v2.cellValue && v1.matchScore == v2.matchScore &&
    v1.passesThreshold == v2.passesThreshold && v1.prevValue == v2.prevValue &&
    v1.phmmVector[0] == v2.phmmVector[0] && v1.phmmVector[1] == v2.phmmVector[1] &&
    v1.phmmVector[2] == v2.phmmVector[2] && v1.phmmVector[3] == v2.phmmVector[3] &&
    v1.symbol == v2.symbol);
}


void addKvToHardwareMap(struct CellCompareKey key, struct CellCompareValue value) {
  hardwareMap.emplace(key, value);
}

void addKvToSoftwareMap(struct CellCompareKey key, struct CellCompareValue value) {
  softwareMap.emplace(key, value);
}


bool compareCellMaps(uint32_t phmmLength, uint32_t sequenceLength, struct HmmSeqPair hmmSeqPair) {
  bool allCellsSuccessful = true;
  for (uint32_t sequenceIndex = 0; sequenceIndex < sequenceLength; sequenceIndex++) {
    for (uint32_t phmmIndex = 0; phmmIndex < phmmLength; phmmIndex++) {
      struct CellCompareKey key = { .phmmIndex = phmmIndex, .globalSequenceIndex = sequenceIndex };


      bool cellRepresentedInBothMaps = true;
      //check to make sure key is in both maps
      auto hwSearch = hardwareMap.find(key);
      if (hwSearch == hardwareMap.end()) {
        //cell not found in hardware map
        std::cout << "TEST FAIL: cell for " << cellKeyToString(key) << " was not represented in the hardwareMap\n";
        bool cellRepresentedInBothMaps = false;
      }

      auto swSearch = softwareMap.find(key);
      if (swSearch == softwareMap.end()) {
        //cell not found in software map
        std::cout << "TEST FAIL: cell for " << cellKeyToString(key) << "was not represented in the softwareMap\n";
        bool cellRepresentedInBothMaps = false;
      }

      if (!cellRepresentedInBothMaps) {
        allCellsSuccessful = false;
      }
      else {
        struct CellCompareValue hardwareValue = hwSearch->second;
        struct CellCompareValue softwareValue = swSearch->second;

        if (!compareCellValues(hardwareValue, softwareValue)) {
          std::cout << "TEST FAIL: cell " << cellKeyToString(hwSearch->first) << "\n\tHW: " << cellValueToString(hardwareValue) <<
            "\n\tSW: " << cellValueToString(softwareValue) << "\n";

          int8_t* phmmVectorPtr = (int8_t*)&hmmSeqPair.phmmList->phmms[0].model.matchEmissionScores[hwSearch->first.phmmIndex];
            std::cout << "actual symbol: " << hmmSeqPair.sequence[hwSearch->first.globalSequenceIndex] << "\tactual phmmVec: ["<<
              (int)(signed char)(phmmVectorPtr[0]) << ", " << (int)(signed char)(phmmVectorPtr[1]) << ", " << (int)(signed char)(phmmVectorPtr[2]) << ", " << (int)(signed char)(phmmVectorPtr[3]) << "]\n";
          allCellsSuccessful = false;
        }
        else {
          #ifdef PRINT_MATCHING_CELLS
          std::cout << "cell match " << cellKeyToString(hwSearch->first) << "->" << cellValueToString(hardwareValue, softwareValue) << "\n";
          #endif
        }

      }
    }
  }
  return allCellsSuccessful;
}


void clearCellComparatorMaps() {
  hardwareMap.clear();
  softwareMap.clear();
}


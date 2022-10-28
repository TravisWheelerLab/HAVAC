#pragma once
#include <cstdint>
#include "../test/generator/hmmSeqGenerator.h"

struct CellCompareKey {
  uint32_t phmmIndex;
  uint32_t globalSequenceIndex;

  bool operator==(const CellCompareKey& k2) const {
    return globalSequenceIndex == k2.globalSequenceIndex && phmmIndex == k2.phmmIndex;
  }

  bool operator<(const CellCompareKey& k2)  const {
    if (globalSequenceIndex == k2.globalSequenceIndex) {
      return phmmIndex < k2.phmmIndex;
    }
    else {
      return globalSequenceIndex > k2.globalSequenceIndex;
    }
  }
};

struct CellCompareValue {
  uint8_t prevValue;
  int8_t matchScore;
  uint8_t cellValue;
  int8_t phmmVector[4];
  uint8_t symbol;
  bool passesThreshold;


};

void addKvToHardwareMap(CellCompareKey key, CellCompareValue value);
void addKvToSoftwareMap(CellCompareKey key, CellCompareValue value);
bool compareCellMaps(uint32_t phmmLength, uint32_t sequenceLength, struct HmmSeqPair hmmSeqPair);
void clearCellComparatorMaps();
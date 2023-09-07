#pragma once
#include <FastaVector.h>
#include <p7HmmReader.h>
#include <cstdint>
#include <string>

class HavacHit {
public:
  HavacHit(const uint64_t sequencePosition, const uint32_t sequenceIndex, const uint32_t phmmPosition, const uint32_t phmmIndex);
  uint64_t sequencePosition;
  uint32_t sequenceIndex;
  uint32_t phmmPosition;
  uint32_t phmmIndex;

  std::string toString();
};

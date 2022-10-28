#pragma once

#include <cstdint>
extern "C" {
#include "../../P7HmmReader/src/p7HmmReader.h"
}

struct HmmSeqPair{
  char *sequence;
  uint32_t sequenceLength;
  P7HmmList *phmmList;
};



struct HmmSeqPair generateRandomHmmSeqPair(const uint32_t seqLength);

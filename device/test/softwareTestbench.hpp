#pragma once

#include <stdint.h>
#include <vector>
#include "../HavacHls.hpp"


uint8_t* generateCompressedSequence(const char* sequence);
std::vector<struct TestbenchHitReport> invokeHardwareSsv(const uint8_t* sequenceAsVectorIndices, size_t sequenceLength,
  uint32_t *emissionsAsUInt32_t, const size_t phmmLength, int8_t& errorCode);
bool compareSsvHitLists(std::vector<struct TestbenchHitReport> hardwareSsvHits, std::vector<struct SoftSsvHit> softwareSsvHits);
struct HmmSeqPair readTestSetFromFiles(const uint32_t testNum);
uint8_t* sequenceAsUnpackedVectorIndices(struct HmmSeqPair& hmmSeqPair);

#pragma once

#include <stdint.h>
#include <vector>
#include "../device/HavacHls.hpp"


uint8_t* generateCompressedSequence(const char* sequence);
std::vector<struct HitReport> invokeHardwareSsv(const uint8_t* sequenceAsVectorIndices, size_t sequenceLength,
  uint32_t *emissionsAsUInt32_t, const size_t phmmLength, int8_t& errorCode);
bool compareSsvHitLists(std::vector<struct HitReport> hardwareSsvHits, std::vector<struct SoftSsvHit> softwareSsvHits);

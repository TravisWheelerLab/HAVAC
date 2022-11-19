#pragma once

#include <stdint.h>
#include <vector>
#include "../device/HavacHls.hpp"


uint8_t* generateCompressedSequence(const char* sequence);
#ifdef USE_HIT_SIEVE
inline std::vector<struct HitReportByGroup> invokeHardwareSsv(const uint8_t* sequenceAsVectorIndices, size_t sequenceLength,
  const int8_t* emissionsAsInt8_t, const size_t phmmLength, int8_t& errorCode);
#else
inline std::vector<struct HitReport> invokeHardwareSsv(const uint8_t* sequenceAsVectorIndices, size_t sequenceLength,
  const int8_t* emissionsAsInt8_t, const size_t phmmLength, int8_t& errorCode);
#endif


#ifdef USE_HIT_SIEVE
bool compareSsvHitLists(std::vector<struct HitReportByGroup> hardwareSsvHits, std::vector<struct SoftSsvHit> softwareSsvHits);
#else
	bool compareSsvHitLists(std::vector<struct HitReport> hardwareSsvHits, std::vector<struct SoftSsvHit> softwareSsvHits);
#endif
bool compareSsvHitLists(std::vector<struct HitReport> hardwareSsvHits, std::vector<struct SoftSsvHit> softwareSsvHits);

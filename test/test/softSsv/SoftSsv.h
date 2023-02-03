#ifndef HAVAC_SOFT_SSV_H
#define HAVAC_SOFT_SSV_H

#include <cstdint>
#include <cstdlib>
#include <vector>

struct SoftSsvHit{
	uint32_t sequencePosition;
	uint32_t phmmPosition;
};

//will return -1 for both positions if not found.
std::vector<struct SoftSsvHit> softSsvThreshold256(const uint8_t* sequenceAsVectorIndices, const uint64_t sequenceLength,
	const int8_t* flattenedPhmmEmissions, size_t phmmLengthInVectors, int8_t& errorCode);


#endif

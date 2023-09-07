#include "SoftSsv.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "../../device/PublicDefines.h"

//if HAVAC_PER_CELL_DATA_TESTING is defined, it's defined in PublicDefines.h
#ifdef HAVAC_PER_CELL_DATA_TESTING
#include "../byCellComparator/byCellComparator.hpp"
#endif


//perform a reference software SSV implementation to compare the hardware results to.
//if HAVAC_PER_CELL_DATA_TESTING is defined, all cells will record their data in the softwareSsv cell hashtable
std::vector<struct SoftSsvHit> softSsvThreshold256(const uint8_t* sequenceAsVectorIndices, const uint64_t sequenceLength,
	const int8_t* flattenedPhmmEmissions, size_t phmmLengthInVectors, int8_t& errorCode) {
	printf("beginning soft ssv\n");
	fflush(stdout);
	errorCode = 0;
	std::vector<struct SoftSsvHit> hitReportList;

	struct SoftSsvHit ssvHitLocation;
	ssvHitLocation.phmmPosition = -1;
	ssvHitLocation.sequencePosition = -1;

	uint8_t* processingRowScores = (uint8_t*)calloc(sequenceLength, sizeof(uint8_t));
	if (processingRowScores == NULL) {
		printf("ERROR: failed to allocate processing row scores in softSsvThreshold256\n");
		fflush(stdout);
		errorCode = -1;
		return hitReportList;
	}
	printf("allocated procrowsc\n");
	fflush(stdout);

	for (uint32_t phmmPosition = 0; phmmPosition < phmmLengthInVectors; phmmPosition++) {
		for (int64_t sequencePosition = sequenceLength - 1; sequencePosition >= 0; sequencePosition--) {

			const int8_t* phmmVector = &flattenedPhmmEmissions[phmmPosition * 4];
			const uint8_t sequenceSymbol = sequenceAsVectorIndices[sequencePosition];
			const int8_t phmmScore = phmmVector[sequenceSymbol];
			const uint8_t prevScore = sequencePosition == 0 ? 0 : processingRowScores[sequencePosition - 1];
			int32_t cellScore = (int32_t)phmmScore + (int32_t)prevScore;
			bool passesThreshold = false;

			if (cellScore < 0) {
				cellScore = 0;
			}
			else if (cellScore >= 256) {
				fflush(stdout);
				cellScore = 0;
				passesThreshold = true;
				ssvHitLocation.phmmPosition = phmmPosition;
				ssvHitLocation.sequencePosition = sequencePosition;
				hitReportList.push_back(ssvHitLocation);
				fflush(stdout);
			}
			// if we're performing per-cell testing, record the cell's input and output data
			#ifdef HAVAC_PER_CELL_DATA_TESTING
			struct CellCompareKey key = { .phmmIndex = phmmPosition, .globalSequenceIndex = (uint32_t)sequencePosition };
			struct CellCompareValue value = { .prevValue = prevScore, .matchScore = phmmScore,
			.cellValue = (uint8_t)cellScore, .phmmVector = {phmmVector[0], phmmVector[1], phmmVector[2], phmmVector[3]},
			.symbol = sequenceSymbol, .passesThreshold = passesThreshold };
			addKvToSoftwareMap(key, value);
			#endif

			processingRowScores[sequencePosition] = cellScore;

		}
	}

	free(processingRowScores);
	return hitReportList;
}

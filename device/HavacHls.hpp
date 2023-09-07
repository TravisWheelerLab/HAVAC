#ifndef HAVAC_HLS_HPP
#define HAVAC_HLS_HPP


#include <stdint.h>
#include <hls_stream.h>
#include <iostream>
#include <ap_int.h>
#include "PublicDefines.h"
#include "HitReporting.hpp"
#include "Sequence.hpp"

#define PHMM_STREAM_DEPTH 64
#define HIT_REPORT_STREAM_DEPTH 64
#define NUM_HITS_STREAM_DEPTH 2
#define SCORES_PER_PHMM_VECTOR 4


class TestbenchHitReport {
public:
	TestbenchHitReport(const uint64_t reportAsU64) {
		uint64_t sequencePartitionBitMask = (1ULL << 14) - 1;
		uint64_t sequenceSegmentBitmask = (1ULL << 40) - 1;
		uint64_t sequencePartitionIndex = reportAsU64 & sequencePartitionBitMask;
		uint64_t sequenceSegmentIndex = (reportAsU64 & sequenceSegmentBitmask) >> 14;
		this->sequencePosition = (reportAsU64 & sequenceSegmentBitmask) >> 14;
		this->sequencePosition = (sequenceSegmentIndex * (12 * 1024)) + sequencePartitionIndex;

		this->phmmPosition = reportAsU64 >> 40;
	}
	uint32_t phmmPosition;
	uint32_t sequencePosition;
};

struct PhmmVector {
	ap_uint<8> scores[SCORES_PER_PHMM_VECTOR];
};

struct CellResult {
	ap_uint<8> cellScore;
	bool passesThreshold;
};


struct MatchScoreList {
	ap_uint<8> scores[NUM_CELL_PROCESSORS];
};

// function prototypes
void HavacKernel(SequenceSegmentWord* sequenceSegmentMemory, uint32_t sequenceLengthInSegments,
	uint32_t* phmmMemory, uint32_t phmmLengthInVectors, uint64_t* hitReportMemory,
	uint32_t* hitReportCountMemory);

void HavacTopLevelDataflow(const seqSegPos_t sequenceLengthInSegments, const phmmPos_t phmmLengthInVectors,
	SequenceSegmentWord* sequenceSegmentMemory, uint32_t* phmmVectorMemory,
	uint64_t* hitReportMemory, uint32_t* hitReportCountMemory);

void HavacMainLoop(const seqSegPos_t sequenceLengthInSegments, const phmmPos_t phmmLengthInVectors,
	SequenceSegmentWord* sequenceSegmentMemory, uint32_t* phmmVectorMemory,
	hls::stream<PositionReport, inputHitReportStreamDepth>& inputPositionReportStream,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_0,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_1,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_2,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_3,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_4,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_5,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_6,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_7,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_8,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_9,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_10,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_11,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_12,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_13,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_14,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_15);

void HavacDataflowFunction(const seqSegPos_t sequenceLengthInSegments, const phmmPos_t phmmLengthInVectors,
	SequenceSegmentWord* sequenceSegmentMemory, uint32_t* phmmVectorMemory,
	seqSegPos_t sequenceSegmentIndex,
	hls::stream<PositionReport, inputHitReportStreamDepth>& inputPositionReportStream,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_0,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_1,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_2,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_3,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_4,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_5,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_6,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_7,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_8,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_9,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_10,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_11,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_12,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_13,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_14,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_15);

void phmmVectorLoop(phmmPos_t phmmLengthInVectors, hls::stream<SequenceSegmentWord, SEQUENCE_STREAM_DEPTH>& sequenceSegmentStream,
	hls::stream<uint32_t, PHMM_STREAM_DEPTH>& phmmStream, bool isFirstSequenceSegment, bool isLastSequenceSegment,
	seqSegPos_t sequenceSegmentIndex,
	hls::stream<PositionReport, inputHitReportStreamDepth>& inputPositionReportStream,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_0,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_1,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_2,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_3,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_4,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_5,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_6,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_7,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_8,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_9,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_10,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_11,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_12,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_13,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_14,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_15);

void computeAllCellProcessors(ap_uint<8> cellScores[NUM_CELL_PROCESSORS], struct MatchScoreList& matchScoreList,
	ap_uint<CELLS_PER_GROUP> cellsPassingThreshold[NUM_CELL_GROUPS], ap_uint<8> leftScoreIn,
	bool isFirstPhmmIndex, ap_uint<8>& lastScoreOut);

// compute the next value of the cell, and determine if it passes the implicit 256 threshold
struct CellResult computeCellProcessor(ap_uint<8> prevScore, ap_uint<8> matchScore);
ap_uint<8> readScoreFromScoreQueue(bool isFirstSequenceSegment, bool isLastPhmmIndex);
void writeScoreToScoreQueue(ap_uint<8> scoreToWrite, bool isLastSequenceSegment, bool isFirstPhmmIndex);
void generateMatchScoreList(struct MatchScoreList& matchScoreList, uint32_t currentPhmmVector, struct SequenceSegment currentSequenceSegment);
void isFirstOrLastPhmmVector(phmmPos_t phmmVectorIndex, phmmPos_t phmmLengthInVectors, bool& isFirstPhmmVector, bool& isLastPhmmVector);
void loadPhmmStream(uint32_t* phmmVectorMemory, hls::stream<uint32_t, PHMM_STREAM_DEPTH>& phmmVectorStream, const phmmPos_t phmmLengthInVectors);
void setPhmmFromStream(uint32_t& phmmVector, hls::stream<uint32_t, PHMM_STREAM_DEPTH>& phmmVectorStream);

#endif

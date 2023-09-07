#ifndef HAVAC_HLS_HPP
#define HAVAC_HLS_HPP


#include <stdint.h>
#include <hls_stream.h>
#include <iostream>
#include <ap_int.h>
#include "PublicDefines.h"
#include "ScoreQueue.hpp"
#include "HitReporting.hpp"
#include "Sequence.hpp"

#define PHMM_STREAM_DEPTH 64
#define HIT_REPORT_STREAM_DEPTH 64
#define NUM_HITS_STREAM_DEPTH 2


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
  uint32_t* phmmMemory, uint32_t phmmLengthInVectors, struct HitReport* hitReportMemory,
  uint32_t *hitReportCountMemory);


void HavacMainLoop(const uint32_t sequenceLengthInSegments, const uint32_t phmmLengthInVectors,
  SequenceSegmentWord* sequenceSegmentMemory, uint32_t* phmmVectorMemory, struct HitReport* hitReportMemory,
  hls::stream<uint32_t, NUM_HITS_STREAM_DEPTH>& numHitsStream, ScoreQueue &scoreQueue);


void HavacDataflowFunction(const uint32_t sequenceLengthInSegments, const uint32_t phmmLengthInVectors,
  SequenceSegmentWord* sequenceSegmentMemory, uint32_t* phmmVectorMemory, struct HitReport* hitReportMemory,
  hls::stream<uint32_t, NUM_HITS_STREAM_DEPTH>& numHitsStream, uint32_t sequenceSegmentIndex,
  ScoreQueue &scoreQueue);

void phmmVectorLoop(uint32_t phmmLengthInVectors, hls::stream<SequenceSegmentWord, SEQUENCE_STREAM_DEPTH>& sequenceSegmentStream,
		hls::stream<uint32_t, PHMM_STREAM_DEPTH>& phmmStream, bool isFirstSequenceSegment, bool isLastSequenceSegment,
		hls::stream<struct HitReportWithTerminator, HIT_REPORT_STREAM_DEPTH>& hitReportStream, uint32_t sequenceSegmentIndex,
		ScoreQueue &scoreQueue);

void loadSequenceSegmentStream(SequenceSegmentWord* sequenceSegmentMemory, hls::stream<SequenceSegmentWord, SEQUENCE_STREAM_DEPTH>& sequenceSegmentStream,
  uint32_t sequenceSegmentIndex);

// Read Data from Global Memory and write into Stream inStream

void computeAllCellProcessors(ap_uint<8> cellScores[NUM_CELL_PROCESSORS], struct MatchScoreList& matchScoreList,
  ap_uint<CELLS_PER_GROUP> cellsPassingThreshold[NUM_CELL_GROUPS], ap_uint<8> leftScoreIn,
  bool isFirstPhmmIndex, ap_uint<8>& lastScoreOut);

// compute the next value of the cell, and determine if it passes the implicit 256 threshold
struct CellResult computeCellProcessor(ap_uint<8> prevScore, ap_uint<8> matchScore);

// takes the phmm/sequence index, and a threshold pass bit from every cell processor, compresses
// to a per-group basis, and enqueues a hit report if any of the groups had a hit in them.
void adjudicateAndWriteHitReport(hls::stream<struct HitReportWithTerminator, HIT_REPORT_STREAM_DEPTH>& hitReportStream,
  ap_uint<CELLS_PER_GROUP> cellsPassingThreshold[NUM_CELL_GROUPS], uint32_t phmmIndex, uint32_t sequenceIndex, bool isFinalReportOfSegment);

void writeHitsToMemory(struct HitReport* hitReportMemory, hls::stream<struct HitReportWithTerminator, HIT_REPORT_STREAM_DEPTH>& hitReportStream,
		hls::stream<uint32_t, NUM_HITS_STREAM_DEPTH>& numHitsStream, bool isFirstSequenceSegment, bool isLastSequenceSegment);

//void writeScoreQueueFifo(hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE>& scoreQueueFifo, hls::stream<ap_uint<8>,
//  SCORE_QUEUE_STREAM_DEPTH> &scoreQueueWriteStream, uint32_t phmmLengthInVectors, bool isLastSequenceSegment);
//
//void readScoreQueueFifo(hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE>& scoreQueueFifo, hls::stream<ap_uint<8>,
//  SCORE_QUEUE_STREAM_DEPTH>& scoreQueueReadStream, uint32_t phmmLengthInVectors, bool isFirstSequenceSegment);

//ap_uint<8> readScoreFromScoreQueue(hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE>& scoreQueue, bool isFirstSequenceSegment, bool isLastPhmmIndex);
ap_uint<8> readScoreFromScoreQueue(ScoreQueue &scoreQueue, bool isFirstSequenceSegment, bool isLastPhmmIndex);
//void writeScoreToScoreQueue(hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE>& scoreQueue, ap_uint<8>& scoreToWrite, bool isLastSequenceSegment, bool isFirstPhmmIndex);
void writeScoreToScoreQueue(ScoreQueue &scoreQueue, ap_uint<8>& scoreToWrite, bool isLastSequenceSegment, bool isFirstPhmmIndex);
void generateMatchScoreList(struct MatchScoreList &matchScoreList, uint32_t currentPhmmVector, struct SequenceSegment currentSequenceSegment);
void isFirstOrLastSequenceSegment(uint32_t sequenceSegmentIndex, uint32_t sequenceLengthInSegments, bool& isFirstSequenceSegment, bool& isLastSequenceSegment);
void isFirstOrLastPhmmVector(uint32_t phmmVectorIndex, uint32_t phmmLengthInVectors, bool& isFirstPhmmVector, bool& isLastPhmmVector);
void loadPhmmStream(uint32_t* phmmVectorMemory, hls::stream<uint32_t, PHMM_STREAM_DEPTH>& phmmVectorStream, const uint32_t phmmLengthInVectors);
//void setSequenceSegment(struct SequenceSegment& currentSequenceSegment, hls::stream<SequenceSegmentWord, SEQUENCE_STREAM_DEPTH>& sequenceSegmentStream);
void setPhmmFromStream(uint32_t& phmmVector, hls::stream<uint32_t, PHMM_STREAM_DEPTH>& phmmVectorStream);
void copyScalarInputs(uint32_t sequenceSegmentLengthInSegments, uint32_t &localSequenceSegmentLengthInSegments,
	uint32_t phmmLengthInVectors, uint32_t &localPhmmLengthInVectors);

#endif

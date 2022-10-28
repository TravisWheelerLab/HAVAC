#ifndef HAVAC_HLS_HPP
#define HAVAC_HLS_HPP


#include <stdint.h>
#include <hls_stream.h>
#include <iostream>
#include <ap_int.h>
#include "PublicDefines.h"

#define SEQUENCE_STREAM_DEPTH 2
#define PHMM_STREAM_DEPTH 512        // 1 ultraram is 288Kb, or 9,000 phmm vectors
#define HIT_REPORT_STREAM_DEPTH 64    // 1 bram is 36Kb, or 4,500 bytes. each hit report is 16 bytes
#define ULTRARAM_BYTES_PER_RAM 36000
#define SCORE_QUEUE_SIZE (1200* ULTRARAM_BYTES_PER_RAM) // 1 ultraram is 288Kb, or  36,000 bytes.

//HAVAC_TESTING is defined in PublicDefines.h, if it's defined at all
#ifdef HAVAC_TESTING
#define SYMBOLS_PER_SEQUENCE_SEGMENT_WORD 1024
#else
#define SYMBOLS_PER_SEQUENCE_SEGMENT_WORD 2048
#endif

#define NUM_SEQUENCE_SEGMENT_WORDS (NUM_CELL_PROCESSORS/SYMBOLS_PER_SEQUENCE_SEGMENT_WORD)

// data structures
//sequence segment is implemented with a wide scalar instead of array of 2-bit values because
//the array implementation blew up utilization like CRAZY just loading the segment.
//once num cell processors is larger than 32*64, this will need to be an array of wide scalars (likely 2 values)
struct SequenceSegment {
  //  ap_uint<NUM_CELL_PROCESSORS * 2> symbols;
  ap_uint<SYMBOLS_PER_SEQUENCE_SEGMENT_WORD * 2> symbols[NUM_SEQUENCE_SEGMENT_WORDS];
};
struct PhmmVector {
  ap_uint<8> scores[SCORES_PER_PHMM_VECTOR];
};
struct HitReport {
  uint32_t phmmIndex;
  uint32_t sequenceIndex;
  ap_uint<NUM_CELL_GROUPS> groupsPassingThreshold;
};

struct CellResult {
  ap_uint<8> cellScore;
  bool passesThreshold;
};


struct MatchScoreList {
  ap_uint<8> scores[NUM_CELL_PROCESSORS];
};

// function prototypes
void HavacKernelTopLevel(struct SequenceSegment* sequenceSegmentMemory, uint32_t sequenceLengthInSegments,
  struct PhmmVector* phmmMemory, uint32_t phmmLengthInVectors, struct HitReport* hitReportMemory, uint32_t& numHits);


void HavacMainLoop(uint32_t sequenceLengthInSegments, hls::stream<struct SequenceSegment, SEQUENCE_STREAM_DEPTH>& sequenceSegmentStream,
  hls::stream<struct PhmmVector, PHMM_STREAM_DEPTH>& phmmStream, uint32_t phmmLengthInVectors,
  hls::stream<struct HitReport, HIT_REPORT_STREAM_DEPTH>& hitReportStream, struct SequenceSegment* sequenceSegmentMemory,
  struct PhmmVector* phmmVectorMemory, struct HitReport* hitReportMemory);


void phmmVectorLoop(uint32_t phmmLengthInVectors, struct SequenceSegment currentSequenceSegment, hls::stream<struct PhmmVector,
  PHMM_STREAM_DEPTH>& phmmStream, bool isFirstSequenceSegment, bool isLastSequenceSegment, hls::stream<struct HitReport,
  HIT_REPORT_STREAM_DEPTH>& hitReportStream, uint32_t sequenceSegmentIndex);

void loadSequenceSegmentStream(struct SequenceSegment* sequenceSegmentMemory, hls::stream<struct SequenceSegment>& sequenceSegmentStream,
  uint32_t sequenceSegmentIndex);

// Read Data from Global Memory and write into Stream inStream

void computeAllCellProcessors(ap_uint<8> cellScores[NUM_CELL_PROCESSORS], struct MatchScoreList& matchScoreList,
  ap_uint<CELLS_PER_GROUP> cellsPassingThreshold[NUM_CELL_GROUPS], ap_uint<8> leftScoreIn,
  bool isFirstPhmmIndex, ap_uint<8>& lastScoreOut);

// compute the next value of the cell, and determine if it passes the implicit 256 threshold
struct CellResult computeCellProcessor(ap_uint<8> prevScore, ap_uint<8> matchScore);

// takes the phmm/sequence index, and a threshold pass bit from every cell processor, compresses
// to a per-group basis, and enqueues a hit report if any of the groups had a hit in them.
void adjudicateAndWriteHitReport(hls::stream<struct HitReport>& hitReportStream,
  ap_uint<CELLS_PER_GROUP> cellsPassingThreshold[NUM_CELL_GROUPS], uint32_t phmmIndex, uint32_t sequenceIndex);

void writeHitsToMemory(struct HitReport* hitReportMemory, hls::stream<struct HitReport, HIT_REPORT_STREAM_DEPTH>& hitReportStream);

//void writeScoreQueueFifo(hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE>& scoreQueueFifo, hls::stream<ap_uint<8>,
//  SCORE_QUEUE_STREAM_DEPTH> &scoreQueueWriteStream, uint32_t phmmLengthInVectors, bool isLastSequenceSegment);
//
//void readScoreQueueFifo(hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE>& scoreQueueFifo, hls::stream<ap_uint<8>,
//  SCORE_QUEUE_STREAM_DEPTH>& scoreQueueReadStream, uint32_t phmmLengthInVectors, bool isFirstSequenceSegment);

ap_uint<8> readScoreFromScoreQueue(hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE>& scoreQueue, bool isFirstSequenceSegment, bool isLastPhmmIndex);
void writeScoreToScoreQueue(hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE>& scoreQueue, ap_uint<8>& scoreToWrite, bool isLastSequenceSegment, bool isFirstPhmmIndex);
void generateMatchScoreList(struct MatchScoreList& matchScoreList, struct PhmmVector& currentPhmmVector, struct SequenceSegment& currentSequenceSegment);
void isFirstOrLastSequenceSegment(uint32_t sequenceSegmentIndex, uint32_t sequenceLengthInSegments, bool& isFirstSequenceSegment, bool& isLastSequenceSegment);
void isFirstOrLastPhmmVector(uint32_t phmmVectorIndex, uint32_t phmmLengthInVectors, bool& isFirstPhmmVector, bool& isLastPhmmVector);
void loadPhmmStream(struct PhmmVector* phmmVectorMemory, hls::stream<struct PhmmVector, PHMM_STREAM_DEPTH>& phmmVectorStream, const uint32_t phmmLengthInVectors);
void setSequenceSegment(struct SequenceSegment& currentSequenceSegment, hls::stream<struct SequenceSegment, SEQUENCE_STREAM_DEPTH>& sequenceSegmentStream);
void setPhmmFromStream(struct PhmmVector& phmmVector, hls::stream<struct PhmmVector, PHMM_STREAM_DEPTH>& phmmVectorStream);

#endif

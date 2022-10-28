#define AP_INT_MAX_W 4096
#include "HavacHls.hpp"

#include "hls_print.h"
#define AP_INT_MAX_W 4096

#ifdef HAVAC_PER_CELL_DATA_TESTING
#include "../byCellComparator/byCellComparator.hpp"
#endif

#ifdef HAVAC_PER_CELL_DATA_TESTING
static uint32_t _phmmIndex;
static uint32_t _sequenceSegmentIndex;
static uint32_t _cellIndexInSegment;
static int8_t _phmmVector[4];
static uint8_t _symbols[NUM_CELL_PROCESSORS];
#endif
static uint32_t numHitReportsSent;

//top level function that is invoked to run HAVAC. This project uses Vitis HLS, please read the documentation for this tech.
//https://docs.xilinx.com/r/en-US/ug1399-vitis-hls/Getting-Started-with-Vitis-HLS
void HavacKernelTopLevel(struct SequenceSegment* sequenceSegmentMemory, uint32_t sequenceLengthInSegments,
  struct PhmmVector* phmmVectorMemory, uint32_t phmmLengthInVectors, struct HitReport* hitReportMemory, uint32_t &numHits) {

  // note: bundle names should be all lowercase
  // depth is just a simulation suggestion, for largest fifo needed for co-simulation
  //seq and phmm lengths are scalars, so they are set via the axilite interface
  //the pointer values are memory partitions in RAM, and use AXI interfaces behind the scenes
  //https://docs.xilinx.com/r/en-US/ug1399-vitis-hls/pragma-HLS-interface
  #pragma HLS INTERFACE PORT = sequenceLengthInSegments mode = s_axilite register
  #pragma HLS INTERFACE PORT = phmmLengthInVectors mode = s_axilite register
  #pragma HLS INTERFACE mode = m_axi port = sequenceSegmentMemory bundle = gmem0 depth =256 num_read_outstanding = 8
  #pragma HLS INTERFACE mode = m_axi port = phmmVectorMemory bundle = gmem1 depth = 356 num_read_outstanding = 8
  #pragma HLS INTERFACE mode = m_axi port = hitReportMemory bundle = gmem0 depth = 65535 num_write_outstanding = 8

  //aggregate the memory so we can access larger portions of the data each cycle
  //https://docs.xilinx.com/r/en-US/ug1399-vitis-hls/pragma-HLS-aggregate
  #pragma HLS aggregate variable = phmmVectorMemory compact = auto
  #pragma HLS aggregate variable = hitReportMemory compact = auto


  //streams to act as producer/consumer models to more effeciently read and write data to/from RAM
  static hls::stream<struct SequenceSegment, SEQUENCE_STREAM_DEPTH> sequenceSegmentStream("sequenceStream");
  static hls::stream<struct PhmmVector, PHMM_STREAM_DEPTH> phmmStream("phmmStream");
  static hls::stream<struct HitReport, HIT_REPORT_STREAM_DEPTH> hitReportStream("hitReportStream");

  //segment stream uses a 2-deep ping-pong buffer so that one segment can be filled while the
  //other one is being used. 
  #pragma HLS STREAM variable= sequenceSegmentStream type=pipo depth=2
  #pragma HLS array_partition variable=sequenceSegmentStream type=complete

  numHitReportsSent = 0;

  HavacMainLoop(sequenceLengthInSegments, sequenceSegmentStream, phmmStream, phmmLengthInVectors,
    hitReportStream, sequenceSegmentMemory, phmmVectorMemory, hitReportMemory);

  numHits = numHitReportsSent;
}


//main computation loop for HAVAC, each iteration computes a wide column down the Dynamic Programming matrix.
void HavacMainLoop(uint32_t sequenceLengthInSegments, hls::stream<struct SequenceSegment, SEQUENCE_STREAM_DEPTH>& sequenceSegmentStream,
  hls::stream<struct PhmmVector, PHMM_STREAM_DEPTH>& phmmStream, uint32_t phmmLengthInVectors,
  hls::stream<struct HitReport, HIT_REPORT_STREAM_DEPTH>& hitReportStream, struct SequenceSegment* sequenceSegmentMemory,
  struct PhmmVector* phmmVectorMemory, struct HitReport* hitReportMemory) {

HavacSequenceSegmentLoop:
  for (uint32_t sequenceSegmentIndex = 0; sequenceSegmentIndex < sequenceLengthInSegments; sequenceSegmentIndex++) {


    #ifdef HAVAC_PER_CELL_DATA_TESTING
    _sequenceSegmentIndex = sequenceSegmentIndex;
    #endif

    //https://docs.xilinx.com/r/en-US/ug1399-vitis-hls/pragma-HLS-dataflow
    #pragma HLS DATAFLOW
    struct SequenceSegment currentSequenceSegment;
    //partition the sequence segment so we can access every character at the same time, 
    // as opposed to only one character per cycle like if we didn't partition it.
    #pragma HLS array_partition variable=currentSequenceSegment type=complete

    //determine if this the first or last sequence segment
    bool isFirstSequenceSegment, isLastSequenceSegment;
    isFirstOrLastSequenceSegment(sequenceSegmentIndex, sequenceLengthInSegments, isFirstSequenceSegment, isLastSequenceSegment);
    loadSequenceSegmentStream(sequenceSegmentMemory, sequenceSegmentStream, sequenceSegmentIndex);
    loadPhmmStream(phmmVectorMemory, phmmStream, phmmLengthInVectors);
    setSequenceSegment(currentSequenceSegment, sequenceSegmentStream);
    phmmVectorLoop(phmmLengthInVectors, currentSequenceSegment, phmmStream, isFirstSequenceSegment,
      isLastSequenceSegment, hitReportStream, sequenceSegmentIndex);

    //write all hits this sequenceSegment pass
    writeHitsToMemory(hitReportMemory, hitReportStream);
  }
}

//computes a full column down the DP matrix, going through the whole phmm.
void phmmVectorLoop(uint32_t phmmLengthInVectors, struct SequenceSegment currentSequenceSegment, hls::stream<struct PhmmVector,
  PHMM_STREAM_DEPTH>& phmmStream, bool isFirstSequenceSegment, bool isLastSequenceSegment, hls::stream<struct HitReport,
  HIT_REPORT_STREAM_DEPTH>& hitReportStream, uint32_t sequenceSegmentIndex) {

  //hardware queue to store the scores from the right side of the dp matrix column, to be given as 
  // prev scores to the leftmost cells for the next column
  static hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE> scoreQueue;
  static ap_uint<8> cellScores[NUM_CELL_PROCESSORS];
  //scoreQueue should be stored in an UltraRAM
  #pragma HLS bind_storage variable = scoreQueue impl=uram type=fifo
  #pragma HLS array_partition variable = cellScores complete

  //we label this function as inline recursive so that all the called functions get inlined.
  // this helps timing because Vitis HLS can put pipeline stage breaks anywhere it wants in the computation pipeline,
  //not just between sub function calls
  #pragma HLS inline recursive

  //create a buffered register for both the score queue read and write.
  //the score queue read gets initialized to zero because the first cell in the column always gets a 0 (from outside the DP matrix)
  //the buffered write also helps with timing, and allows us to read and write in a seemly inverted way.
  //each phmmVector, we write whatever was buffered in bufferedScoreQueueWrite, then compute the cells, and 
  //only after we read what the next scoreQueueRead should be.
  ap_uint<8> bufferedScoreQueueRead(0);
  ap_uint<8> bufferedScoreQueueWrite;

HavacPhmmVectorLoop:
  for (uint32_t phmmIndex = 0; phmmIndex < phmmLengthInVectors; phmmIndex++) {

    #ifdef HAVAC_PER_CELL_DATA_TESTING
    _phmmIndex = phmmIndex;
    #endif

    //https://docs.xilinx.com/r/en-US/ug1399-vitis-hls/pragma-HLS-pipeline
    #pragma HLS PIPELINE II=1
    //create a match score list, i.e., for each cell, what score should they be given from the phmm vector?
    struct MatchScoreList matchScoreList;
    struct PhmmVector currentPhmmVector;
    // vector of bools where a 1 indicates the the corresponding cell passed its threshold this cycle.
    ap_uint<CELLS_PER_GROUP> cellsPassingThreshold[NUM_CELL_GROUPS];
    #pragma HLS array_partition type=complete variable=matchScoreList
    #pragma HLS array_partition type=complete variable=currentPhmmVector
    #pragma HLS array_partition variable=cellsPassingThreshold type=complete

    bool isFirstPhmmIndex, isLastPhmmIndex;
    isFirstOrLastPhmmVector(phmmIndex, phmmLengthInVectors, isFirstPhmmIndex, isLastPhmmIndex);

    setPhmmFromStream(currentPhmmVector, phmmStream);

    #ifdef HAVAC_PER_CELL_DATA_TESTING
    _phmmVector[0] = currentPhmmVector.scores[0];
    _phmmVector[1] = currentPhmmVector.scores[1];
    _phmmVector[2] = currentPhmmVector.scores[2];
    _phmmVector[3] = currentPhmmVector.scores[3];
    #endif

    generateMatchScoreList(matchScoreList, currentPhmmVector, currentSequenceSegment);

    writeScoreToScoreQueue(scoreQueue, bufferedScoreQueueWrite, isLastSequenceSegment, isFirstPhmmIndex);
    computeAllCellProcessors(cellScores, matchScoreList, cellsPassingThreshold, bufferedScoreQueueRead,
      isFirstPhmmIndex, bufferedScoreQueueWrite);
    bufferedScoreQueueRead = readScoreFromScoreQueue(scoreQueue, isFirstSequenceSegment, isLastPhmmIndex);
    //send hit data to adjudicate hit
    adjudicateAndWriteHitReport(hitReportStream, cellsPassingThreshold, phmmIndex, sequenceSegmentIndex);
  }
}


void computeAllCellProcessors(ap_uint<8> cellScores[NUM_CELL_PROCESSORS], struct MatchScoreList& matchScoreList,
		ap_uint<CELLS_PER_GROUP> cellsPassingThreshold[NUM_CELL_GROUPS], ap_uint<8> leftScoreIn, bool isFirstPhmmIndex, ap_uint<8>& lastScoreOut) {


  #ifdef HAVAC_PER_CELL_DATA_TESTING
  _cellIndexInSegment = NUM_CELL_PROCESSORS - 1;
  #endif


  //handle the final cell, since that'll need to save its score to lastScoreOut
  //this seems to need to be broken out to meet II=1. otherwise, vitis tries to schedule the second access to the final cell
  //on the next clock cycle
  ap_uint<8> thisCellPrevScore = isFirstPhmmIndex ? ap_uint<8>(0) : cellScores[NUM_CELL_PROCESSORS - 2];
  struct CellResult result = computeCellProcessor(thisCellPrevScore, matchScoreList.scores[NUM_CELL_PROCESSORS - 1]);
  lastScoreOut = result.cellScore;
  cellsPassingThreshold[NUM_CELL_GROUPS - 1].set_bit(CELLS_PER_GROUP- 1, result.passesThreshold);

  //iterate through the cell processors in reverse, since we need to use cell index i-1, and write to index i.
CellProcessorLoop:
  for (uint32_t cellIndex = NUM_CELL_PROCESSORS - 2; cellIndex != 0; cellIndex--) {

    #ifdef HAVAC_PER_CELL_DATA_TESTING
    _cellIndexInSegment = cellIndex;
    #endif
    #pragma HLS unroll
    ap_uint<8> thisCellPrevScore = isFirstPhmmIndex ? ap_uint<8>(0) : cellScores[cellIndex - 1];
    struct CellResult result = computeCellProcessor(thisCellPrevScore, matchScoreList.scores[cellIndex]);
    cellScores[cellIndex] = result.cellScore;
    cellsPassingThreshold[cellIndex / CELLS_PER_GROUP].set_bit(cellIndex % CELLS_PER_GROUP, result.passesThreshold);
  }

  #ifdef HAVAC_PER_CELL_DATA_TESTING
  _cellIndexInSegment = 0;
  #endif

  //handle the first cell processor, since it's an edge case that gets it's prev score from the score queue
  thisCellPrevScore = isFirstPhmmIndex ? ap_uint<8>(0) : leftScoreIn;
  result = computeCellProcessor(thisCellPrevScore, matchScoreList.scores[0]);
  cellScores[0] = result.cellScore;
  cellsPassingThreshold[0].set_bit(0, result.passesThreshold);
}


// Read Data from Global Memory and write into Stream inStream
void loadPhmmStream(struct PhmmVector* phmmVectorMemory, hls::stream<struct PhmmVector, PHMM_STREAM_DEPTH>& phmmVectorStream, const uint32_t phmmLengthInVectors) {
phmmReadInnerLoop:
  for (int phmmIndex = 0; phmmIndex < phmmLengthInVectors; phmmIndex++) {
    phmmVectorStream << phmmVectorMemory[phmmIndex];
  }
}


struct CellResult computeCellProcessor(ap_uint<8> prevScore, ap_uint<8> matchScore) {
  #pragma HLS INLINE
  //this function uses the overflow bit of an 8-bit addition to determine if a threshold hit happened,
  //or if the value needs to be reset to 0 because of a hit or an underflow to a negative number.
  //doing it this way saves a comparator inside the cell processes, dramatically cutting down on utilization.
  ap_uint<9> putativeSum = prevScore + matchScore;
  bool matchScoreSign = matchScore[7]; // matchScore is still only 8 bits
  bool sumCarryBit = putativeSum[8];   // putative sum got promoted to 9 bits
  // if there was a carry, but the matchScore was positive, we have an overflow, hence threshold hit.
  // if we had a negative match score but it didn't carry,  it's an underflow and we need to reset to 0.
  // if either are true, cellScore resets to 0.
  bool requiresReset = sumCarryBit != matchScoreSign;

  struct CellResult result;
  result.cellScore = requiresReset ? 0 : putativeSum.range(7, 0); // return 8 bits (without carry)
  result.passesThreshold = sumCarryBit && !matchScoreSign;

  #ifdef HAVAC_PER_CELL_DATA_TESTING
	  uint8_t _prevScore = prevScore.to_uint64();
	  int8_t _matchScore = matchScore.to_int64();
	  uint16_t _sum = putativeSum.to_uint64();

	  struct CellCompareKey key = { .phmmIndex = _phmmIndex,
	  .globalSequenceIndex = (_sequenceSegmentIndex * NUM_CELL_PROCESSORS) + _cellIndexInSegment };
	  struct CellCompareValue value = { .prevValue = _prevScore, .matchScore = _matchScore,
	  .cellValue = (uint8_t)result.cellScore.to_uint(), .phmmVector = {_phmmVector[0], _phmmVector[1], _phmmVector[2], _phmmVector[3]},
	  .symbol=_symbols[_cellIndexInSegment], .passesThreshold = result.passesThreshold};
	  addKvToHardwareMap(key, value);
  #endif


  return result;
}

//reads the 'cellsPassingThreshold' array of bools, groups them into hit groups, and if any of the cells in any hit group asserted,
//create a hit report and write it to the hit report stream.
void adjudicateAndWriteHitReport(hls::stream<struct HitReport>& hitReportStream,
  ap_uint<CELLS_PER_GROUP> cellsPassingThreshold[NUM_CELL_GROUPS], uint32_t phmmIndex, uint32_t sequenceIndex) {
	ap_uint<NUM_CELL_GROUPS> groupsPassingThreshold;
  for (uint32_t i = 0; i < NUM_CELL_GROUPS; i++) {
    #pragma HLS unroll
    groupsPassingThreshold.set_bit(i, cellsPassingThreshold[i].or_reduce());
  }

  if (groupsPassingThreshold.or_reduce()) {
    struct HitReport hitReport = { phmmIndex, sequenceIndex, groupsPassingThreshold };
    if (!hitReportStream.full()) {
      //perform a non-blocking write, so the write doesn't potentially stall the pipeline.
      //if we have so many hits to report that it blocks (unlikely, given size of stream), we can safely ignore the hit
      //because homology is almost certainly guaranteed if it's reporting that much.
      hitReportStream.write_nb(hitReport);
    }
    else {
      hls::print("hit report stream was full!\n");
    }
  }
}

//todo: this is causing warning HLS 200-1449 because of the while loop relying on the count from the prev thing in the dataflow
void writeHitsToMemory(struct HitReport* hitReportMemory, hls::stream<struct HitReport, HIT_REPORT_STREAM_DEPTH>& hitReportStream) {
  struct HitReport tempHitReport;
  while (!hitReportStream.empty()) {
    tempHitReport = hitReportStream.read();
    hitReportMemory[numHitReportsSent] = tempHitReport;
    numHitReportsSent++;
  }
}


void loadSequenceSegmentStream(struct SequenceSegment* sequenceSegmentMemory, hls::stream<struct SequenceSegment>& sequenceSegmentStream,
  uint32_t sequenceSegmentIndex) {
  sequenceSegmentStream << sequenceSegmentMemory[sequenceSegmentIndex];
}


//generates the list of match scores for all cells in the pipeline. 
//That is, for every cell processor, find the value from the phmm to add for that cell.
void generateMatchScoreList(struct MatchScoreList& matchScoreList, struct PhmmVector& currentPhmmVector, struct SequenceSegment& currentSequenceSegment) {
  #pragma HLS INLINE
  matchScoreListGenLoop :
  for (uint32_t sequenceSymbolIndex = 0; sequenceSymbolIndex < NUM_CELL_PROCESSORS; sequenceSymbolIndex++) {
    #pragma HLS unroll
    //    ap_uint<2> matchSymbol = currentSequenceSegment.symbols.range((sequenceSymbolIndex * 2) + 2 - 1, sequenceSymbolIndex * 2);
    uint32_t sequenceSegmentWordIndex = sequenceSymbolIndex / SYMBOLS_PER_SEQUENCE_SEGMENT_WORD;
    uint32_t symbolInSequenceSegmentWord = sequenceSymbolIndex % SYMBOLS_PER_SEQUENCE_SEGMENT_WORD;
    ap_uint<2> matchSymbol = currentSequenceSegment.symbols[sequenceSegmentWordIndex].range((symbolInSequenceSegmentWord * 2) + 2 - 1, symbolInSequenceSegmentWord * 2);
    matchScoreList.scores[sequenceSymbolIndex] = currentPhmmVector.scores[matchSymbol];

	#ifdef HAVAC_PER_CELL_DATA_TESTING
    	_symbols[sequenceSymbolIndex] = matchSymbol;
	#endif
  }
}


ap_uint<8> readScoreFromScoreQueue(hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE>& scoreQueue, bool isFirstSequenceSegment, bool isLastPhmmIndex) {
  if (isFirstSequenceSegment || isLastPhmmIndex) {
    return ap_uint<8>(0);
  }
  else {
    return scoreQueue.read();
  }
}

void writeScoreToScoreQueue(hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE>& scoreQueue, ap_uint<8>& scoreToWrite, bool isLastSequenceSegment, bool isFirstPhmmIndex) {
  #pragma HLS pipeline II=1
  ap_uint<8> bufferedScore = scoreToWrite;
  if (!isLastSequenceSegment && !isFirstPhmmIndex) {
    scoreQueue.write(bufferedScore);
  }
}


//this is the single dumbest function I have every written, and I hope it brings pain to anyone who reads it.
//all code in a hls dataflow must be either a variable declaration or a void-returning function, so to
// determine if we're on the first segment, we need to check for zero inside a function, and return by reference.
void isFirstOrLastSequenceSegment(uint32_t sequenceSegmentIndex, uint32_t sequenceLengthInSegments, bool& isFirstSequenceSegment, bool& isLastSequenceSegment) {
  isFirstSequenceSegment = sequenceSegmentIndex == 0;
  isLastSequenceSegment = sequenceSegmentIndex == (sequenceLengthInSegments - 1);
}


//this is also dumb, see makeIsFirstSequenceSegment
void isFirstOrLastPhmmVector(uint32_t phmmVectorIndex, uint32_t phmmLengthInVectors, bool& isFirstPhmmVector, bool& isLastPhmmVector) {
  isFirstPhmmVector = phmmVectorIndex == 0;
  isLastPhmmVector = phmmVectorIndex == (phmmLengthInVectors - 1);
}


void setSequenceSegment(struct SequenceSegment& currentSequenceSegment, hls::stream<struct SequenceSegment, SEQUENCE_STREAM_DEPTH>& sequenceSegmentStream) {
  currentSequenceSegment = sequenceSegmentStream.read();
}


void setPhmmFromStream(struct PhmmVector& phmmVector, hls::stream<struct PhmmVector, PHMM_STREAM_DEPTH>& phmmVectorStream) {
  phmmVector = phmmVectorStream.read();
}
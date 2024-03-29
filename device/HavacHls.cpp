#include "HavacHls.hpp"

#ifdef HAVAC_PER_CELL_DATA_TESTING
#include "../test/byCellComparator/byCellComparator.hpp"
#endif

#ifdef HAVAC_PER_CELL_DATA_TESTING
static uint32_t _phmmIndex;
static uint32_t _sequenceSegmentIndex;
static uint32_t _cellIndexInSegment;
static int8_t _phmmVector[4];
static uint8_t _symbols[NUM_CELL_PROCESSORS];
#endif


hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE> scoreQueueStream("scoreQueueStream");

//top level function that is invoked to run HAVAC. This project uses Vitis HLS, please read the documentation for this tech.
//https://docs.xilinx.com/r/en-US/ug1399-vitis-hls/Getting-Started-with-Vitis-HLS
void HavacKernel(SequenceSegmentWord* sequenceSegmentMemory,
	uint32_t sequenceLengthInSegments, uint32_t* phmmVectorMemory,
	uint32_t phmmLengthInVectors, uint64_t* hitReportMemory,
	uint32_t* hitReportCountMemory) {
	// note: bundle names should be all lowercase
	// depth is just a simulation suggestion, for largest fifo needed for co-simulation
	//seq and phmm lengths are scalars, so they are set via the axilite interface
	//the pointer values are memory partitions in RAM, and use AXI interfaces behind the scenes
	//https://docs.xilinx.com/r/en-US/ug1399-vitis-hls/pragma-HLS-interface

	constexpr uint32_t PHMM_MEM_DEPTH = NUM_CELL_PROCESSORS * TEST_NUM_SEQUENCE_SEGMENTS;
	constexpr uint32_t SEQUENCE_MEM_DEPTH = (NUM_CELL_PROCESSORS / SYMBOLS_PER_SEQUENCE_SEGMENT_WORD) * TEST_NUM_SEQUENCE_SEGMENTS;
	constexpr uint32_t HIT_REPORT_MEM_DEPTH = 256;
	#pragma HLS INTERFACE mode = m_axi port = sequenceSegmentMemory	bundle = gmem0 depth = SEQUENCE_MEM_DEPTH	num_read_outstanding = 16 name=SequenceSegmentMemory
	#pragma HLS INTERFACE mode = m_axi port = phmmVectorMemory 		bundle = gmem1 depth = PHMM_MEM_DEPTH		num_read_outstanding = 128 name=PhmmVectorMemory
	#pragma HLS INTERFACE mode = m_axi port = hitReportMemory 		bundle = gmem0 depth = HIT_REPORT_MEM_DEPTH num_write_outstanding = 16 name=HitReportMemory
	#pragma HLS INTERFACE mode = m_axi port = hitReportCountMemory 	bundle = gmem0 depth = 1 num_write_outstanding = 16 name=HitReportCountMemory

	//massage the length types into their smaller bit-length versions, might be unnecessary, but here for explicity
	seqSegPos_t sequenceLengthAsPos_t = sequenceLengthInSegments;
	phmmPos_t phmmLengthAsPos_t = phmmLengthInVectors;
	HavacTopLevelDataflow(sequenceLengthAsPos_t, phmmLengthAsPos_t, sequenceSegmentMemory,
		phmmVectorMemory, hitReportMemory, hitReportCountMemory);
}

void HavacTopLevelDataflow(const seqSegPos_t sequenceLengthInSegments, const phmmPos_t phmmLengthInVectors,
	SequenceSegmentWord* sequenceSegmentMemory, uint32_t* phmmVectorMemory,
	uint64_t* hitReportMemory, uint32_t* hitReportCountMemory) {

	#pragma HLS DATAFLOW

	//hit reporting streams. There are here because they need run asynchronously from the main design.
	hls::stream<PositionReport, inputHitReportStreamDepth> inputPositionReportStream("inputPositionReportStream");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream0("inputHitReportGroupStream0");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream1("inputHitReportGroupStream1");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream2("inputHitReportGroupStream2");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream3("inputHitReportGroupStream3");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream4("inputHitReportGroupStream4");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream5("inputHitReportGroupStream5");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream6("inputHitReportGroupStream6");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream7("inputHitReportGroupStream7");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream8("inputHitReportGroupStream8");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream9("inputHitReportGroupStream9");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream10("inputHitReportGroupStream10");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream11("inputHitReportGroupStream11");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream12("inputHitReportGroupStream12");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream13("inputHitReportGroupStream13");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream14("inputHitReportGroupStream14");
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream15("inputHitReportGroupStream15");

	hls::stream<PositionReport, intermediateHitReportStreamDepth> hitReportPositionReportStream0("hitReportPositionReportStream0");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_0("hitReportGroupBitsStream0_0");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_1("hitReportGroupBitsStream0_1");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_2("hitReportGroupBitsStream0_2");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_3("hitReportGroupBitsStream0_3");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_4("hitReportGroupBitsStream0_4");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_5("hitReportGroupBitsStream0_5");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_6("hitReportGroupBitsStream0_6");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_7("hitReportGroupBitsStream0_7");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_8("hitReportGroupBitsStream0_8");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_9("hitReportGroupBitsStream0_9");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_10("hitReportGroupBitsStream0_10");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_11("hitReportGroupBitsStream0_11");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_12("hitReportGroupBitsStream0_12");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_13("hitReportGroupBitsStream0_13");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_14("hitReportGroupBitsStream0_14");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_15("hitReportGroupBitsStream0_15");

	hls::stream<PositionReportTier1, intermediateHitReportStreamDepth> hitReportPositionReportStream1("hitReportPositionReportStream1");
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream1("hitReportGroupBitsStream1");
	hls::stream<HitReportTier2, intermediateHitReportStreamDepth> hitReportTier2Stream("hitReportTier2Stream");
	hls::stream<HitReportTier3, intermediateHitReportStreamDepth> hitReportTier3Stream("hitReportTier3Stream");


	HavacMainLoop(sequenceLengthInSegments, phmmLengthInVectors,
		sequenceSegmentMemory, phmmVectorMemory,
		inputPositionReportStream, inputHitReportGroupStream0, inputHitReportGroupStream1,
		inputHitReportGroupStream2, inputHitReportGroupStream3, inputHitReportGroupStream4,
		inputHitReportGroupStream5, inputHitReportGroupStream6, inputHitReportGroupStream7,
		inputHitReportGroupStream8, inputHitReportGroupStream9, inputHitReportGroupStream10,
		inputHitReportGroupStream11, inputHitReportGroupStream12, inputHitReportGroupStream13,
		inputHitReportGroupStream14, inputHitReportGroupStream15);


	filterInputHitReports(inputPositionReportStream, inputHitReportGroupStream0, inputHitReportGroupStream1,
		inputHitReportGroupStream2, inputHitReportGroupStream3, inputHitReportGroupStream4,
		inputHitReportGroupStream5, inputHitReportGroupStream6, inputHitReportGroupStream7,
		inputHitReportGroupStream8, inputHitReportGroupStream9, inputHitReportGroupStream10,
		inputHitReportGroupStream11, inputHitReportGroupStream12, inputHitReportGroupStream13,
		inputHitReportGroupStream14, inputHitReportGroupStream15, hitReportPositionReportStream0,
		hitReportGroupBitsStream0_0, hitReportGroupBitsStream0_1, hitReportGroupBitsStream0_2,
		hitReportGroupBitsStream0_3, hitReportGroupBitsStream0_4, hitReportGroupBitsStream0_5,
		hitReportGroupBitsStream0_6, hitReportGroupBitsStream0_7, hitReportGroupBitsStream0_8,
		hitReportGroupBitsStream0_9, hitReportGroupBitsStream0_10, hitReportGroupBitsStream0_11,
		hitReportGroupBitsStream0_12, hitReportGroupBitsStream0_13, hitReportGroupBitsStream0_14,
		hitReportGroupBitsStream0_15);
	filterHitReportTier0(hitReportPositionReportStream0,
		hitReportGroupBitsStream0_0, hitReportGroupBitsStream0_1, hitReportGroupBitsStream0_2,
		hitReportGroupBitsStream0_3, hitReportGroupBitsStream0_4, hitReportGroupBitsStream0_5,
		hitReportGroupBitsStream0_6, hitReportGroupBitsStream0_7, hitReportGroupBitsStream0_8,
		hitReportGroupBitsStream0_9, hitReportGroupBitsStream0_10, hitReportGroupBitsStream0_11,
		hitReportGroupBitsStream0_12, hitReportGroupBitsStream0_13, hitReportGroupBitsStream0_14,
		hitReportGroupBitsStream0_15, hitReportPositionReportStream1, hitReportGroupBitsStream1);
	filterHitReportTier1(hitReportPositionReportStream1, hitReportGroupBitsStream1,
		hitReportTier2Stream);
	filterHitReportTier2(hitReportTier2Stream, hitReportTier3Stream);
	filterHitReportTier3(hitReportTier3Stream, hitReportMemory, hitReportCountMemory);

}

//main computation loop for HAVAC, each iteration computes a wide column down the Dynamic Programming matrix.
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
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_15) {

HavacSequenceSegmentLoop: for (seqSegPos_t sequenceSegmentIndex = 0;
	sequenceSegmentIndex < sequenceLengthInSegments; sequenceSegmentIndex++) {
	HavacDataflowFunction(sequenceLengthInSegments, phmmLengthInVectors,
		sequenceSegmentMemory, phmmVectorMemory, sequenceSegmentIndex,
		inputPositionReportStream, inputHitReportGroupStream_0, inputHitReportGroupStream_1,
		inputHitReportGroupStream_2, inputHitReportGroupStream_3, inputHitReportGroupStream_4,
		inputHitReportGroupStream_5, inputHitReportGroupStream_6, inputHitReportGroupStream_7,
		inputHitReportGroupStream_8, inputHitReportGroupStream_9, inputHitReportGroupStream_10,
		inputHitReportGroupStream_11, inputHitReportGroupStream_12, inputHitReportGroupStream_13,
		inputHitReportGroupStream_14, inputHitReportGroupStream_15);
}
}



void HavacDataflowFunction(const seqSegPos_t sequenceLengthInSegments, const phmmPos_t phmmLengthInVectors,
	SequenceSegmentWord* sequenceSegmentMemory, uint32_t* phmmVectorMemory, seqSegPos_t sequenceSegmentIndex,
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
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_15) {
	#ifdef HAVAC_PER_CELL_DATA_TESTING
	_sequenceSegmentIndex = sequenceSegmentIndex;
	#endif

	//https://docs.xilinx.com/r/en-US/ug1399-vitis-hls/pragma-HLS-dataflow
	#pragma HLS DATAFLOW

	//streams to act as producer/consumer models to more efficiently read and write data to/from RAM
	hls::stream<SequenceSegmentWord, SEQUENCE_STREAM_DEPTH> sequenceSegmentStream("sequenceStream");
	hls::stream<uint32_t, PHMM_STREAM_DEPTH> phmmStream("phmmStream");
	#pragma HLS DISAGGREGATE variable= sequenceSegmentStream
	#pragma HLS array_partition variable=sequenceSegmentStream type=complete

	//determine if this the first or last sequence segment
	loadSequenceSegmentStream(sequenceSegmentMemory, sequenceSegmentStream,
		sequenceSegmentIndex);
	loadPhmmStream(phmmVectorMemory, phmmStream, phmmLengthInVectors);

	bool isFirstSequenceSegment, isLastSequenceSegment;
	isFirstOrLastSequenceSegment(sequenceSegmentIndex, sequenceLengthInSegments,
		isFirstSequenceSegment, isLastSequenceSegment);

	phmmVectorLoop(phmmLengthInVectors, sequenceSegmentStream, phmmStream,
		isFirstSequenceSegment, isLastSequenceSegment,
		sequenceSegmentIndex, inputPositionReportStream,
		inputHitReportGroupStream_0, inputHitReportGroupStream_1, inputHitReportGroupStream_2,
		inputHitReportGroupStream_3, inputHitReportGroupStream_4, inputHitReportGroupStream_5,
		inputHitReportGroupStream_6, inputHitReportGroupStream_7, inputHitReportGroupStream_8,
		inputHitReportGroupStream_9, inputHitReportGroupStream_10, inputHitReportGroupStream_11,
		inputHitReportGroupStream_12, inputHitReportGroupStream_13, inputHitReportGroupStream_14,
		inputHitReportGroupStream_15);
}


//computes a full column down the DP matrix, going through the whole phmm.
void phmmVectorLoop(phmmPos_t phmmLengthInVectors,
	hls::stream<SequenceSegmentWord, SEQUENCE_STREAM_DEPTH>& sequenceSegmentStream,
	hls::stream<uint32_t, PHMM_STREAM_DEPTH>& phmmStream, bool isFirstSequenceSegment,
	bool isLastSequenceSegment, seqSegPos_t sequenceSegmentIndex,
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
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth>& inputHitReportGroupStream_15) {

	static ap_uint<8> cellScores[NUM_CELL_PROCESSORS];
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

	struct SequenceSegment currentSequenceSegment;

	sequenceConstructionLoop: for (uint8_t sequenceWordIndex = 0; sequenceWordIndex < NUM_SEQUENCE_SEGMENT_WORDS;
		sequenceWordIndex++) {
		currentSequenceSegment.words[sequenceWordIndex] = sequenceSegmentStream.read();
	}

	HavacPhmmVectorLoop: for (phmmPos_t phmmIndex = 0; phmmIndex < phmmLengthInVectors; phmmIndex++) {
	#pragma HLS PIPELINE II=1
	//https://docs.xilinx.com/r/en-US/ug1399-vitis-hls/pragma-HLS-pipeline

	#ifdef HAVAC_PER_CELL_DATA_TESTING
	_phmmIndex = phmmIndex;
	#endif

	const phmmPos_t phmmIndexLocalCopy = phmmIndex;
	//create a match score list, i.e., for each cell, what score should they be given from the phmm vector?
	struct MatchScoreList matchScoreList;
	uint32_t prefetchPhmmVector;
	// vector of bools where a 1 indicates the the corresponding cell passed its threshold this cycle.
	ap_uint<CELLS_PER_GROUP> cellsPassingThreshold[NUM_CELL_GROUPS];
	#pragma HLS array_partition variable=matchScoreList type=complete
	#pragma HLS array_partition variable=cellsPassingThreshold type=complete

	bool isFirstPhmmIndex, isLastPhmmIndex;
	isFirstOrLastPhmmVector(phmmIndexLocalCopy, phmmLengthInVectors,
		isFirstPhmmIndex, isLastPhmmIndex);
	setPhmmFromStream(prefetchPhmmVector, phmmStream);

	bool hitReportTerminator = isLastPhmmIndex && isLastSequenceSegment;

	//make a copy of the phmm vector, hopefully this will help with fanout.
	uint32_t currentPhmmVector = prefetchPhmmVector;

	#ifdef HAVAC_PER_CELL_DATA_TESTING
	ap_uint<32> vectorAsApUint(currentPhmmVector);
	_phmmVector[0] = vectorAsApUint(8 - 1, 0);
	_phmmVector[1] = vectorAsApUint(16 - 1, 8);
	_phmmVector[2] = vectorAsApUint(24 - 1, 16);
	_phmmVector[3] = vectorAsApUint(32 - 1, 24);
	#endif

	generateMatchScoreList(matchScoreList, currentPhmmVector, currentSequenceSegment);

	writeScoreToScoreQueue(bufferedScoreQueueWrite, isLastSequenceSegment, isFirstPhmmIndex);
	computeAllCellProcessors(cellScores, matchScoreList, cellsPassingThreshold,
		bufferedScoreQueueRead, isFirstPhmmIndex, bufferedScoreQueueWrite);

	bufferedScoreQueueRead = readScoreFromScoreQueue(isFirstSequenceSegment, isLastPhmmIndex);

	//in order to reduce the amount of work done on a single cycle, we enqueue all hit reports, even those with no hits
	// to an input queue. this queue is then read, and any reports that actually contain a hit (or terminator) are passed
	//along down the filter chain.

	enqueueHitReportToInputQueue(sequenceSegmentIndex, phmmIndexLocalCopy, cellsPassingThreshold, hitReportTerminator,
		inputPositionReportStream, inputHitReportGroupStream_0, inputHitReportGroupStream_1,
		inputHitReportGroupStream_2, inputHitReportGroupStream_3, inputHitReportGroupStream_4,
		inputHitReportGroupStream_5, inputHitReportGroupStream_6, inputHitReportGroupStream_7,
		inputHitReportGroupStream_8, inputHitReportGroupStream_9, inputHitReportGroupStream_10,
		inputHitReportGroupStream_11, inputHitReportGroupStream_12, inputHitReportGroupStream_13,
		inputHitReportGroupStream_14, inputHitReportGroupStream_15);
	}
}

void computeAllCellProcessors(ap_uint<8> cellScores[NUM_CELL_PROCESSORS], struct MatchScoreList& matchScoreList,
	ap_uint<CELLS_PER_GROUP> cellsPassingThreshold[NUM_CELL_GROUPS], ap_uint<8> leftScoreIn, bool isFirstPhmmIndex,
	ap_uint<8>& lastScoreOut) {

	#ifdef HAVAC_PER_CELL_DATA_TESTING
	_cellIndexInSegment = NUM_CELL_PROCESSORS - 1;
	#endif

	//handle the final cell, since that'll need to save its score to lastScoreOut
	//this seems to need to be broken out to meet II=1. otherwise, vitis tries to schedule the second access to the final cell
	//on the next clock cycle
	ap_uint<8> thisCellPrevScore = isFirstPhmmIndex ? ap_uint<8>(0) : cellScores[NUM_CELL_PROCESSORS - 2];
	struct CellResult result = computeCellProcessor(thisCellPrevScore, matchScoreList.scores[NUM_CELL_PROCESSORS - 1]);
	lastScoreOut = result.cellScore;
	cellsPassingThreshold[NUM_CELL_GROUPS - 1].set_bit(CELLS_PER_GROUP - 1, result.passesThreshold);

	//iterate through the cell processors in reverse, since we need to use cell index i-1, and write to index i.
	CellProcessorLoop: for (uint32_t cellIndex = NUM_CELL_PROCESSORS - 2; cellIndex != 0; cellIndex--) {

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
void loadPhmmStream(uint32_t* phmmVectorMemory, hls::stream<uint32_t, PHMM_STREAM_DEPTH>& phmmVectorStream,
	const phmmPos_t phmmLengthInVectors) {
	phmmReadInnerLoop: for (phmmPos_t phmmIndex = 0; phmmIndex < phmmLengthInVectors; phmmIndex++) {
		phmmVectorStream << phmmVectorMemory[phmmIndex];
	}
}

struct CellResult computeCellProcessor(ap_uint<8> prevScore,
	ap_uint<8> matchScore) {
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
	.symbol = _symbols[_cellIndexInSegment], .passesThreshold = result.passesThreshold };
	addKvToHardwareMap(key, value);
	#endif

	return result;
}


//generates the list of match scores for all cells in the pipeline.
//That is, for every cell processor, find the value from the phmm to add for that cell.
void generateMatchScoreList(struct MatchScoreList& matchScoreList, uint32_t currentPhmmVector,
	struct SequenceSegment currentSequenceSegment) {
	#pragma HLS INLINE

	ap_uint<32> phmmAsApUint(currentPhmmVector);
	ap_uint<32> phmmDuplicateList[NUM_CELL_PROCESSORS];
	phmmDuplicateListGenLoop: for (uint32_t i = 0; i < NUM_CELL_PROCESSORS; i++) {
		phmmDuplicateList[i] = phmmAsApUint;
	}

	matchScoreListGenLoop: for (uint32_t sequenceSymbolIndex = 0;
		sequenceSymbolIndex < NUM_CELL_PROCESSORS; sequenceSymbolIndex++) {
		#pragma HLS unroll

		uint32_t sequenceSegmentWordIndex = sequenceSymbolIndex / SYMBOLS_PER_SEQUENCE_SEGMENT_WORD;
		uint32_t symbolInSequenceSegmentWord = sequenceSymbolIndex % SYMBOLS_PER_SEQUENCE_SEGMENT_WORD;

		ap_uint<2> matchSymbol = currentSequenceSegment.words[sequenceSegmentWordIndex].range(
			(symbolInSequenceSegmentWord * 2) + 1, symbolInSequenceSegmentWord * 2);

		ap_int<8> matchScore;

		switch (matchSymbol) {
		case 0:
			matchScore = phmmDuplicateList[sequenceSymbolIndex](8 - 1, 0);
			break;
		case 1:
			matchScore = phmmDuplicateList[sequenceSymbolIndex](16 - 1, 8);
			break;
		case 2:
			matchScore = phmmDuplicateList[sequenceSymbolIndex](24 - 1, 16);
			break;
		case 3:
			matchScore = phmmDuplicateList[sequenceSymbolIndex](32 - 1, 24);
			break;
		}
		matchScoreList.scores[sequenceSymbolIndex] = matchScore;

		#ifdef HAVAC_PER_CELL_DATA_TESTING
		_symbols[sequenceSymbolIndex] = matchSymbol;
		#endif
	}
}

ap_uint<8> readScoreFromScoreQueue(bool isFirstSequenceSegment, bool isLastPhmmIndex) {
	if (isFirstSequenceSegment || isLastPhmmIndex) {
		return ap_uint<8>(0);
	}
	else {
		return scoreQueueStream.read();
	}
}

void writeScoreToScoreQueue(ap_uint<8> scoreToWrite, bool isLastSequenceSegment, bool isFirstPhmmIndex) {
	#pragma HLS pipeline II=1
	if ((!isLastSequenceSegment && !isFirstPhmmIndex)) {
		scoreQueueStream.write(scoreToWrite);
	}
}

//this is also dumb, see makeIsFirstSequenceSegment
void isFirstOrLastPhmmVector(phmmPos_t phmmVectorIndex, phmmPos_t phmmLengthInVectors,
	bool& isFirstPhmmVector, bool& isLastPhmmVector) {
	isFirstPhmmVector = phmmVectorIndex == 0;
	isLastPhmmVector = phmmVectorIndex == (phmmLengthInVectors - 1);
}


void setPhmmFromStream(uint32_t& phmmVector,
	hls::stream<uint32_t, PHMM_STREAM_DEPTH>& phmmVectorStream) {
	phmmVector = phmmVectorStream.read();
}


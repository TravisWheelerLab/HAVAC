#ifndef HAVAC_HIT_REPORTING_HPP
#define HAVAC_HIT_REPORTING_HPP

#include "PublicDefines.h"
#include <string>
#include <ap_int.h>
#include <hls_stream.h>

#include <stdint.h>

#define inputHitReportStreamDepth 32
#define intermediateHitReportStreamDepth 4


struct PositionReport{
	seqSegPos_t sequenceSegmentIndex;
	phmmPos_t phmmPosition;
	bool terminator;
};

struct PositionReportTier1{
	PositionReport position;
	ap_uint<4> partitionIndex;
};

struct PositionReportTier2{
	PositionReport position;
	ap_uint<8> partitionIndex;
};
struct PositionReportTier3{
	PositionReport position;
	ap_uint<12> partitionIndex;
};
struct PositionReportTier4{
	PositionReport position;
	ap_uint<14> partitionIndex;
};

struct HitReportTier2{
	PositionReportTier2 position;
	ap_uint<CELLS_PER_GROUP/16> hitBits;
};

struct HitReportTier3{
	PositionReportTier3 position;
	ap_uint<CELLS_PER_GROUP/(16*16)> hitBits;
};

struct HitReportTier4{
	PositionReportTier4 position;
	ap_uint<1> hitBit;
};
uint64_t HitReportTier4ToUint64_t(const HitReportTier4 hitReport);


void enqueueHitReportToInputQueue(const seqSegPos_t sequenceSegmentIndex, const phmmPos_t phmmPosition,
	const ap_uint<CELLS_PER_GROUP> hitBits[NUM_CELL_GROUPS], const bool terminator,
	hls::stream<PositionReport, inputHitReportStreamDepth> &inputPositionReportStream,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_0,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_1,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_2,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_3,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_4,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_5,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_6,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_7,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_8,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_9,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_10,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_11,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_12,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_13,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_14,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_15);

void filterInputHitReports(hls::stream<PositionReport, inputHitReportStreamDepth> &inputPositionReportStream,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_0,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_1,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_2,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_3,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_4,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_5,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_6,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_7,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_8,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_9,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_10,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_11,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_12,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_13,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_14,
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_15,
	hls::stream<PositionReport, intermediateHitReportStreamDepth> &hitReportPositionReportStream0,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_0,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_1,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_2,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_3,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_4,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_5,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_6,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_7,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_8,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_9,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_10,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_11,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_12,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_13,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_14,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_15);

void filterHitReportTier0(hls::stream<PositionReport, intermediateHitReportStreamDepth> &hitReportPositionReportStream0,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_0,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_1,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_2,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_3,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_4,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_5,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_6,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_7,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_8,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_9,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_10,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_11,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_12,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_13,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_14,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_15,
	hls::stream<PositionReportTier1, intermediateHitReportStreamDepth> &hitReportPositionReportStream1,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream1);

void filterHitReportTier1(hls::stream<PositionReportTier1, intermediateHitReportStreamDepth> &hitReportPositionReportStream1,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream1,
	hls::stream<HitReportTier2, intermediateHitReportStreamDepth> &hitReportTier2Stream);

void filterHitReportTier2(hls::stream<HitReportTier2, intermediateHitReportStreamDepth> &hitReportTier2Stream,
	hls::stream<HitReportTier3, intermediateHitReportStreamDepth> &hitReportTier3Stream);

void filterHitReportTier3(hls::stream<HitReportTier3, intermediateHitReportStreamDepth> &hitReportTier3Stream,
	uint64_t *hitReportMemory, uint32_t *hitReportCountMemory);

//class HitReporter{
//public:
//	HitReporter(uint64_t *hitReportMemory, uint32_t *hitReportCountMemory);
//	void writeNumHitReports();
//	void enqueueHitReportToInputQueue(const seqSegPos_t sequenceSegmentIndex,
//		const phmmPos_t phmmPosition, const ap_uint<CELLS_PER_GROUP> hitBits[NUM_CELL_GROUPS], const bool terminator);
//
//	void filterInputHitReports();
//	void filterHitReportTier0();
//	void filterHitReportTier1();
//	void filterHitReportTier2();
//	void filterHitReportTier3();
//private:
//	uint32_t *hitReportCountMemory;
//	uint64_t *hitReportMemory;
//	uint32_t numHitReports;
//	hls::stream<PositionReport, inputHitReportStreamDepth> inputPositionReportStream;
//	hls::stream<PositionReport, intermediateHitReportStreamDepth> hitReportPositionReportStream0;
//	hls::stream<PositionReportTier1, intermediateHitReportStreamDepth> hitReportPositionReportStream1;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream1;
//	hls::stream<HitReportTier2, intermediateHitReportStreamDepth> hitReportTier2Stream;
//	hls::stream<HitReportTier3, intermediateHitReportStreamDepth> hitReportTier3Stream;
//
//
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_0;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_1;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_2;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_3;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_4;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_5;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_6;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_7;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_8;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_9;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_10;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_11;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_12;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_13;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_14;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> inputHitReportGroupStream_15;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_0;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_1;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_2;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_3;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_4;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_5;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_6;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_7;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_8;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_9;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_10;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_11;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_12;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_13;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_14;
//	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> hitReportGroupBitsStream0_15;
//};

//
//void filterInputHitReports(HitReporter &hitReporter);
//void filterHitReportTier0(HitReporter &hitReporter);
//void filterHitReportTier1(HitReporter &hitReporter);
//void filterHitReportTier2(HitReporter &hitReporter);
//void filterHitReportTier3(HitReporter &hitReporter);
//void writeHitReports(HitReporter &hitReporter);

//void initNumHitReports();
//void enqueueHitReportToInputQueue(hls::stream<HitReportPartition0, inputHitReportStreamDepth> &inputHitReportStream,
//		const seqSegPos_t sequenceSegmentIndex, const phmmPos_t phmmPosition,
//	const ap_uint<CELLS_PER_GROUP> hitBits[NUM_CELL_GROUPS], const bool terminator);
//void filterPartition0(hls::stream<PositionReport, inputHitReportStreamDepth> &inputPositionReportStream,
//		hls::stream<CELLS_PER_GROUP, inputHitReportStreamDepth> (&inputHitReportGroupStreams)[NUM_CELL_GROUPS],
//		hls::stream<PositionReport, intermediateHitReportStreamDepth> &hitReportPositionReportStream0,
//		hls::stream<CELLS_PER_GROUP, intermediateHitReportStreamDepth> &hitReportBGroupStream0);
//void filterPartition1(hls::stream<PositionReport, intermediateHitReportStreamDepth> &hitReportPositionReportStream0,
//		hls::stream<CELLS_PER_GROUP, intermediateHitReportStreamDepth> (&hitReportStream0Grou0)[NUM_CELL_GROUPS/4],
//		hls::stream<PositionReport, intermediateHitReportStreamDepth> &hitReportPositionReportStream1,
//		hls::stream<HitReportPartition1, intermediateHitReportStreamDepth> &hitReportStream1);
//void filterPartition2(hls::stream<PositionReport, intermediateHitReportStreamDepth> &hitReportPositionReportStream1,
//		hls::stream<HitReportPartition1, intermediateHitReportStreamDepth> &hitReportStream1,
//		hls::stream<HitReportPartition2, intermediateHitReportStreamDepth> &hitReportStream2);
//void filterPartition3(hls::stream<HitReportPartition2, intermediateHitReportStreamDepth> &hitReportStream2,
//		hls::stream<HitReportPartition3, intermediateHitReportStreamDepth> &hitReportStream3);
//void filterPartition4(hls::stream<HitReportPartition3, intermediateHitReportStreamDepth> &hitReportStream3,
//		hls::stream<HitReportPartition4, intermediateHitReportStreamDepth> &hitReportStream4);
//void filterPartition5(hls::stream<HitReportPartition4, intermediateHitReportStreamDepth> &hitReportStream4,
//		hls::stream<HitReportPartition5, intermediateHitReportStreamDepth> &hitReportStream5);
//void filterPartition6(hls::stream<HitReportPartition5, intermediateHitReportStreamDepth> &hitReportStream5,
//		hls::stream<HitReportPartition6, intermediateHitReportStreamDepth> &hitReportStream6);
//void filterPartition7(hls::stream<HitReportPartition6, intermediateHitReportStreamDepth> &hitReportStream6,
//		hls::stream<HitReportPartition7, intermediateHitReportStreamDepth> &hitReportStream7);
//void writeFilteredHitReportsToMemory(hls::stream<HitReportPartition7, intermediateHitReportStreamDepth> &hitReportStream7,
//		uint64_t *hitReportMemory);
//void writeNumHitsToMemory(uint32_t *numHitsMemory);

#endif

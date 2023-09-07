#include "HitReporting.hpp"
#include <hls_stream.h>

#define reportPartitionIndexBitWidth 14
#define reportSeqeunceSegmentIndexBitWidth 26
#define reportPhmmPositionBitWidth 24
#define reportPartitionIndexEndBit reportPartitionIndexBitWidth
#define reportSequenceSegmentIndexEndBit (reportPartitionIndexEndBit + reportSeqeunceSegmentIndexBitWidth)
#define reportPhmmPositionEndBit (reportSequenceSegmentIndexEndBit + reportPhmmPositionBitWidth)


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
	hls::stream<ap_uint<CELLS_PER_GROUP>, inputHitReportStreamDepth> &inputHitReportGroupStream_15){
#pragma HLS PIPELINE II=1
	PositionReport pr;
	pr.phmmPosition = phmmPosition;
	pr.sequenceSegmentIndex = sequenceSegmentIndex;
	pr.terminator = terminator;
	inputPositionReportStream.write(pr);

	//manually unrolled
	inputHitReportGroupStream_0.write(hitBits[0]);
	inputHitReportGroupStream_1.write(hitBits[1]);
	inputHitReportGroupStream_2.write(hitBits[2]);
	inputHitReportGroupStream_3.write(hitBits[3]);
	inputHitReportGroupStream_4.write(hitBits[4]);
	inputHitReportGroupStream_5.write(hitBits[5]);
	inputHitReportGroupStream_6.write(hitBits[6]);
	inputHitReportGroupStream_7.write(hitBits[7]);
	inputHitReportGroupStream_8.write(hitBits[8]);
	inputHitReportGroupStream_9.write(hitBits[9]);
	inputHitReportGroupStream_10.write(hitBits[10]);
	inputHitReportGroupStream_11.write(hitBits[11]);
	inputHitReportGroupStream_12.write(hitBits[12]);
	inputHitReportGroupStream_13.write(hitBits[13]);
	inputHitReportGroupStream_14.write(hitBits[14]);
	inputHitReportGroupStream_15.write(hitBits[15]);
}


void filterInputHitReports(	hls::stream<PositionReport, inputHitReportStreamDepth> &inputPositionReportStream,
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
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_15){
	bool hasSeenTerminator = false;
	ap_uint<CELLS_PER_GROUP> allHitBits[NUM_CELL_GROUPS];

	while(!hasSeenTerminator){
		#pragma HLS PIPELINE II=1
		PositionReport pr =inputPositionReportStream.read();
		hasSeenTerminator = pr.terminator;

		ap_uint<NUM_CELL_GROUPS> groupsContainingHits(0);


		//manually unrolled
		allHitBits[0] = inputHitReportGroupStream_0.read();
		allHitBits[1] = inputHitReportGroupStream_1.read();
		allHitBits[2] = inputHitReportGroupStream_2.read();
		allHitBits[3] = inputHitReportGroupStream_3.read();
		allHitBits[4] = inputHitReportGroupStream_4.read();
		allHitBits[5] = inputHitReportGroupStream_5.read();
		allHitBits[6] = inputHitReportGroupStream_6.read();
		allHitBits[7] = inputHitReportGroupStream_7.read();
		allHitBits[8] = inputHitReportGroupStream_8.read();
		allHitBits[9] = inputHitReportGroupStream_9.read();
		allHitBits[10] = inputHitReportGroupStream_10.read();
		allHitBits[11] = inputHitReportGroupStream_11.read();
		allHitBits[12] = inputHitReportGroupStream_12.read();
		allHitBits[13] = inputHitReportGroupStream_13.read();
		allHitBits[14] = inputHitReportGroupStream_14.read();
		allHitBits[15] = inputHitReportGroupStream_15.read();

		hitReportFIlteringGroupSetBitLoop:for(uint8_t i = 0; i < NUM_CELL_GROUPS;i++){
			#pragma HLS UNROLL
			groupsContainingHits.set_bit(i, allHitBits[i].or_reduce());
		}

		bool containsHit = groupsContainingHits.or_reduce();
		if(containsHit || pr.terminator){

			hitReportPositionReportStream0.write(pr);
			hitReportGroupBitsStream0_0.write(allHitBits[0]);
			hitReportGroupBitsStream0_1.write(allHitBits[1]);
			hitReportGroupBitsStream0_2.write(allHitBits[2]);
			hitReportGroupBitsStream0_3.write(allHitBits[3]);
			hitReportGroupBitsStream0_4.write(allHitBits[4]);
			hitReportGroupBitsStream0_5.write(allHitBits[5]);
			hitReportGroupBitsStream0_6.write(allHitBits[6]);
			hitReportGroupBitsStream0_7.write(allHitBits[7]);
			hitReportGroupBitsStream0_8.write(allHitBits[8]);
			hitReportGroupBitsStream0_9.write(allHitBits[9]);
			hitReportGroupBitsStream0_10.write(allHitBits[10]);
			hitReportGroupBitsStream0_11.write(allHitBits[11]);
			hitReportGroupBitsStream0_12.write(allHitBits[12]);
			hitReportGroupBitsStream0_13.write(allHitBits[13]);
			hitReportGroupBitsStream0_14.write(allHitBits[14]);
			hitReportGroupBitsStream0_15.write(allHitBits[15]);
		}
	}
}


void filterHitReportTier0(
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
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream0_15,
	hls::stream<PositionReportTier1, intermediateHitReportStreamDepth> &hitReportPositionReportStream1,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream1){
	bool hasSeenTerminator = false;
	while(!hasSeenTerminator){
		#pragma HLS PIPELINE off
		PositionReport pr = hitReportPositionReportStream0.read();
		hasSeenTerminator = pr.terminator;

		//manually unrolled
		ap_uint<CELLS_PER_GROUP> hitBits_0 = hitReportGroupBitsStream0_0.read();
		if(hitBits_0.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 0;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_0);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_1 = hitReportGroupBitsStream0_1.read();
		if(hitBits_1.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 1;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_1);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_2 = hitReportGroupBitsStream0_2.read();
		if(hitBits_2.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 2;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_2);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_3 = hitReportGroupBitsStream0_3.read();
		if(hitBits_3.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 3;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_3);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_4 = hitReportGroupBitsStream0_4.read();
		if(hitBits_4.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 4;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_4);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_5 = hitReportGroupBitsStream0_5.read();
		if(hitBits_5.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 5;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_5);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_6 = hitReportGroupBitsStream0_6.read();
		if(hitBits_6.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 6;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_6);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_7 = hitReportGroupBitsStream0_7.read();
		if(hitBits_7.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 7;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_7);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_8 = hitReportGroupBitsStream0_8.read();
		if(hitBits_8.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 8;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_8);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_9 = hitReportGroupBitsStream0_9.read();
		if(hitBits_9.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 9;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_9);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_10 = hitReportGroupBitsStream0_10.read();
		if(hitBits_10.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 10;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_10);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_11 = hitReportGroupBitsStream0_11.read();
		if(hitBits_11.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 11;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_11);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_12 = hitReportGroupBitsStream0_12.read();
		if(hitBits_12.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 12;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_12);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_13 = hitReportGroupBitsStream0_13.read();
		if(hitBits_13.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 13;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_13);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_14 = hitReportGroupBitsStream0_14.read();
		if(hitBits_14.or_reduce()){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = false;
			prTier1.partitionIndex = 14;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_14);
		}
		ap_uint<CELLS_PER_GROUP> hitBits_15 = hitReportGroupBitsStream0_15.read();
		if(hitBits_15.or_reduce()|| pr.terminator){
			PositionReportTier1 prTier1;
			prTier1.position.sequenceSegmentIndex = pr.sequenceSegmentIndex;
			prTier1.position.phmmPosition = pr.phmmPosition;
			prTier1.position.terminator = pr.terminator;
			prTier1.partitionIndex = 15;
			hitReportPositionReportStream1.write(prTier1);
			hitReportGroupBitsStream1.write(hitBits_15);
		}
	}
}


void filterHitReportTier1(
	hls::stream<PositionReportTier1, intermediateHitReportStreamDepth> &hitReportPositionReportStream1,
	hls::stream<ap_uint<CELLS_PER_GROUP>, intermediateHitReportStreamDepth> &hitReportGroupBitsStream1,
	hls::stream<HitReportTier2, intermediateHitReportStreamDepth> &hitReportTier2Stream){
	bool hasSeenTerminator = false;
	while(!hasSeenTerminator){
		PositionReportTier1 prTier1 = hitReportPositionReportStream1.read();
		hasSeenTerminator = prTier1.position.terminator;
		HitReportTier2 hrTier2;
		hrTier2.position.position.sequenceSegmentIndex = prTier1.position.sequenceSegmentIndex;
		hrTier2.position.position.phmmPosition = prTier1.position.phmmPosition;
		ap_uint<CELLS_PER_GROUP> hitBits = hitReportGroupBitsStream1.read();
		for(uint8_t i = 0; i < 16; i++){
			#pragma HLS PIPELINE off
			hrTier2.position.partitionIndex = (prTier1.partitionIndex, ap_uint<4>(i));
			bool thisTerminator = prTier1.position.terminator && (i == 15);
			hrTier2.position.position.terminator = thisTerminator;
			ap_uint<CELLS_PER_GROUP/16> partitionBits = hitBits.range((CELLS_PER_GROUP/16)*(i+1), (CELLS_PER_GROUP/16)*i);
			if(partitionBits.or_reduce() || thisTerminator){
				hrTier2.hitBits = partitionBits;
				hitReportTier2Stream.write(hrTier2);
			}
		}
	}
}

void filterHitReportTier2(
	hls::stream<HitReportTier2, intermediateHitReportStreamDepth> &hitReportTier2Stream,
	hls::stream<HitReportTier3, intermediateHitReportStreamDepth> &hitReportTier3Stream){
	bool hasSeenTerminator = false;
	while(!hasSeenTerminator){
		HitReportTier2 hrTier2 = hitReportTier2Stream.read();
		hasSeenTerminator = hrTier2.position.position.terminator;
		HitReportTier3 hrTier3;
		hrTier3.position.position.sequenceSegmentIndex = hrTier2.position.position.sequenceSegmentIndex;
		hrTier3.position.position.phmmPosition = hrTier2.position.position.phmmPosition;
		for (ap_uint<5> i = 0; i < 16; i++) {
			#pragma HLS PIPELINE off
			bool thisTerminator = hrTier2.position.position.terminator && (i == 15);
			hrTier3.position.partitionIndex = (hrTier2.position.partitionIndex, ap_uint<4>(i));
			hrTier3.position.position.terminator = thisTerminator;
			ap_uint<CELLS_PER_GROUP/(16*16)> partitionBits = hrTier2.hitBits.range((CELLS_PER_GROUP/(16*16))*(i+1), (CELLS_PER_GROUP/(16*16))*i);
			if(partitionBits.or_reduce() || thisTerminator){
				hrTier3.hitBits = partitionBits;
				hitReportTier3Stream.write(hrTier3);
			}
		}
	}
}

void filterHitReportTier3(
	hls::stream<HitReportTier3, intermediateHitReportStreamDepth> &hitReportTier3Stream,
	uint64_t *hitReportMemory, uint32_t *hitReportCountMemory){
	uint32_t numHitReports = 0;
	bool hasSeenTerminator = false;
	while(!hasSeenTerminator){
		HitReportTier3 hrTier3 = hitReportTier3Stream.read();
		HitReportTier4 hrTier4;
		hasSeenTerminator = hrTier3.position.position.terminator;

		hrTier4.position.position.sequenceSegmentIndex = hrTier3.position.position.sequenceSegmentIndex;
		hrTier4.position.position.phmmPosition = hrTier3.position.position.phmmPosition;
		for (ap_uint<2> i = 0; i < 3; i++) {
			#pragma HLS PIPELINE off
			hrTier4.position.partitionIndex = hrTier3.position.partitionIndex;
			hrTier4.position.partitionIndex = (hrTier4.position.partitionIndex << 1) + hrTier4.position.partitionIndex;//multiply by 3
			hrTier4.position.partitionIndex += i;
			bool hitBit = hrTier3.hitBits[i];
			hrTier4.hitBit = hitBit;
			if(hitBit){
				hitReportMemory[numHitReports++] = HitReportTier4ToUint64_t(hrTier4);
			}
		}
	}
	hitReportCountMemory[0] = numHitReports;
}



uint64_t HitReportTier4ToUint64_t(const HitReportTier4 hitReport){
#pragma HLS INLINE
	ap_uint<64> bits;
	bits(14,0) = hitReport.position.partitionIndex;
	bits(40,14) = hitReport.position.position.sequenceSegmentIndex;
	bits(64,40) = hitReport.position.position.phmmPosition;


	return bits.to_uint64();
}

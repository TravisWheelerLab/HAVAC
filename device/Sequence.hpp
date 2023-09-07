#ifndef HAVAC_SEQUENCE_HPP
#define HAVAC_SEQUENCE_HPP

#include <stdint.h>
#include <ap_int.h>
#include "PublicDefines.h"
#include <hls_stream.h>


//HAVAC_TESTING is defined in PublicDefines.h, if it's defined at all
#ifdef HAVAC_TESTING
#define SYMBOLS_PER_SEQUENCE_SEGMENT_WORD 32
#else
#define SYMBOLS_PER_SEQUENCE_SEGMENT_WORD 512
#endif

#define NUM_SEQUENCE_SEGMENT_WORDS (NUM_CELL_PROCESSORS/SYMBOLS_PER_SEQUENCE_SEGMENT_WORD)
#if ((NUM_CELL_PROCESSORS *2) % SYMBOLS_PER_SEQUENCE_SEGMENT_WORD) != 0
#error "bits for whole sequence segment not evenly divisible by the num seq segment words"
#endif

#define SEQUENCE_STREAM_DEPTH (2*NUM_SEQUENCE_SEGMENT_WORDS)

typedef ap_uint<SYMBOLS_PER_SEQUENCE_SEGMENT_WORD*2> SequenceSegmentWord;
struct SequenceSegment {
	SequenceSegmentWord words[NUM_SEQUENCE_SEGMENT_WORDS];
};


void loadSequenceSegmentStream(SequenceSegmentWord* sequenceSegmentMemory, hls::stream<SequenceSegmentWord,
	SEQUENCE_STREAM_DEPTH>& sequenceSegmentStream, seqSegPos_t sequenceSegmentIndex);

void isFirstOrLastSequenceSegment(seqSegPos_t sequenceSegmentIndex, seqSegPos_t sequenceLengthInSegments, bool& isFirstSequenceSegment, bool& isLastSequenceSegment);

#endif

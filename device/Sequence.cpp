#include "Sequence.hpp"

//this is the single dumbest function I have every written, and I hope it brings pain to anyone who reads it.
//all code in a hls dataflow must be either a variable declaration or a void-returning function, so to
// determine if we're on the first segment, we need to check for zero inside a function, and return by reference.
void isFirstOrLastSequenceSegment(seqSegPos_t sequenceSegmentIndex,
  seqSegPos_t sequenceLengthInSegments, bool &isFirstSequenceSegment,
  bool &isLastSequenceSegment) {

	isFirstSequenceSegment = sequenceSegmentIndex == 0;
  isLastSequenceSegment = sequenceSegmentIndex == (sequenceLengthInSegments - 1);
}

void loadSequenceSegmentStream(SequenceSegmentWord *sequenceSegmentMemory,
		hls::stream<SequenceSegmentWord, SEQUENCE_STREAM_DEPTH> &sequenceSegmentStream,
		seqSegPos_t sequenceSegmentIndex) {
	uint32_t baseSequenceWordOffset = (uint32_t)sequenceSegmentIndex
			* NUM_SEQUENCE_SEGMENT_WORDS;
	for (uint8_t sequenceSegmentWordIndex = 0;
			sequenceSegmentWordIndex < NUM_SEQUENCE_SEGMENT_WORDS;
			sequenceSegmentWordIndex++) {
		uint32_t sequenceOffset = baseSequenceWordOffset
				+ sequenceSegmentWordIndex;
		sequenceSegmentStream << sequenceSegmentMemory[sequenceOffset];
	}
}


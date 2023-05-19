#ifndef HAVAC_SCORE_QUEUE_HPP
#define HAVAC_SCORE_QUEUE_HPP

#include <stdint.h>
#include <hls_stream.h>
#include <ap_int.h>

//#define ULTRARAM_BYTES_PER_RAM 36000ULL  // 1 bram is 288Kb, or 36000 bytes.
//#define NUM_ULTRARAMS 120ULL

//#define SCORE_QUEUE_SIZE (512*1024)	//this is used for Csim/cosim, since the real size blows the software stack size
#define SCORE_QUEUE_SIZE (1024*1024)

class ScoreQueue{
public:
	ScoreQueue();
	ap_uint<8> read();
	void write(const ap_uint<8> score);
private:
	ap_uint<8> memory[SCORE_QUEUE_SIZE];
	uint32_t scoreQueueWritePosition;
	uint32_t scoreQueueReadPosition;
};


class FifoScoreQueue{
public:
	FifoScoreQueue();
	ap_uint<8> read();
	void write(const ap_uint<8> score);
private:
	hls::stream<ap_uint<8>, SCORE_QUEUE_SIZE> stream;
};

#endif

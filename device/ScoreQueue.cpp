#include "ScoreQueue.hpp"
#include "HavacHls.hpp"


ScoreQueue::ScoreQueue(){
#pragma HLS bind_storage variable=memory type=RAM_S2P impl=BRAM
	this->scoreQueueWritePosition = 0;
	this->scoreQueueReadPosition = 0;
}

ap_uint<8> ScoreQueue::read(){
#pragma HLS pipeline II=1
	ap_uint<8> readScore = this->memory[this->scoreQueueReadPosition];
	this->scoreQueueReadPosition = (this->scoreQueueReadPosition + 1) % SCORE_QUEUE_SIZE;
	return readScore;
}

void ScoreQueue::write(const ap_uint<8> score){
#pragma HLS pipeline II=1
	this->memory[this->scoreQueueWritePosition] = score;
	this->scoreQueueWritePosition = (this->scoreQueueWritePosition + 1) % SCORE_QUEUE_SIZE;
}

FifoScoreQueue::FifoScoreQueue():
	stream("FifoScoreQueueStream"){
	  #pragma HLS bind_storage variable=stream type=fifo impl=bram
}

ap_uint<8> FifoScoreQueue::read(){
	ap_uint<8> readValue;
	this->stream.read_nb(readValue);
	return readValue;
}

void FifoScoreQueue::write(const ap_uint<8> score){
	this->stream.write_nb(score);
}

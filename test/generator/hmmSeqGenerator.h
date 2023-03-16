#pragma once

#include <cstdint>
extern "C" {
#include "../../P7HmmReader/src/p7HmmReader.h"
}

struct HmmSeqPair{
  char *sequence;
  uint32_t sequenceLength;
  P7HmmList *phmmList;
};


///@brief generates a hmm for a random sequence of the given length, and a mutated 
///version of that sequence, and returns them in a HmmSeqPair
///@param seqLength length of the sequence to generate
///@return hmm and sequence in-memory
struct HmmSeqPair generateRandomHmmSeqPair(const uint32_t seqLength);

///@brief  generates a hmm for a random sequence of the given length, and a mutated
///version of that sequence, and saves them to the files at the given sources.
///@param seqLength lenght of the sequence to generate
///@param seqSrc location to save the .fasta sequence file to.
///@param phmmSrc location to save the .hmm model file to.
void generateRandomHmmSeqPairToFiles(const uint32_t seqLength, const char* seqSrc, const char* phmmSrc);
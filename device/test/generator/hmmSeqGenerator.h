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


struct HmmSeqPair readPregeneratedHmmSeqPair();
struct HmmSeqPair generateRandomHmmSeqPair(const uint32_t seqLength);
void writeSequenceToFile(struct HmmSeqPair* hmmSeqPair, const char *referenceSequenceSrc);
char *readSequenceFromFile(uint32_t sequenceLength, const char *referenceSequenceSrc);
struct P7HmmList* readHmm(const char* outputHmmFileSrc);

void generateSequencesToFile(const uint32_t seqLength, char *sequenceFileSrc, char *mutatedSeqFileSrc, const float sequenceMutationSubProbability);
void generatePhmmToFile(char *seqFileSrc, char *phmmFileSrc);

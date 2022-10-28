#include "hmmSeqGenerator.h"
extern "C" {
#include "../../P7HmmReader/src/p7HmmReader.h"
#include "../../P7HmmReader/src/p7ProfileHmm.h"
#include "../../P7HmmReader/src/p7HmmReaderLog.h"
#include "../../FastaVector/src/FastaVector.h"
}
//#include <p7HmmReader.h>
//#include <p7ProfileHmm.h>
//#include <p7HmmReaderLog.h>
//#include <FastaVector.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

const float sequenceMutationSubProbability = 0.2f;

const char* tmpFastaFileSrc = "tmpHavacTest.fasta";
const char* outputHmmFileSrc = "tmpHavacTest.hmm";
const char nucs[4] = { 'a', 'c', 'g', 't' };



//private function prototypes
void randomlyAssignSequence(char* sequence, const uint32_t seqLength);
void writeSequenceToFasta(char* sequence, const char* fastaFileSrc);
void generateHmmFromFasta(const char* fastaFileSrc, const char* outputHmmFileSrc);
struct P7HmmList* readHmm(const char* outputHmmFileSrc);
void mutateSequence(char* sequence, uint32_t sequenceLength, float probability);
void hmmSequenceGeneratorCleanupTempFiles(void);


struct HmmSeqPair generateRandomHmmSeqPair(const uint32_t seqLength) {
  //initialize the RNG seed
  srand(time(0));
  struct HmmSeqPair hmmSeqPair;
  hmmSeqPair.sequenceLength = seqLength;
  hmmSeqPair.sequence = (char*)malloc(seqLength + 1); //alloc an additiona byte for null terminator
  if(hmmSeqPair.sequence == NULL){
    printf("ALLOCATION FAILURE: could not allocate memory for sequence in hmmSeqPair\n");
  }
  printf("generating randomly assigned sequence w/ esl...\n");
  fflush(stdout);
  randomlyAssignSequence(hmmSeqPair.sequence, seqLength);
  writeSequenceToFasta(hmmSeqPair.sequence, tmpFastaFileSrc);
  printf("generating hmm from sequence...\n");
  fflush(stdout);
  generateHmmFromFasta(tmpFastaFileSrc, outputHmmFileSrc);
  hmmSeqPair.phmmList = readHmm(outputHmmFileSrc);
  printf("mutating sequence with %f sub rate\n", sequenceMutationSubProbability);
  fflush(stdout);
  mutateSequence(hmmSeqPair.sequence, hmmSeqPair.sequenceLength, sequenceMutationSubProbability);
  printf("cleaning up temp files\n");
  fflush(stdout);
  hmmSequenceGeneratorCleanupTempFiles();

  return hmmSeqPair;
}


void randomlyAssignSequence(char* sequence, const uint32_t seqLength) {
  for (size_t i = 0; i < seqLength; i++) {
    sequence[i] = nucs[rand() % 4];
  }
  //null terminate the string
  sequence[seqLength] = 0;
}


void writeSequenceToFasta(char* sequence, const char* fastaFileSrc) {
  FILE* f = fopen(fastaFileSrc, "w");
  if (!f) {
    printf("ERROR: could not open fasta file %s for writing\n", fastaFileSrc);
    exit(2);
  }
  int result = fputs(">generatedSequence\n", f);
  if (result < 0) {
    printf("ERROR: could not write header to fasta file\n");
    fclose(f);
    exit(3);
  }
  result = fputs(sequence, f);
  if (result < 0) {
    printf("ERROR: could not write sequence to fasta file\n");
    fclose(f);
    exit(3);
  }
  fclose(f);
}


void generateHmmFromFasta(const char* fastaFileSrc, const char* outputHmmFileSrc) {
  char buffer[4096];
  sprintf(buffer, "hmmbuild %s %s", outputHmmFileSrc, fastaFileSrc);
  system(buffer);
}

struct P7HmmList* readHmm(const char* outputHmmFileSrc) {
  struct P7HmmList *phmmList = (struct P7HmmList*)malloc(sizeof(struct P7HmmList));
  if(phmmList == NULL){
    printf("Error: could not allocate memory for P7HmmList struct\n");
    exit(5);
  }
  printf("invoking readP7Hmm\n");
  enum P7HmmReturnCode rc = readP7Hmm(outputHmmFileSrc, phmmList);

  if (rc != p7HmmSuccess){
    printf("ERROR: could not read phmm file, returned code %i\n", rc);
    fflush(stdout);
    exit(4);
  }
  else{
	  printf("\nreadP7Hmm finished\n");
	    fflush(stdout);
  }
  return phmmList;
}


void mutateSequence(char* sequence, uint32_t sequenceLength, float probability){
  for (uint32_t i = 0; i < sequenceLength; i++) {
    double randProb = (double)rand() / (double)RAND_MAX;
    if(randProb < probability){
      sequence[i] = nucs[rand() % 4];
    }
  }
}

void hmmSequenceGeneratorCleanupTempFiles(void) {
  remove(tmpFastaFileSrc);
  remove(outputHmmFileSrc);
}

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

const float sequenceMutationSubProbability = 0.3f;

const char* tmpFastaFileSrc = "tmpHavacTest.fasta";
const char* outputHmmFileSrc = "tmpHavacTest.hmm";
const char* referenceSequenceSrc = "ref.txt";
const char nucs[4] = { 'a', 'c', 'g', 't' };



//private function prototypes
void randomlyAssignSequence(char* sequence, const uint32_t seqLength);
void writeSequenceToFasta(char* sequence, const char* fastaFileSrc);
void generateHmmFromFasta(const char* fastaFileSrc, const char* outputHmmFileSrc);
void mutateSequence(char* sequence, uint32_t sequenceLength, float probability);
void hmmSequenceGeneratorCleanupTempFiles(void);


struct HmmSeqPair readPregeneratedHmmSeqPair(){
	  struct HmmSeqPair hmmSeqPair;
	  hmmSeqPair.phmmList = readHmm(outputHmmFileSrc);
	  printf("hmm read\n");
	  hmmSeqPair.sequenceLength = hmmSeqPair.phmmList->phmms[0].header.modelLength;
	  hmmSeqPair.sequence =  readSequenceFromFile(hmmSeqPair.sequenceLength , referenceSequenceSrc);
	  return hmmSeqPair;
}

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
  printf("writing to file...\n");
  writeSequenceToFile(&hmmSeqPair, referenceSequenceSrc);
  printf("cleaning up temp files\n");
  fflush(stdout);

  return hmmSeqPair;
}

void generateSequencesToFile(const uint32_t seqLength, char *sequenceFileSrc, char *mutatedSeqFileSrc, const float sequenceMutationSubProbability){
	char *sequence = (char*)malloc((seqLength+1) * sizeof(char));
  randomlyAssignSequence(sequence, seqLength);
  writeSequenceToFasta(sequence, sequenceFileSrc);
  mutateSequence(sequence, seqLength, sequenceMutationSubProbability);
  writeSequenceToFasta(sequence, mutatedSeqFileSrc);
  free(sequence);


}
void generatePhmmToFile(char *seqFileSrc, char *phmmFileSrc){
  generateHmmFromFasta(seqFileSrc, phmmFileSrc);
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



void writeSequenceToFile(struct HmmSeqPair* hmmSeqPair, const char *referenceSequenceSrc){
	printf("writing seq to file\n");
	fflush(stdout);
	FILE *sequenceFile = fopen(referenceSequenceSrc, "w+");
	if(!sequenceFile){
		printf("ERROR: could not open ref sequence file for writing");
		exit(100);
	}
	printf("file opened for writing\n");
	fflush(stdout);
	fprintf(sequenceFile, "%.*s", hmmSeqPair->sequenceLength, hmmSeqPair->sequence);
	printf("seq written\n");
	fflush(stdout);
	fclose(sequenceFile);
}

char *readSequenceFromFile(uint32_t sequenceLength, const char *referenceSequenceSrc){
	FILE *sequenceFile = fopen(referenceSequenceSrc, "r");
	if(!sequenceFile){
		printf("ERROR: could not open ref sequence file for reading");
		exit(102);
	}
	char *sequence = (char*)malloc(sequenceLength * sizeof(char));
	if(!sequence){
		fclose(sequenceFile);
		printf("ERROR: could not allocate sequence for reading the sequence from file\n");
		exit(103);
	}
	fgets(sequence, sequenceLength, sequenceFile);
	fclose(sequenceFile);
	return sequence;

}

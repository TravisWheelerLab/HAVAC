#include "hmmSeqGenerator.h"
#include <time.h>
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
#include <sstream>
#include <time.h>
#include <vector>

// const float sequenceMutationSubProbability = 0.3f;

const char* tmpFastaFileSrc = "tmpHavacTest.fasta";
const char* outputHmmFileSrc = "tmpHavacTest.hmm";
const char nucs[4] = { 'a', 'c', 'g', 't' };
const uint32_t SEQUENCE_FLANK_LENGTH = 50;


//private function prototypes
void randomlyAssignSequence(char* sequence, const uint32_t seqLength);
void writeSequenceToFasta(char* sequence, const char* fastaFileSrc);
void generateHmmFromFasta(const char* fastaFileSrc, const char* outputHmmFileSrc);
struct P7HmmList* readHmm(const char* outputHmmFileSrc);
void mutateSequence(char* sequence, uint32_t sequenceLength, float probability);
void mutateSequenceWithFlanks(char** sequence, uint32_t seqLength, float sequenceMutationSubProbability, uint32_t flankLengths);
void mutateSequenceWithIndels(char** sequence, const uint32_t seqLength, const float subProbability, const float indelProbability);
void hmmSequenceGeneratorCleanupTempFiles(void);


struct HmmSeqPair generateRandomHmmSeqPair(const uint32_t seqLength, const bool addFlankingRegions, float sequenceMutationSubProbability) {
  struct HmmSeqPair hmmSeqPair;
  hmmSeqPair.sequenceLength = seqLength;
  hmmSeqPair.sequence = (char*)malloc(seqLength + 1); //alloc an additiona byte for null terminator
  if (hmmSeqPair.sequence == NULL) {
    printf("ALLOCATION FAILURE: could not allocate memory for sequence in hmmSeqPair\n");
  }
  fflush(stdout);
  randomlyAssignSequence(hmmSeqPair.sequence, seqLength);
  writeSequenceToFasta(hmmSeqPair.sequence, tmpFastaFileSrc);
  fflush(stdout);
  generateHmmFromFasta(tmpFastaFileSrc, outputHmmFileSrc);
  hmmSeqPair.phmmList = readHmm(outputHmmFileSrc);
  if (addFlankingRegions) {
    mutateSequenceWithFlanks(&hmmSeqPair.sequence, hmmSeqPair.sequenceLength, sequenceMutationSubProbability, SEQUENCE_FLANK_LENGTH);
  }
  else {
    mutateSequence(hmmSeqPair.sequence, hmmSeqPair.sequenceLength, sequenceMutationSubProbability);
  }
  fflush(stdout);
  hmmSequenceGeneratorCleanupTempFiles();

  return hmmSeqPair;
}

void generateRandomHmmSeqPairToFiles(const uint32_t seqLength, const char* seqSrc, const char* phmmSrc,
  const bool addFlankingRegions, float sequenceMutationSubProbability) {
  char* sequence = (char*)malloc(seqLength + 1); //alloc an additiona byte for null terminator
  if (sequence == NULL) {
    printf("ALLOCATION FAILURE: could not allocate memory for sequence in hmmSeqPair\n");
  }
  printf("generating randomly assigned sequence w/ esl...\n");
  fflush(stdout);
  randomlyAssignSequence(sequence, seqLength);
  writeSequenceToFasta(sequence, seqSrc);
  printf("generating hmm from sequence...\n");
  fflush(stdout);
  generateHmmFromFasta(seqSrc, phmmSrc);
  printf("mutating sequence with %f sub rate\n", sequenceMutationSubProbability);
  fflush(stdout);
  if (addFlankingRegions) {
    float indelProbability = sequenceMutationSubProbability / 2.0f;
    mutateSequenceWithIndels(&sequence, seqLength, sequenceMutationSubProbability, indelProbability);
  }
  else {
    float indelProbability = sequenceMutationSubProbability / 4.0f;
    mutateSequenceWithIndels(&sequence, seqLength, sequenceMutationSubProbability, indelProbability);
  }
  writeSequenceToFasta(sequence, seqSrc);
  fflush(stdout);
  hmmSequenceGeneratorCleanupTempFiles();
  free(sequence);

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
  sprintf(buffer, "hmmbuild --dna %s %s  > /dev/null", outputHmmFileSrc, fastaFileSrc);
  system(buffer);
}

struct P7HmmList* readHmm(const char* outputHmmFileSrc) {
  struct P7HmmList* phmmList = (struct P7HmmList*)malloc(sizeof(struct P7HmmList));
  if (phmmList == NULL) {
    printf("Error: could not allocate memory for P7HmmList struct\n");
    exit(5);
  }
  printf("invoking readP7Hmm\n");
  enum P7HmmReturnCode rc = readP7Hmm(outputHmmFileSrc, phmmList);

  if (rc != p7HmmSuccess) {
    printf("ERROR: could not read phmm file, returned code %i\n", rc);
    fflush(stdout);
    exit(4);
  }
  else {
    printf("\nreadP7Hmm finished\n");
    fflush(stdout);
  }
  return phmmList;
}


void mutateSequence(char* sequence, uint32_t sequenceLength, float probability) {
  for (uint32_t i = 0; i < sequenceLength; i++) {
    double randProb = (double)rand() / (double)RAND_MAX;
    if (randProb < probability) {
      sequence[i] = nucs[rand() % 4];
    }
  }
}


//adds random pre and post flanking regions to reduce the chance that hits stride the sequence boundary
void mutateSequenceWithFlanks(char** sequence, uint32_t seqLength, float sequenceMutationSubProbability, uint32_t flankLengths) {
  char* newSequence = (char*)malloc(seqLength + (2 * flankLengths) * sizeof(char));
  if (newSequence == NULL) {
    printf("ALLOCATION FAILURE: could not allocate memory for new sequence with flanks\n");
    exit(-7);
  }
  for (uint32_t i = 0; i < flankLengths;i++) {
    newSequence[i] = nucs[rand() % 4];
  }
  for (uint32_t i = 0; i < seqLength; i++) {
    double randProb = (double)rand() / (double)RAND_MAX;
    char thisChar;
    if (randProb < sequenceMutationSubProbability) {
      thisChar = nucs[rand() % 4];
    }
    else {
      thisChar = (*sequence)[i];
    }
    newSequence[i + flankLengths] = thisChar;
  }
  memcpy(newSequence + flankLengths, *sequence, seqLength);
  for (uint32_t i = 0; i < flankLengths;i++) {
    newSequence[flankLengths + seqLength + i] = nucs[rand() % 4];
  }

  //dealloc the old sequence, and change the ptr to the new sequence
  free(*sequence);
  *sequence = newSequence;
}

void hmmSequenceGeneratorCleanupTempFiles(void) {
  remove(tmpFastaFileSrc);
  remove(outputHmmFileSrc);
}

void mutateSequenceWithIndels(char** sequence, const uint32_t seqLength, const float subProbability, const float indelProbability) {
  std::vector<char> newSequence;
  uint32_t sequencePosition = 0;
  while (sequencePosition < seqLength) {
    double indelChance = (double)rand() / (double)RAND_MAX;
    if (indelChance < indelProbability) {
      //randomize whether its a insert or delete
      if (rand() % 2) {
        newSequence.push_back(nucs[rand() % 4]);  //a new character is inserted
        continue;
      }
      else {
        sequencePosition++;//the character from the sequence is deleted
        continue;
      }
    }

    double subChance = (double)rand() / (double)RAND_MAX;
    if (subChance < subProbability) {
      newSequence.push_back(nucs[rand() % 4]);
      sequencePosition++;
    }
    else {
      newSequence.push_back((*sequence)[sequencePosition]);
      sequencePosition++;
    }
  }
  newSequence.push_back('\0');

  free(*sequence);
  *sequence = (char*)malloc(newSequence.size());
  memcpy(*sequence, newSequence.data(), newSequence.size());
}
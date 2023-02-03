#include "softwareTestbench.hpp"
#include "../device/HavacHls.hpp"
#include "../device/PublicDefines.h"
extern "C" {
  #include "../P7HmmReader/src/p7ProfileHmm.h"
  #include "../FastaVector/src/FastaVector.h"
  #include "../FastaVector/src/FastaVectorMetadataVector.h"
  #include "../FastaVector/src/FastaVectorString.h"
  #include "../P7HmmReader/src/p7HmmReader.h"
  #include "../P7HmmReader/src/p7HmmReaderLog.h"
}
#include "generator/hmmSeqGenerator.h"
#include "softSsv/SoftSsv.h"
#include "../PhmmReprojection/PhmmReprojection.h"

#include <cstdint>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <unistd.h>
#include <time.h>
#include <unistd.h>

//defined in PublicDefines.h, if defined at all
#ifdef HAVAC_PER_CELL_DATA_TESTING
#include "byCellComparator/byCellComparator.hpp"
#endif

int mockTestbench();

/* This file contains the testbench that ensure that HAVAC will return the same results as a traditional software SSV implementation.
*/
uint8_t* sequenceAsUnpackedVectorIndices(struct HmmSeqPair& hmmSeqPair);
void packSequenceVectorIndices(uint8_t* sequenceAsVectorIndices, size_t sequenceLength);
// int8_t* makeEmissionsAsInt8_tList(struct HmmSeqPair& hmmSeqPair);
uint8_t nucToEncoding(const uint8_t nuc);
void generateTestSet(const uint32_t numTests, const uint32_t sequenceLengthInSegments, const float sequenceMutProbability);
void makePhmmFileSrcFromTestNum(const uint32_t testNum, char *buffer);
void makeSequenceFileSrcFromTestNum(const uint32_t testNum, char *buffer);
void makeMutatedSequenceFileSrcFromTestNum(const uint32_t testNum, char *buffer);
struct HmmSeqPair readTestSetFromFiles(const uint32_t testNum);

const uint32_t numTests = 1;
const float requiredPValue = 0.05f;
//const bool REQUIRES_INPUT_GENERATION = true;
const char *inputGenFlagFile = "/home/ta112817/tmp/inputGen.flag";

int main() {
  srand(time(0));

  printf("BEGINNING TESTBENCH...\n");
  bool allTestsPassed = true;

  const float sequenceMutProbability = 0.3f;

  //generate inputs if the input file flag was not found.
  FILE *inputGenFile = fopen(inputGenFlagFile, "r");
  if (inputGenFile) {
      // file exists
	  fclose(inputGenFile);
  } else {
      // file doesn't exist
	  printf("generating input sequences/phmms...\n");
	  generateTestSet(numTests, TEST_NUM_SEQUENCE_SEGMENTS, sequenceMutProbability);
	  inputGenFile = fopen(inputGenFlagFile, "w");
	  fputs("exists\n", inputGenFile);
	  fclose(inputGenFile);
  }


//  if(REQUIRES_INPUT_GENERATION){
//  }

  for (size_t i = 1; i <= numTests; i++) {
    #ifdef HAVAC_PER_CELL_DATA_TESTING
    clearCellComparatorMaps();
    #endif

    printf(" Test #%u/%u \n", i, numTests);

    struct HmmSeqPair hmmSeqPair = readTestSetFromFiles(i);

    uint32_t phmmLength = hmmSeqPair.phmmList->phmms[0].header.modelLength;
    uint32_t modelCardinality = 4;
    int8_t *projectedPhmm = new int8_t[phmmLength * modelCardinality];
    if(!projectedPhmm){
      printf("ERROR: could not allocate memory for the projected phmm\n");
      exit(-200);
    }
    printf("reprojecting scores...\n");
    fflush(stdout);
    // reproject the scores
    //the hmm/seq generator makes a phmm with only one model in it, so we can use phmms[0] to reference it in this function.
    if(!hmmSeqPair.phmmList){
    	printf("phmmList was null!\n");
    	fflush(stdout);
    }
    if(!hmmSeqPair.phmmList->phmms){
    	printf("phmms were null!\n");
    	fflush(stdout);
    }

    p7HmmProjectForThreshold256(&hmmSeqPair.phmmList->phmms[0], requiredPValue, projectedPhmm);
    // int8_t* emissionsAsInt8_t = makeEmissionsAsInt8_tList(hmmSeqPair);

    uint8_t* sequenceAsVectorIndices = sequenceAsUnpackedVectorIndices(hmmSeqPair);


    printf("gen sequence done\n");
    fflush(stdout);


    int8_t errorCode = 0;
    std::vector<struct SoftSsvHit> softwareSsvHits = softSsvThreshold256(sequenceAsVectorIndices,
      hmmSeqPair.sequenceLength, projectedPhmm, hmmSeqPair.phmmList->phmms[0].header.modelLength, errorCode);
    if (errorCode != 0) {
      printf("Error: soft ssv returned error code %i\n", errorCode);
      exit(8);
    }


    printf("softSsv ran\n");
    fflush(stdout);
    //pack the sequence into 2-bit values for the hardware HAVAC implementation.
    packSequenceVectorIndices(sequenceAsVectorIndices, hmmSeqPair.sequenceLength);

    uint32_t *phmmAsUint32 = new uint32_t[hmmSeqPair.phmmList->phmms[0].header.modelLength];
	printf("memcpying\n");
	fflush(stdout);
	memcpy(phmmAsUint32, projectedPhmm, hmmSeqPair.phmmList->phmms[0].header.modelLength * sizeof(uint32_t));

    printf("invoking hardware ssv...\n");
    fflush(stdout);
    std::vector<struct HitReport> hardwareSsvHits = invokeHardwareSsv(sequenceAsVectorIndices, hmmSeqPair.sequenceLength,
    	phmmAsUint32, hmmSeqPair.phmmList->phmms[0].header.modelLength, errorCode);
    delete[] phmmAsUint32;
    printf("hardware finished\n");

    printf("soft ssv returned %zu hits\n", softwareSsvHits.size());
    printf("hard ssv returned %zu hits\n", hardwareSsvHits.size());

    #ifdef HAVAC_PER_CELL_DATA_TESTING
    printf("comparing data on per-cell basis\n");
    compareCellMaps(hmmSeqPair.phmmList->phmms[0].header.modelLength, hmmSeqPair.sequenceLength, hmmSeqPair);
    #endif


    fflush(stdout);
    bool comparesSuccessfully = compareSsvHitLists(hardwareSsvHits, softwareSsvHits);
    allTestsPassed &= comparesSuccessfully;
    if (!comparesSuccessfully) {
      printf("Test Fail: test number %zu did not compare successfully.\n", i);
    }

    p7HmmListDealloc(hmmSeqPair.phmmList);
    delete[] sequenceAsVectorIndices;
    delete[] projectedPhmm;
    free(hmmSeqPair.sequence);
  }

  printf("testing finished.\n");
  if (allTestsPassed) {
    printf("all tests passed\n");
    return 0;
  }
  else {
    printf("FAILURE in at least one test.\n");
    return 1;
  }
}


/*invoke the hardware HAVAC toplevel function, and return a vector of the hit reports that HAVAC generated*/
std::vector<struct HitReport> invokeHardwareSsv(const uint8_t* sequenceAsVectorIndices, size_t sequenceLength,
	uint32_t *emissionsAsUInt32_t, const size_t phmmLength, int8_t& errorCode) {
  errorCode = 0;
  std::vector<struct HitReport> hitReportList;
  const size_t sequenceLengthInSegments = sequenceLength / (NUM_CELL_PROCESSORS);
  if (sequenceLength % NUM_CELL_PROCESSORS != 0) {
    printf("ERROR: seq length was %zu was not a multiple of num_cell_processors %zu!\n", sequenceLength, NUM_CELL_PROCESSORS);
    errorCode = -1;
  }


  const uint32_t NUM_HIT_REPORTS_SUPPORTED = phmmLength * sequenceLengthInSegments;
  struct HitReport* hitReportMemory = (struct HitReport*)malloc(sizeof(struct HitReport) * NUM_HIT_REPORTS_SUPPORTED);
  if (hitReportMemory == NULL) {
    printf("ERROR: could not allocate memory for hit report buffer\n");
    errorCode = -2;
    return hitReportList;
  }

  uint32_t numHits;
  HavacKernelTopLevel((SequenceSegmentWord*)sequenceAsVectorIndices, sequenceLengthInSegments,
		  emissionsAsUInt32_t, phmmLength, hitReportMemory, numHits);

  	  printf("    havac kernel finished\n");
  	  fflush(stdout);

  if (numHits > NUM_HIT_REPORTS_SUPPORTED) {
    printf("ERROR: num hits %zu was greater than max supported of %zu\n", numHits, NUM_HIT_REPORTS_SUPPORTED);
    errorCode = -3;
    return hitReportList;
  }
  // transfer over the hit reports
  for (size_t i = 0; i < numHits; i++) {
    hitReportList.push_back(hitReportMemory[i]);
  }

  free(hitReportMemory);

  return hitReportList;
}


/*compares the hardware ssv hits to the software reference. */
bool compareSsvHitLists(std::vector<struct HitReport> hardwareSsvHits, std::vector<struct SoftSsvHit> softwareSsvHits) {
  bool allTestsPass = true;
  //iterate over the software ssv hits, and check that each one can be explained by a hardware hit.
  for (size_t i = 0; i < softwareSsvHits.size();i++) {
    //get the hit data, and create the positions that would correspond to the hw hit
    size_t softwarePhmmIndex = softwareSsvHits[i].phmmPosition;
    size_t sequenceIndex = softwareSsvHits[i].sequencePosition;
    size_t sequenceSegmentIndex = sequenceIndex / NUM_CELL_PROCESSORS;
    size_t sequencePositionInSegment = sequenceIndex - (sequenceSegmentIndex * NUM_CELL_PROCESSORS);
    size_t softwareGroupIndex = sequencePositionInSegment / CELLS_PER_GROUP;

    //iterate through the hardware hits, and check if any of them explain the hit.
    bool foundMatchingHitReport = false;
    for (auto& hardwareHit : hardwareSsvHits) {
      //check to see if the phmm index and sequence pass match
      if ((hardwareHit.phmmIndex == softwarePhmmIndex) && (hardwareHit.sequenceIndex == sequenceSegmentIndex)) {
        //check to see if the group agrees, too
        if (hardwareHit.groupsPassingThreshold[softwareGroupIndex]) {
          //the hit checks out, we're good!
          foundMatchingHitReport = true;
          break;
        }
        else {
          printf("TEST FAIL: software hit @ p=%zu, s=%zu had a corresponding hardware hit, but the group %zu did not represent a hit\n",
            softwarePhmmIndex, sequenceIndex, softwareGroupIndex);
          allTestsPass = false;
          break;
        }
        break;
      }
    }
    if (!foundMatchingHitReport) {
      printf("TEST FAIL: software hit @ p=%zu, s=%zu did not have a corresponding hardware hit report\n", softwarePhmmIndex, sequenceIndex);
      allTestsPass = false;
    }
  }

  //now, we can iterate over the hardware hits, and see if they're all represented 
  for (auto& hardwareHit : hardwareSsvHits) {
    //assertion check to make sure the hardware hit actually has reported hit groups.
    if (hardwareHit.groupsPassingThreshold.or_reduce() == 0) {
      printf("TEST FAILURE: hardware hit report @ s=%zu p=%zu had no groups asserted.\n", hardwareHit.sequenceIndex, hardwareHit.phmmIndex);
      allTestsPass = false;
    }


    //any group in this hitreport that asserted a hit should check for a corresponding software hit.
    for (size_t groupIndex = 0; groupIndex < NUM_CELL_GROUPS; groupIndex++) {
      if (hardwareHit.groupsPassingThreshold[groupIndex]) {

        bool foundMatchingSoftwareHit = false;
        for (auto& softwareHit : softwareSsvHits) {
          size_t sequenceIndex = softwareHit.sequencePosition;
          size_t sequenceSegmentIndex = sequenceIndex / NUM_CELL_PROCESSORS;
          size_t sequencePositionInSegment = sequenceIndex - (sequenceSegmentIndex * NUM_CELL_PROCESSORS);
          size_t softwareGroupIndex = sequencePositionInSegment / CELLS_PER_GROUP;

          if ((softwareHit.phmmPosition == hardwareHit.phmmIndex) && (sequenceSegmentIndex == hardwareHit.sequenceIndex) &&
            (softwareGroupIndex == groupIndex)) {
            foundMatchingSoftwareHit = true;
            break;
          }
        }

        if (!foundMatchingSoftwareHit) {
          printf("TEST FAILURE: hardware hit report s=%zu p=%zu, group %zu was not represented in the software hit list.\n",
            hardwareHit.sequenceIndex, hardwareHit.phmmIndex, groupIndex);
          allTestsPass = false;
        }
      }
    }
  }

  return allTestsPass;
}

// //cast the phmm's match emission scores to int8_t so that HAVAC can work with them.
// int8_t* makeEmissionsAsInt8_tList(struct HmmSeqPair& hmmSeqPair) {
//   int8_t* emissionsAsInt8_t = (int8_t*)malloc(hmmSeqPair.sequenceLength * 4 * sizeof(int8_t));
//   if (emissionsAsInt8_t == NULL) {
//     printf("could not allocate memory for emission array as int8_t array\n");
//     exit(6);
//   }
//   for (size_t i = 0; i < hmmSeqPair.sequenceLength * 4; i++) {
//     emissionsAsInt8_t[i] = (int8_t)hmmSeqPair.phmmList->phmms[0].model.matchEmissionScores[i];
//   }

//   return emissionsAsInt8_t;
// }


//converts the sequence from ascii characters to integer indexes (a=0, c=1, g=2, c=3).
uint8_t* sequenceAsUnpackedVectorIndices(struct HmmSeqPair& hmmSeqPair) {
	uint32_t allocSize = hmmSeqPair.sequenceLength * sizeof(uint8_t);
	uint8_t* sequenceAsVectorIndices =  new uint8_t[hmmSeqPair.sequenceLength];
	if(!sequenceAsVectorIndices){
		printf("ERROR: could not allocate memory for sequenceAsVectorIndices array. requested size\n" );
		fflush(stdout);
		exit(7);
	}
  for (uint32_t i = 0; i < hmmSeqPair.sequenceLength; i++) {
    sequenceAsVectorIndices[i] = nucToEncoding(hmmSeqPair.sequence[i]);
  }
  return sequenceAsVectorIndices;
}

//packs the integer index sequence from sequenceAsUnpackedVectorIndices() into 2-bit values, 4 per byte
//this is the required format for the hardware HAVAC implementation
void packSequenceVectorIndices(uint8_t* sequenceAsVectorIndices, size_t sequenceLength) {
  for (size_t i = 0; i < sequenceLength; i++) {
    uint8_t symbolEncoding = sequenceAsVectorIndices[i];
    uint_fast8_t bitShiftAmount = (i % 4) * 2;
    size_t byteToSetEncoding = i / 4;
    uint8_t bitmask = 0x3 << bitShiftAmount;

    //bitmask the previous bits away from where we want to write
    sequenceAsVectorIndices[byteToSetEncoding] &= ~bitmask;
    //write the new bits there
    sequenceAsVectorIndices[byteToSetEncoding] |= (symbolEncoding << bitShiftAmount);
  }
}

uint8_t* generateCompressedSequence(const char* sequence) {
  const uint32_t seqLength = strlen(sequence);
  const uint32_t compressedSeqByteLength = (seqLength + 3) / 4; // rounded up
  uint8_t* compressedSeq = (uint8_t*)calloc(compressedSeqByteLength, 1);

  for (size_t i = 0; i < seqLength; i++) {
    uint32_t byteIndex = i / 4;
    uint8_t bitPosition = (i % 4) * 2;
    compressedSeq[byteIndex] |= (nucToEncoding(sequence[i])) << bitPosition;
  }

  return compressedSeq;
}

uint8_t nucToEncoding(const uint8_t nuc) {
  switch (nuc) {
  case 'a':
    return 0;
  case 'c':
    return 1;
  case 'g':
    return 2;
  default:
    return 3;
  }
}


void generateTestSet(const uint32_t numTests, const uint32_t sequenceLengthInSegments, const float sequenceMutProbability){
	const float sequenceMutationProbability = .3f;
	for(uint32_t testIndex = 1; testIndex <= numTests; testIndex++){
		char phmmSrc[256];
		char seqSrc[256];
		char mutSeqSrc[256];

		makePhmmFileSrcFromTestNum(testIndex, phmmSrc);
		makeSequenceFileSrcFromTestNum(testIndex, seqSrc);
		makeMutatedSequenceFileSrcFromTestNum(testIndex, mutSeqSrc);

	    const uint32_t seqLength = NUM_CELL_PROCESSORS * sequenceLengthInSegments;
		generateSequencesToFile( seqLength, seqSrc, mutSeqSrc, sequenceMutationProbability);
		generatePhmmToFile(seqSrc, phmmSrc);
	}
}

struct HmmSeqPair readTestSetFromFiles(const uint32_t testNum){
	char phmmSrc[128];
	char mutSeqSrc[128];
	makePhmmFileSrcFromTestNum(testNum,phmmSrc);
	makeMutatedSequenceFileSrcFromTestNum(testNum, mutSeqSrc);

	//read the fastaSequence
	struct FastaVector fastaVector;
	enum FastaVectorReturnCode returnCode = fastaVectorInit(&fastaVector);
	if(returnCode != FASTA_VECTOR_OK){
		printf("ERROR: could not create fasta vector, returned error code %d\n", returnCode);
		exit(300);
	}

	returnCode = fastaVectorReadFasta(mutSeqSrc, &fastaVector);
	if(returnCode != FASTA_VECTOR_OK){
		printf("ERROR: could not open fasta file %s, returned error code %d\n", mutSeqSrc , returnCode);
		exit(301);
	}

	struct HmmSeqPair hmmSeqPair;
	hmmSeqPair.phmmList = readHmm(phmmSrc);
	const uint32_t modelLength = hmmSeqPair.phmmList->phmms[0].header.modelLength;
	hmmSeqPair.sequenceLength = modelLength;
	hmmSeqPair.sequence = (char*)malloc(modelLength * sizeof(char));
	memcpy(hmmSeqPair.sequence, fastaVector.sequence.charData, modelLength);
	fastaVectorDealloc(&fastaVector);

	return hmmSeqPair;

}
void makePhmmFileSrcFromTestNum(const uint32_t testNum, char *buffer){
	sprintf(buffer, "/home/ta112817/tmp/test_%d.phmm", testNum);

}
void makeSequenceFileSrcFromTestNum(const uint32_t testNum, char *buffer){
	sprintf(buffer, "/home/ta112817/tmp/test_%d.fasta", testNum);
}

void makeMutatedSequenceFileSrcFromTestNum(const uint32_t testNum, char *buffer){
	sprintf(buffer, "/home/ta112817/tmp/testMut_%d.fasta", testNum);
}

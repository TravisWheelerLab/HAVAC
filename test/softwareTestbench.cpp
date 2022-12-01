#include "softwareTestBench.h"
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

//defined in PublicDefines.h, if defined at all
#ifdef HAVAC_PER_CELL_DATA_TESTING
#include "byCellComparator/byCellComparator.hpp"
#endif



/* This file contains the testbench that ensure that HAVAC will return the same results as a traditional software SSV implementation.
*/
uint8_t* sequenceAsUnpackedVectorIndices(struct HmmSeqPair& hmmSeqPair);
void packSequenceVectorIndices(uint8_t* sequenceAsVectorIndices, size_t sequenceLength);
int8_t* makeEmissionsAsInt8_tList(struct HmmSeqPair& hmmSeqPair);
uint8_t nucToEncoding(const uint8_t nuc);

const uint32_t numTests = 4;
const float requiredPValue = 0.05f;

int main() {
  srand(time(0));
  printf("BEGINNING TESTBENCH...\n");
  bool allTestsPassed = true;

  for (size_t i = 1; i <= numTests; i++) {
    #ifdef HAVAC_PER_CELL_DATA_TESTING
    clearCellComparatorMaps();
    #endif
    //generate the sequence and phmm for this test. The sequence will be mutated with a substitution rate specified
    //in the hmmSeqGenerator file
    printf(" Test #%u/%u \n", i, numTests);
    const uint32_t numSequenceSegments = (rand() % 3) + 1;
    const uint32_t seqLength = NUM_CELL_PROCESSORS * numSequenceSegments;
    printf("generating hmm seq pair...\n");
    fflush(stdout);
    struct HmmSeqPair hmmSeqPair = generateRandomHmmSeqPair(seqLength);

    printf("reprojecting scores...\n");
    fflush(stdout);
    // reproject the scores
    //the hmm/seq generator makes a phmm with only one model in it, so we can use phmms[0] to reference it in this function.
    p7HmmProjectForThreshold256(&hmmSeqPair.phmmList->phmms[0], requiredPValue);
    int8_t* emissionsAsInt8_t = makeEmissionsAsInt8_tList(hmmSeqPair);
    uint8_t* sequenceAsVectorIndices = sequenceAsUnpackedVectorIndices(hmmSeqPair);

    printf("invoking software ssv...\n");
    fflush(stdout);
    int8_t errorCode = 0;
    std::vector<struct SoftSsvHit> softwareSsvHits = softSsvThreshold256(sequenceAsVectorIndices,
      hmmSeqPair.sequenceLength, emissionsAsInt8_t, hmmSeqPair.phmmList->phmms[0].header.modelLength, errorCode);
    if (errorCode != 0) {
      printf("Error: soft ssv returned error code %i\n", errorCode);
      exit(8);
    }
    //pack the sequence into 2-bit values for the hardware HAVAC implementation.
    packSequenceVectorIndices(sequenceAsVectorIndices, hmmSeqPair.sequenceLength);


    printf("invoking hardware ssv...\n");
    fflush(stdout);
#ifdef USE_HIT_SIEVE
    std::vector<struct HitReportByGroup> hardwareSsvHits = invokeHardwareSsv(sequenceAsVectorIndices, hmmSeqPair.sequenceLength,
      emissionsAsInt8_t, hmmSeqPair.phmmList->phmms[0].header.modelLength, errorCode);
#else
    std::vector<struct HitReport> hardwareSsvHits = invokeHardwareSsv(sequenceAsVectorIndices, hmmSeqPair.sequenceLength,
      emissionsAsInt8_t, hmmSeqPair.phmmList->phmms[0].header.modelLength, errorCode);
#endif


    printf("soft ssv returned %zu hits\n", softwareSsvHits.size());
    printf("hard ssv returned %zu hits\n", hardwareSsvHits.size());
    fflush(stdout);

    #ifdef HAVAC_PER_CELL_DATA_TESTING
    printf("comparing data on per-cell basis\n");
    compareCellMaps(hmmSeqPair.phmmList->phmms[0].header.modelLength, hmmSeqPair.sequenceLength, hmmSeqPair);
    #endif

    printf("comparing hit lists...\n");
    fflush(stdout);
    bool comparesSuccessfully = compareSsvHitLists(hardwareSsvHits, softwareSsvHits);
    allTestsPassed &= comparesSuccessfully;
    if (!comparesSuccessfully) {
      printf("Test Fail: test number %zu did not compare successfully.\n", i);
    }

    p7HmmListDealloc(hmmSeqPair.phmmList);
    free(emissionsAsInt8_t);
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
#ifdef USE_HIT_SIEVE
inline std::vector<struct HitReportByGroup> invokeHardwareSsv(const uint8_t* sequenceAsVectorIndices, size_t sequenceLength,
  const int8_t* emissionsAsInt8_t, const size_t phmmLength, int8_t& errorCode) {
#else
inline std::vector<struct HitReport> invokeHardwareSsv(const uint8_t* sequenceAsVectorIndices, size_t sequenceLength,
  const int8_t* emissionsAsInt8_t, const size_t phmmLength, int8_t& errorCode) {
#endif
	errorCode = 0;
#ifdef USE_HIT_SIEVE
	std::vector<struct HitReportByGroup> hitReportList;
#else
  std::vector<struct HitReport> hitReportList;
#endif
  const size_t sequenceLengthInSegments = sequenceLength / (NUM_CELL_PROCESSORS);
  if (sequenceLength % NUM_CELL_PROCESSORS != 0) {
    printf("ERROR: seq length was %zu was not a multiple of num_cell_processors %zu!\n", sequenceLength, NUM_CELL_PROCESSORS);
    fflush(stdout);
    errorCode = -1;
  }


  const uint32_t NUM_HIT_REPORTS_SUPPORTED = phmmLength * sequenceLengthInSegments;
#ifdef USE_HIT_SIEVE
  struct HitReportByGroup* hitReportMemory = (struct HitReportByGroup*)malloc(sizeof(struct HitReportByGroup) * NUM_HIT_REPORTS_SUPPORTED);
#else
  struct HitReport* hitReportMemory = (struct HitReport*)malloc(sizeof(struct HitReport) * NUM_HIT_REPORTS_SUPPORTED);
#endif
  if (hitReportMemory == NULL) {
    printf("ERROR: could not allocate memory for hit report buffer\n");
    errorCode = -2;
    return hitReportList;
  }
  printf("invoking havac kernel\n");
  fflush(stdout);
  uint32_t numHits = 0;
  HavacKernelTopLevel((struct SequenceSegment*)sequenceAsVectorIndices, sequenceLengthInSegments,
    (struct PhmmVector*)emissionsAsInt8_t, phmmLength, hitReportMemory, numHits);
  printf("havac kernel finished\n");
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
#ifdef USE_HIT_SIEVE
bool compareSsvHitLists(std::vector<struct HitReportByGroup> hardwareSsvHits, std::vector<struct SoftSsvHit> softwareSsvHits) {
#else
bool compareSsvHitLists(std::vector<struct HitReport> hardwareSsvHits, std::vector<struct SoftSsvHit> softwareSsvHits) {
#endif

  bool allTestsPass = true;
  //iterate over the software ssv hits, and check that each one can be explained by a hardware hit.
  for (size_t i = 0; i < softwareSsvHits.size();i++) {
    //get the hit data, and create the positions that would correspond to the hw hit
    size_t softwarePhmmIndex = softwareSsvHits[i].phmmPosition;
    size_t sequenceIndex = softwareSsvHits[i].sequencePosition;
    size_t sequenceSegmentIndex = sequenceIndex / NUM_CELL_PROCESSORS;
    size_t sequencePositionInSegment = sequenceIndex - (sequenceSegmentIndex * NUM_CELL_PROCESSORS);
    size_t softwareGroupIndex = sequencePositionInSegment / CELLS_PER_GROUP;

#ifdef USE_HIT_SIEVE
    size_t expectedGroupIndex = sequencePositionInSegment / HIT_REPORT_BY_GROUP_MIN_BIT_WIDTH;
    size_t expectedBitInGroup = sequencePositionInSegment % HIT_REPORT_BY_GROUP_MIN_BIT_WIDTH;
#endif

    //iterate through the hardware hits, and check if any of them explain the hit.
    bool foundMatchingHitReport = false;
    for (auto& hardwareHit : hardwareSsvHits) {
      //check to see if the phmm index and sequence pass match
#ifdef USE_HIT_SIEVE
      if ((hardwareHit.phmmIndex == softwarePhmmIndex) && (hardwareHit.sequenceIndex == sequenceSegmentIndex) && (hardwareHit.groupIndex == expectedGroupIndex)) {
#else
      if ((hardwareHit.phmmIndex == softwarePhmmIndex) && (hardwareHit.sequenceIndex == sequenceSegmentIndex)) {

#endif
    	  //check to see if the group agrees, too
#ifdef USE_HIT_SIEVE
        if (hardwareHit.groupBits[expectedBitInGroup]) {
#else
        if (hardwareHit.groupsPassingThreshold[softwareGroupIndex]) {
#endif

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
#ifdef USE_HIT_SIEVE
	    if (hardwareHit.groupBits.or_reduce() == 0) {
#else
    if (hardwareHit.groupsPassingThreshold.or_reduce() == 0) {
#endif
      printf("TEST FAILURE: hardware hit report @ s=%zu p=%zu had no groups asserted.\n", hardwareHit.sequenceIndex, hardwareHit.phmmIndex);
      allTestsPass = false;
    }

#ifdef USE_HIT_SIEVE
    for(uint32_t cellIndexInGroup = 0; cellIndexInGroup < HIT_REPORT_BY_GROUP_MIN_BIT_WIDTH; cellIndexInGroup++){
    	if(hardwareHit.groupBits[cellIndexInGroup]){
    		size_t sequencePosition = (hardwareHit.groupIndex * HIT_REPORT_BY_GROUP_MIN_BIT_WIDTH) + (hardwareHit.sequenceIndex * NUM_CELL_PROCESSORS) + cellIndexInGroup;
    		bool foundMatchingSoftwareHit = false;
    		for(auto& softwareHit: softwareSsvHits){
    			if(softwareHit.phmmPosition = hardwareHit.phmmIndex && softwareHit.sequencePosition == sequencePosition){
    				foundMatchingSoftwareHit = true;
    				break;
    			}
    		}
    		if(!foundMatchingSoftwareHit){
    	          printf("TEST FAILURE: hardware hit report s=%u p=%u, bit %u was not represented in the software hit list.\n",
    	            hardwareHit.sequenceIndex, hardwareHit.phmmIndex, cellIndexInGroup);
    	          allTestsPass = false;
    		}
    	}
    }
#else
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
#endif
  }

  return allTestsPass;
}

//cast the phmm's match emission scores to int8_t so that HAVAC can work with them.
int8_t* makeEmissionsAsInt8_tList(struct HmmSeqPair& hmmSeqPair) {
  int8_t* emissionsAsInt8_t = (int8_t*)malloc(hmmSeqPair.sequenceLength * 4 * sizeof(int8_t));
  if (emissionsAsInt8_t == NULL) {
    printf("could not allocate memory for emission array as int8_t array\n");
    exit(6);
  }
  for (size_t i = 0; i < hmmSeqPair.sequenceLength * 4; i++) {
    emissionsAsInt8_t[i] = (int8_t)hmmSeqPair.phmmList->phmms[0].model.matchEmissionScores[i];
  }

  return emissionsAsInt8_t;
}


//converts the sequence from ascii characters to integer indexes (a=0, c=1, g=2, c=3).
uint8_t* sequenceAsUnpackedVectorIndices(struct HmmSeqPair& hmmSeqPair) {
  uint8_t* sequenceAsVectorIndices = (uint8_t*)malloc(hmmSeqPair.sequenceLength * sizeof(uint8_t));
  if (sequenceAsVectorIndices == NULL) {
    printf("ERROR: could not allocate memory for sequenceAsVectorIndices array\n");
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

// std::vector<struct HitReport> softSsv(const char* sequence, const uint32_t sequenceLength,
//   const uint8_t* phmmMatchScores, const uint32_t phmmLengthInVectors,
//   const uint32_t cellsPerGroup, const uint32_t numGroups) {

//   const uint32_t numCellProcessors = cellsPerGroup / numGroups;
//   std::vector<struct HitReport> hitReportList;
//   uint8_t* segmentScores = (uint8_t*)malloc(numCellProcessors);
//   if (segmentScores == NULL) {
//     printf("ERROR: could not allocate memory for the soft ssv segment scores\n");
//     return hitReportList;
//   }

//   uint8_t* prevScoreColumn = (uint8_t*)malloc(phmmLengthInVectors);
//   if (prevScoreColumn == NULL) {
//     printf("ERROR: could not allocate memory for the soft ssv prev score column\n");
//     return hitReportList;
//   }

//   const uint32_t numSequenceSegments = seqLength / numCellProcessors;
//   for (size_t sequenceSegmentIndex = 0; sequenceSegmentIndex < numSequenceSegments; sequenceSegmentIndex++) {
//     memset(segmentScores, 0, numCellProcessors);

//     const uint32_t baseCellOffset = sequenceSegmentIndex * numCellProcessors;
//     for (size_t phmmIndex = 0; phmmIndex < phmmLengthInVectors; phmmIndex++) {
//       uint64_t groupedHits = 0;
//       for (int64_t cellInSegment = numCellProcessors - 1; cellInSegment > 0; cellInSegment--) {
//         uint8_t prevScore;
//         if (phmmIndex == 0) {
//           prevScore = 0;
//         }
//         else {
//           prevScore = cellInSegment == 0 ? prevScoreColumn[phmmIndex - 1] : segmentScores[phmmIndex - 1];
//         }
//         uint8_t sequenceSymbolOffset;
//         char sequenceSymbol = sequence[baseCellOffset + cellInSegment];
//         switch (sequenceSymbol) {
//         case 'a':
//           sequenceSymbolOffset = 0;
//           break;
//         case 'c':
//           sequenceSymbolOffset = 1;
//           break;
//         case 'g':
//           sequenceSymbolOffset = 2;
//           break;
//         default:
//           sequenceSymbolOffset = 3;
//         }

//         uint8_t matchScore = phmmMatchScores[phmmIndex * 4 + sequenceSymbolOffset];
//         int32_t prevScoreAs32 = prevScore;
//         int32_t matchScoreAs32 = matchScore;
//         int32_t currentScore = prevScoreAs32 + matchScoreAs32;

//         if (currentScore < 0) {
//           currentScore = 0;
//         }
//         else if (currentScore >= 256) {
//           currentScore = 0;
//           groupedHits |= 1 << cellInSegment / cellsPerGroup;
//         }
//       }

//       if (groupedHits != 0) {
//         // create and puch_back the hit
//         struct HitReport hitReport;
//         hitReport.phmmIndex = phmmIndex;
//         hitReport.sequenceIndex = sequenceSegmentIndex;
//         hitReport.groupedHits = groupedHits;

//         hitReportList.push_back(hitReport);
//       }
//     }
//   }
//   return hitReportList;
// }

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

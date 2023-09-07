#include "Havac.hpp"
#include "types/HavacHit.hpp"
#include "../PhmmReprojection/PhmmReprojection.h"
#include "phmm/PhmmPreprocessor.hpp"
#include "sequence/SequencePreprocessor.hpp"

#include <stdexcept>
#include <exception>
#include <array>
#include <chrono>
#include <tuple>

struct PhmmLocalPosition {
  uint32_t phmmIndex;
  uint32_t phmmPosition;
};

PhmmLocalPosition phmmPrefixSumsBinarySearch(uint32_t phmmGlobalPosition, vector<uint32_t>& prefixSums);

Havac::Havac(const uint32_t deviceIndex, const float requiredPValue, const std::string xclbinSrc)
  :deviceIndex(deviceIndex),
  requiredPValue(requiredPValue),
  havacXclbinFileSrc(xclbinSrc) {
  this->hwClient = std::make_shared<HavacHwClient>(this->havacXclbinFileSrc, this->havacKernelName, deviceIndex);
  this->fastaVector = (FastaVector*)malloc(sizeof(FastaVector));
  this->p7HmmList = (P7HmmList*)malloc(sizeof(P7HmmList));
  enum FastaVectorReturnCode rc = fastaVectorInit(this->fastaVector);
  if (rc == FASTA_VECTOR_ALLOCATION_FAIL) {
    throw std::bad_alloc();
  }
}

Havac::~Havac() {
  fastaVectorDealloc(this->fastaVector);
  p7HmmListDealloc(this->p7HmmList);
  free(this->fastaVector);
  free(this->p7HmmList);

}


void Havac::loadPhmm(const std::string phmmSrc) {
  enum P7HmmReturnCode rc = readP7Hmm(phmmSrc.c_str(), this->p7HmmList);
  if(rc == p7HmmFileNotFound)
  if (rc == p7HmmAllocationFailure) {
    throw std::bad_alloc();
  }
  else if (rc == p7HmmFormatError) {
    throw std::runtime_error("Phmm file was not formatted correctly.");
  }
  shared_ptr<PhmmPreprocessor> phmmPreprocessor = std::make_shared<PhmmPreprocessor>(this->p7HmmList, this->requiredPValue);
  this->compressedPhmmScores = phmmPreprocessor->getProcessedPhmmData();
  this->hwClient->writePhmm(compressedPhmmScores);
  this->phmmLoadedToDevice = true;
}

void Havac::loadSequence(const std::string fastaSrc) {
  enum FastaVectorReturnCode rc = fastaVectorReadFasta(fastaSrc.c_str(), fastaVector);
  if (rc == FASTA_VECTOR_ALLOCATION_FAIL) {
    throw std::bad_alloc();
  }
  else if (rc == FASTA_VECTOR_FILE_OPEN_FAIL) {
    throw std::runtime_error("Could not open fasta file for reading.");
  }
  else if (rc == FASTA_VECTOR_FILE_READ_FAIL) {
    throw std::runtime_error("Error while reading from the opened fasta file.");
  }

  shared_ptr<SequencePreprocessor> sequencePreprocessor =
    std::make_shared<SequencePreprocessor>(this->fastaVector);

  vector<uint8_t>& compressedSeq = sequencePreprocessor->getCompressedSequenceBuffer();

  size_t seqSize = compressedSeq.size();
  hwClient->writeSequence(compressedSeq);
  this->sequenceLoadedToDevice = true;
}


void Havac::runHardwareClient() {
  this->runHardwareClientAsync();
  this->waitHardwareClientAsync();
}

void Havac::runHardwareClientAsync() {
  if (!this->phmmLoadedToDevice) {
    throw std::logic_error("Phmm was not loaded to device before hardware was requested to run.");
  }
  if (!this->sequenceLoadedToDevice) {
    throw std::logic_error("Sequence was not loaded to device before hardware was requested to run.");
  }

  this->hwClient->invokeHavacSsvAsync();
}

void Havac::waitHardwareClientAsync() {
  this->hwClient->waitForHavacSsvAsync();
}

void Havac::abortHardwareClient() {
  this->hwClient->abort();
}

vector<uint32_t> Havac::generatePhmmLenPrefixSums(){
  vector<uint32_t> phmmPrefixSums;
  phmmPrefixSums.reserve(this->p7HmmList->count+1);
  phmmPrefixSums.push_back(0);
  uint32_t modelStartingPosition = 0;

  for(uint32_t i = 0; i < this->p7HmmList->count; i++){
    modelStartingPosition += this->p7HmmList->phmms[i].header.modelLength;
    phmmPrefixSums.push_back(modelStartingPosition);
  }

  return phmmPrefixSums;
}

//finds the lowest prefix sum that's greater than the global position
PhmmLocalPosition phmmPrefixSumsBinarySearch(uint32_t phmmGlobalPosition, vector<uint32_t> &prefixSums){
  int32_t leftBound = 0;
  int32_t rightBound = prefixSums.size() - 1;
  int32_t result = -1; // Initialize result to an invalid value

  while (leftBound <= rightBound) {
    uint32_t pivotIndex = leftBound + (rightBound - leftBound) / 2;
    if (prefixSums[pivotIndex] <= phmmGlobalPosition) {
      result = pivotIndex;
      leftBound = pivotIndex + 1;
    }
    else {
      rightBound = pivotIndex - 1;
    }
  }

  PhmmLocalPosition localPosition;
  localPosition.phmmIndex = result;
  if(result != -1){
    fflush(stdout);
    localPosition.phmmPosition = phmmGlobalPosition - prefixSums[result];
  }
  return localPosition;
}


vector<HavacHit> Havac::getHitsFromFinishedRun() {
  vector<uint64_t> rawHitsAsU64 = this->hwClient->getHitList();
  vector<uint32_t> phmmPrefixSums = generatePhmmLenPrefixSums();

  vector<HavacHit> resolvedHavacHits;
  resolvedHavacHits.reserve(rawHitsAsU64.size());
  //try to resolve the HavacHit
  for(uint64_t i = 0; i < rawHitsAsU64.size(); i++){
    
    uint64_t hitAsU64 = rawHitsAsU64[i];
    const uint64_t sequencePartitionBitmaskEndBit = 14;
    const uint64_t sequenceSegmentBitmaskEndBit = 40;
    uint64_t sequencePartitionBitMask = (1ULL << sequencePartitionBitmaskEndBit) - 1;
    uint64_t sequenceSegmentBitmask = (1ULL << sequenceSegmentBitmaskEndBit) - 1;
    uint64_t sequencePartitionIndex = hitAsU64 & sequencePartitionBitMask;
    uint64_t sequenceSegmentIndex = (hitAsU64 & sequenceSegmentBitmask) >> 14;
    uint64_t globalSequencePosition = (hitAsU64 & sequenceSegmentBitmask) >> 14;
    globalSequencePosition = (sequenceSegmentIndex * (12 * 1024)) + sequencePartitionIndex;
    uint32_t globalPhmmPosition = hitAsU64 >> sequenceSegmentBitmaskEndBit;

    struct FastaVectorLocalPosition localSequencePosition;
    bool fastaVectorSuccess = fastaVectorGetLocalSequencePositionFromGlobal(
      fastaVector, globalSequencePosition, &localSequencePosition);

    if (!fastaVectorSuccess) {
      //the sequence position was outside the bounds of the fasta sequences, 
      //so it was a hit in the padding data, we can ignore this hit
      continue;
    }

    PhmmLocalPosition phmmLocalPosition = phmmPrefixSumsBinarySearch(globalPhmmPosition, phmmPrefixSums);

    if (phmmLocalPosition.phmmIndex == -1) {
      std::cerr<< "ERROR: could not resolve phmm position for raw hit report #"<< i << "\n" << std::endl;
      continue;
    }

    resolvedHavacHits.push_back(HavacHit(localSequencePosition.positionInSequence, localSequencePosition.sequenceIndex,
      phmmLocalPosition.phmmPosition, phmmLocalPosition.phmmIndex));
  }

  return resolvedHavacHits;
}


enum havac_cmd_state Havac::currentHardwareState() {
  return (havac_cmd_state) hwClient->getHwState();
}

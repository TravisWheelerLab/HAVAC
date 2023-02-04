#include "Havac.hpp"
#include "../PhmmReprojection/PhmmReprojection.h"
#include <exception>
#include <stdexcept>
#include "exceptions/FileReadException.hpp"

Havac::Havac(const uint32_t deviceIndex = 0)
  :deviceIndex(deviceIndex) {
  this->hwClient = std::shared_ptr<HavacHwClient>(new HavacHwClient(this->havacXclbinFileSrc, this->havacKernelName, deviceIndex));
  this->fastaVector = new FastaVector;
  this->p7HmmList = new P7HmmList;

  enum FastaVectorReturnCode rc = fastaVectorInit(&this->fastaVector);
  if (rc == FASTA_VECTOR_INIT_FAIL) {
    throw std::bad_alloc();
  }
}
Havac::~Havac() {
  fastaVectorDealloc(this->fastaVector);
  p7HmmListDealloc(this->p7HmmList);
  delete (this->fastaVector);
  delete (this->p7HmmList);

}


void Havac::loadPhmm(const std::string phmmSrc, const float desiredPValue) {

  enum P7HmmReturnCode rc = readP7Hmm(phmmSrc.c_str(), this->p7HmmList);
  if (rc == p7HmmAllocationFailure) {
    throw std::bad_alloc();
  }
  else if (rc == p7HmmFormatError) {
    throw std::runtime_error("Phmm file was not formatted correctly.");
  }

  std::shared_ptr<PhmmPreprocessor> phmmPreprocessor = 
    std::shared_ptr<PhmmPreprocessor>(new PhmmPreprocessor(this->p7HmmList, desiredPValue));
  const int8_t* phmmDataPtr = phmmPreprocessor->getProcessedPhmmData().get();
  const uint32_t phmmLengthInBytes = phmmPreprocessor->getPhmmLengthInBytes();

  this->hwClient->writePhmm(phmmDataPtr, phmmLengthInBytes);
  this->phmmLoadedToDevice = true;
}

void Havac::loadSequence(const std::string fastaSrc) {
  enum FastaVectorReturnCode rc = fastaVectorReadFasta(fastaSrc, this->fastaVector);
  if (rc == FASTA_VECTOR_ALLOCATION_FAIL) {
    throw std::bad_alloc();
  }
  else if (rc == FASTA_VECTOR_FILE_OPEN_FAIL) {
    throw FileReadException("Could not open fasta file for reading.");
  }
  else if (rc == FASTA_VECTOR_FILE_READ_FAIL) {
    throw FileReadException("Error while reading from the opened fasta file.");
  }

  std::shared_ptr<SequencePreprocessor> sequencePreprocessor = 
    std::shared_ptr<SequencePreprocessor>(new SequencePreprocessor(this->fastaVector));


  uint8_t* sequenceAsBytes = sequencePreprocessor->getCompressedSequenceBuffer.get();
  uint32_t sequenceLengthInBytes = sequencePreprocessor->getCompressedSquenceLengthInBytes()
  hwClient->writeSequence(sequenceAsBytes, sequenceLengthInBytes);

  this->sequenceLoadedToDevice = true;

}


void Havac::runHardwareClient(){
  this->runHardwareClientAsync();
  this->waitHardwareClientAsync();
}

void Havac::runHardwareClientAsync() {
  if(!this->phmmLoadedToDevice){
    throw std::logic_error("Phmm was not loaded to device before hardware was requested to run.");
  }
  if (!this->sequenceLoadedToDevice) {
    throw std::logic_error("Sequence was not loaded to device before hardware was requested to run.");
  }
  this->hwClient->invokeHavacSsvAsync();
}

void Havac::waitHardwareClientAsync(){
  this->hwClient->waitForHavacSsvAsync();
}

void Havac::abortHardwareClient(){
  this->hwClient->abort();
}

std::shared_ptr<std::vector<HavacHitLocation>> getHitsFromFinishedRun(){
  std::shared_ptr<HitVerifier> hitVerifier = 
  std::shared_ptr<HitVerifier>(new HitVerifier(this->fastaVector, this->p7HmmmList));
  hitVerifer->verify();
}

bool Havac::isRunning() {
  //not yet implemented
}


// void Havac::generateProjectedPhmmMatchScores(const float pValue = 0.05f) {
//   //determine how long the projectedScoresArray needs to be.
//   uint32_t projectedScoresLength = 0;
//   for (uint32_t phmmIndex = 0; phmmIndex < this->p7HmmList->count; phmmIndex++) {
//     projectedScoresLength += this->phmmList->phmms[phmmIndex].header.modelLength * 4;
//   }

//   try {
//     projectedPhmmMatchScores = make_shared<int8_t[]>(projectedScoresLength);
//   }
//   (const std::bad_alloc & e) {
//     std::cerr << "ERROR: could not allocate memory for reprojected phmm scores.\n" << e.what() << std::endl;
//   }
//   uint32_t phmmOffset = 0
//     for (uint32_t phmmIndex = 0; phmmIndex < this->phmmList->count; phmmIndex++) {
//       int8_t* projectedPhmmScorePtr = &projectedPhmmMatchScores[phmmOffset * 4];
//       p7HmmProjectForThreshold256(&this->phmmList->phmms[phmmIndex], pValue, projectedPhmmScorePtr);
//       phmmOffset += this->phmmList->phmms[phmmIndex].header.modelLength;
//     }
// }
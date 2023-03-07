#include "Havac.hpp"
#include "../PhmmReprojection/PhmmReprojection.h"
#include "exceptions/FileReadException.hpp"
#include "phmm/PhmmPreprocessor.hpp"
#include "sequence/SequencePreprocessor.hpp"

#include <stdexcept>
#include <exception>

Havac::Havac(const uint32_t deviceIndex)
  :deviceIndex(deviceIndex) {
  this->hwClient = shared_ptr<HavacHwClient>(new HavacHwClient(this->havacXclbinFileSrc, this->havacKernelName, deviceIndex));
  this->fastaVector = shared_ptr<FastaVector>(new FastaVector);
  this->p7HmmList = shared_ptr<P7HmmList>(new P7HmmList);
  //here, 
  enum FastaVectorReturnCode rc = fastaVectorInit(this->fastaVector.get());
  if (rc == FASTA_VECTOR_ALLOCATION_FAIL) {
    throw std::bad_alloc();
  }
}
Havac::~Havac() {
  fastaVectorDealloc(this->fastaVector.get());
  p7HmmListDealloc(this->p7HmmList.get());

}


void Havac::loadPhmm(const std::string phmmSrc, const float desiredPValue) {

  enum P7HmmReturnCode rc = readP7Hmm(phmmSrc.c_str(), this->p7HmmList.get());
  if (rc == p7HmmAllocationFailure) {
    throw std::bad_alloc();
  }
  else if (rc == p7HmmFormatError) {
    throw std::runtime_error("Phmm file was not formatted correctly.");
  }

  shared_ptr<PhmmPreprocessor> phmmPreprocessor = std::make_shared<PhmmPreprocessor>(this->p7HmmList, desiredPValue);

  shared_ptr<vector<int8_t>> phmmDataPtr = phmmPreprocessor->getProcessedPhmmData();

  this->hwClient->writePhmm(phmmDataPtr);
  this->phmmLoadedToDevice = true;
}

void Havac::loadSequence(const std::string fastaSrc) {
  enum FastaVectorReturnCode rc = fastaVectorReadFasta(fastaSrc.c_str(), this->fastaVector.get());
  if (rc == FASTA_VECTOR_ALLOCATION_FAIL) {
    throw std::bad_alloc();
  }
  else if (rc == FASTA_VECTOR_FILE_OPEN_FAIL) {
    throw "Could not open fasta file for reading.";
  }
  else if (rc == FASTA_VECTOR_FILE_READ_FAIL) {
    throw "Error while reading from the opened fasta file.";
  }

  shared_ptr<SequencePreprocessor> sequencePreprocessor =
    std::make_shared<SequencePreprocessor>(this->fastaVector.get());


  hwClient->writeSequence(sequencePreprocessor->getCompressedSequenceBuffer());
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

shared_ptr<vector<VerifiedHit>> Havac::getHitsFromFinishedRun() {
  shared_ptr<vector<HardwareHitReport>> unverifiedHitList = this->hwClient->getHitReportList();
  shared_ptr<HitVerifier> hitVerifier = std::make_shared<HitVerifier>(this->fastaVector,
    this->p7HmmList);

  return hitVerifier->verify(unverifiedHitList);
}


ert_cmd_state Havac::currentHardwareState() {
  return hwClient->getHwState();
}

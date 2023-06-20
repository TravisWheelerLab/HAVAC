#include "Havac.hpp"
#include "../PhmmReprojection/PhmmReprojection.h"
#include "phmm/PhmmPreprocessor.hpp"
#include "sequence/SequencePreprocessor.hpp"

#include <stdexcept>
#include <exception>
#include <array>
#include <chrono>

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

shared_ptr<vector<VerifiedHit>> Havac::getHitsFromFinishedRun() {
  shared_ptr<vector<HardwareHitReport>> unverifiedHitList = this->hwClient->getHitReportList();
  shared_ptr<HitVerifier> hitVerifier = std::make_shared<HitVerifier>(this->fastaVector,
    this->p7HmmList, this->compressedPhmmScores);

  auto verifiedHits = hitVerifier->verify(unverifiedHitList, this->requiredPValue);

  return verifiedHits;
}


enum havac_cmd_state Havac::currentHardwareState() {
  return (havac_cmd_state) hwClient->getHwState();
}

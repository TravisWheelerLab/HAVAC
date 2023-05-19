#include "HavacHwClient.hpp"
#include "HavacSSV.hpp"
#include "types/HardwareHitReport.hpp"
#include "../device/PublicDefines.h"
#include <string>
#include <cassert>
extern "C" {
  #include <FastaVector.h>
  #include <p7HmmReader.h>
}


using std::shared_ptr;
using std::vector;
using std::tuple;

//group IDs are used by the XRT framework to identify what input the data goes to
#define HAVAC_SEQUENCE_BUFFER_GROUP_ID 0
#define HAVAC_PHMM_BUFFER_GROUP_ID 2
#define HAVAC_HIT_REPORT_BUFFER_GROUP_ID 4
#define HAVAC_HIT_REPORT_COUNT_BUFFER_GROUP_ID 5


HavacHwClient::HavacHwClient(const std::string& xclbinFileSrc, const std::string& havacKernelName, const uint32_t deviceIndex) {
  //identify the device and load the xclbin file.
  this->havacDevice = xrt::device(deviceIndex);

  this->generateKernel(xclbinFileSrc, havacKernelName);
  this->allocateBuffers();
}

void HavacHwClient::generateKernel(const std::string& xclbinFileSrc, const std::string& havacKernelName) {
  if (!havacDevice) {
    throw std::logic_error("The havac device was not initialized prior to generating the kernel.");
  }

  try {
    this->havacUuid = havacDevice->load_xclbin(xclbinFileSrc);
  }
  catch (const std::exception& e) {
    std::cerr << "ERROR: HavacHwClient failed to load xclbin file " << xclbinFileSrc << ". does this file exist?\n" << e.what() << std::endl;
    throw;
  }

  //now load the kernel with the xclbin we just setup
  try {
    this->havacKernel = xrt::kernel(*havacDevice, *havacUuid, havacKernelName);
  }
  catch (const std::exception& e) {
    std::cerr << "ERROR: HavacHwClient was unable to build the kernel. perhaps the kernel name is wrong?\n" << e.what() << std::endl;
  }
}

void HavacHwClient::allocateBuffer(const int argumentIndex, boost::optional<xrt::bo>& buffer, const uint64_t sizeInBytes) {
  int bankGroup;
  try {
    bankGroup = this->havacKernel->group_id(argumentIndex);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: could not generate bank group for argument index " << argumentIndex << "." << std::endl;
    throw;
  }
  try {
    buffer = xrt::bo(this->havacDevice.get(), sizeInBytes, bankGroup);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: failure to allocate memory for device buffer on argument index " << argumentIndex << std::endl;
    throw;
  }
}

void HavacHwClient::allocateBuffers() {
  this->allocateBuffer(HAVAC_SEQUENCE_BUFFER_GROUP_ID, this->sequenceBuffer, this->sequenceAllocationSizeInBytes);
  this->allocateBuffer(HAVAC_PHMM_BUFFER_GROUP_ID, this->phmmBuffer, this->phmmAllocationSizeInBytes);
  this->allocateBuffer(HAVAC_HIT_REPORT_BUFFER_GROUP_ID, this->hitReportBuffer, this->hitReportAllocationSizeInBytes);
  this->allocateBuffer(HAVAC_HIT_REPORT_COUNT_BUFFER_GROUP_ID, this->hitReportCountBuffer, sizeof(uint32_t));
}

void HavacHwClient::writeSequence(const vector<uint8_t>& compressedSequence) {
  //get the length of the sequence, in segments. this will throw if the 
  //compressed sequence isn't divisble into even segments. 
  uint32_t numSymbolsInSequence = compressedSequence.size() * 4;  //2-bit symbols = 4 symbols/byte
  uint32_t numSequenceSegments = numSymbolsInSequence / NUM_CELL_PROCESSORS;
  if (numSymbolsInSequence != (numSequenceSegments * NUM_CELL_PROCESSORS)) {
    std::ostringstream stringStream;
    stringStream << "sequence length must be a multiple of the sequence segment length " <<
      std::to_string(NUM_CELL_PROCESSORS) << " but sequence of length " << std::to_string(numSymbolsInSequence) <<
      " was not \n";
    throw std::length_error(stringStream.str());
  }
  this->sequenceLengthInSegments = numSequenceSegments;

  try {
    const size_t bufferOffset = 0;
    const size_t compressedSequenceSize = compressedSequence.size();
    sequenceBuffer->write(compressedSequence.data(), compressedSequence.size(), bufferOffset);
    sequenceBuffer->sync(XCL_BO_SYNC_BO_TO_DEVICE);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: failed to write sequence data to fpga client.\n" << e.what() << std::endl;
    throw;
  }
}
void HavacHwClient::writePhmm(std::shared_ptr<vector<int8_t>> phmmAsFlattenedArray) {
  uint32_t phmmLengthInVectors = phmmAsFlattenedArray->size() / 4;
  if (phmmAsFlattenedArray->size() != phmmLengthInVectors * 4) {
    std::ostringstream stringStream;
    stringStream << "phmm of length length " << phmmAsFlattenedArray->size() <<
      " was given, but was not divisible by 4. only "
      "nucleotide phmms with 4 scores/position is supported.\n";
    throw std::length_error(stringStream.str());
  }
  this->phmmLengthInVectors = phmmLengthInVectors;

  try {
    const size_t bufferOffset = 0;
    //dereference the phmm shared_ptr to get the actual array data. 
    phmmBuffer->write(phmmAsFlattenedArray->data(), phmmAsFlattenedArray->size(), bufferOffset);
    phmmBuffer->sync(XCL_BO_SYNC_BO_TO_DEVICE);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: failed to write phmm to fpga client.\n" << e.what() << std::endl;
    throw;
  }
}


void HavacHwClient::invokeHavacSsvAsync() {
  if (this->sequenceLengthInSegments == 0) {
    throw std::length_error("sequence length in segments cannot be 0, but 0 was given to the client.");
  }
  if (this->phmmLengthInVectors == 0) {
    throw std::length_error("phmm length in vectors cannot be 0, but 0 was given to the client.");
  }

  this->havacRunObject = this->havacKernel.get()(this->sequenceBuffer.get(), this->sequenceLengthInSegments, this->phmmBuffer.get(),
    this->phmmLengthInVectors, this->hitReportBuffer.get(), this->hitReportCountBuffer.get());
}

ert_cmd_state HavacHwClient::waitForHavacSsvAsync(
  const std::chrono::milliseconds& timeout) {

  return this->havacRunObject->wait(timeout);
}

ert_cmd_state HavacHwClient::abort() {
  return havacRunObject->abort();
}

ert_cmd_state HavacHwClient::getHwState() {
  if (havacRunObject) {
    return havacRunObject->state();
  }
  else {
    throw std::logic_error("run object was not initialized. run function invokeHavacSsvAsync to initialize this object.");
  }
}

uint32_t HavacHwClient::getNumHits() {
  uint32_t numHits = -1;
  try {
    this->hitReportCountBuffer->sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    this->hitReportCountBuffer->read(&numHits, sizeof(uint32_t), 0);
  }
  catch (std::exception&) {
    throw std::runtime_error("ERROR: unable to read hit report count from hardware client.");
  }
  if (numHits == (uint32_t)-1) {
    throw std::runtime_error("num hits was not set by the client!");
  }
  return numHits;
}

shared_ptr<vector<HardwareHitReport>> HavacHwClient::getHitReportList() {
  uint32_t numHits = this->getNumHits();
  shared_ptr<vector<HardwareHitReport>> hitReportList =
    std::make_shared<vector<HardwareHitReport>>(numHits);

  try {
    this->hitReportBuffer->sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    hitReportBuffer->read(hitReportList->data(), sizeof(HardwareHitReport) * numHits, 0);
  }
  catch (std::exception&) {
    std::cerr << "ERROR: unable to read hit reports from hardware client" << std::endl;
    throw;
  }

  return hitReportList;
}
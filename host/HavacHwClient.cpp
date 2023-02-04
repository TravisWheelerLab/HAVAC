#include "HavacHwClient.hpp"
#include "HavacSSV.hpp"
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

#define HAVAC_SEQUENCE_BUFFER_GROUP_ID 0
#define HAVAC_SEQUENCE_LENGTH_GROUP_ID 1
#define HAVAC_PHMM_BUFFER_GROUP_ID 2
#define HAVAC_PHMM_LENGTH_GROUP_ID 3
#define HAVAC_HIT_REPORT_BUFFER_GROUP_ID 4
#define HAVAC_HIT_REPORT_LENGTH_GROUP_ID 5

//constructor
HavacHwClient::HavacHwClient(const std::string& xclbinFileSrc, const std::string& havacKernelName, const uint32_t deviceIndex = 0) {
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

void HavacHwClient::allocateBuffers() {
  assert((uint64_t)this->sequenceAllocationSizeInBytes <= 4294967295ULL && "Assert fail: sequence buffer length must be less than 4GiB.");
  assert((uint64_t)this->phmmAllocationSizeInBytes <= 4294967295ULL && "Assert fail: phmm buffer length must be less than 4GiB.");

  int sequenceBankGroup, sequenceLenBankGroup, phmmBankGroup, phmmLenBankGroup, hitReportBankGroup;
  try {
    auto sequenceBankGroup = havacKernel->group_id(HAVAC_SEQUENCE_BUFFER_GROUP_ID);
    auto phmmBankGroup = havacKernel->group_id(HAVAC_PHMM_BUFFER_GROUP_ID);
    auto hitReportBankGroup = havacKernel->group_id(HAVAC_HIT_REPORT_BUFFER_GROUP_ID);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: could not generate bank group index for one of the buffers\n" << e.what() << std::endl;
    throw;
  }
  try {
    this->sequenceBuffer = xrt::bo(*havacDevice, this->sequenceAllocationSizeInBytes, sequenceBankGroup);
    this->phmmBuffer = xrt::bo(*havacDevice, this->phmmAllocationSizeInBytes, sequenceBankGroup);
    this->hitReportBuffer = xrt::bo(*havacDevice, this->hitReportAllocationSizeInBytes, hitReportBankGroup);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: could not allocate memory for one or more hardware-side buffers\n" << e.what() << std::endl;
    throw;
  }
}
void HavacHwClient::writeSequence(const uint8_t* sequenceAsEncodedBytes, const uint32_t sequenceLengthInBytes) {
  try {
    const size_t bufferOffset = 0;
    phmmBuffer->write(sequenceAsEncodedBytes, sequenceLengthInBytes, bufferOffset);
    phmmBuffer->sync(XCL_BO_SYNC_BO_TO_DEVICE);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: failed to write sequence data to fpga client.\n" << e.what() << std::endl;
    throw;
  }
}
void HavacHwClient::writePhmm(const int8_t* phmmAsFlattenedArray, const uint32_t phmmLengthInBytes) {
  try {
    const size_t bufferOffset = 0;
    phmmBuffer->write(phmmAsFlattenedArray.get(), phmmLengthInBytes, bufferOffset);
    phmmBuffer->sync(XCL_BO_SYNC_BO_TO_DEVICE);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: failed to write phmm to fpga client.\n" << e.what() << std::endl;
    throw;
  }
}


void HavacHwClient::invokeHavacSsvAsync() {
  //reset the num hits, just in case. this may also help prevent errors.
  this->numHits = 0;

  this->havacRunObject = (*this->havacKernel)(*this->sequenceBuffer, this->sequenceLengthInSegments, *this->phmmBuffer,
    this->phmmLengthInVectors, *this->hitReportBuffer, this->numHits);
}

void HavacHwClient::abort() {
  havacRunObject->abort();
}

ert_cmd_state HavacHwClient::getHwState() {
  if (havacRunObject) {
    return havacRunObject->state();
  }
  else {
    throw std::logic_error("run object was not initialized. run function invokeHavacSsvAsync to initialize this object.");
  }
}

shared_ptr<vector<HavacHardwareHitReport>> HavacHwClient::getHitReportList() {
  shared_ptr<HavacHardwareHitReport[]> hitReportVector(new vector<HavacHardwareHitReport>(numHits));


  try {
    hitReportBuffer->read(hitReportVector->data(), sizeof(HavacHardwareHitReport) * this->numHits, 0);
  }
  catch (std::exception&) {
    std::cerr << "ERROR: unable to read hit reports from hardware client" << std::endl;
    throw;
  }

  return hitReportVector;
}
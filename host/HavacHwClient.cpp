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

#define HAVAC_SEQUENCE_BUFFER_GROUP_ID 0
#define HAVAC_SEQUENCE_LENGTH_GROUP_ID 1
#define HAVAC_PHMM_BUFFER_GROUP_ID 2
#define HAVAC_PHMM_LENGTH_GROUP_ID 3
#define HAVAC_HIT_REPORT_BUFFER_GROUP_ID 4
#define HAVAC_HIT_REPORT_LENGTH_GROUP_ID 5

//constructor
HavacHwClient::HavacHwClient(const std::string& xclbinFileSrc, const std::string& havacKernelName, const uint32_t deviceIndex) {
  //identify the device and load the xclbin file.
  this->havacDevice = xrt::device(deviceIndex);
  this->generateKernel(xclbinFileSrc, havacKernelName);
  this->allocateBuffers();
  this->numHits = -1;
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

void HavacHwClient::writeSequence(std::shared_ptr<vector<uint8_t>> compressedSequence) {
  try {
    const size_t bufferOffset = 0;
    phmmBuffer->write(compressedSequence->data(), sizeof(compressedSequence), bufferOffset);
    phmmBuffer->sync(XCL_BO_SYNC_BO_TO_DEVICE);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: failed to write sequence data to fpga client.\n" << e.what() << std::endl;
    throw;
  }
}
void HavacHwClient::writePhmm(std::shared_ptr<vector<int8_t>> phmmAsFlattenedArray) {
  try {
    const size_t bufferOffset = 0;
    //dereference the phmm shared_ptr to get the actual array data. 
    phmmBuffer->write(phmmAsFlattenedArray->data(), sizeof(phmmAsFlattenedArray), bufferOffset);
    phmmBuffer->sync(XCL_BO_SYNC_BO_TO_DEVICE);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: failed to write phmm to fpga client.\n" << e.what() << std::endl;
    throw;
  }
}


void HavacHwClient::invokeHavacSsvAsync() {
  //reset the num hits, just in case. this allows us to check to make sure it was set by the hardware.
  //I considered using a std::optional<> type, but getting the runObject to work with an optional_T would be a pain.
  this->numHits = -1;

  this->havacRunObject = (*this->havacKernel)(*this->sequenceBuffer, this->sequenceLengthInSegments, *this->phmmBuffer,
    this->phmmLengthInVectors, *this->hitReportBuffer, this->numHits);
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

shared_ptr<vector<HardwareHitReport>> HavacHwClient::getHitReportList() {
  if (this->numHits < 0) {
    throw std::logic_error("number of hits was not initialized, did the hardare fail to run?");
  }
  shared_ptr<vector<HardwareHitReport>> hitReportList = std::make_shared<vector<HardwareHitReport>>(numHits);
  // shared_ptr<HavacHardwareHitReport[]> hitReportList(new HavacHardwareHitReport[numHits]);

  try {
    hitReportBuffer->read(hitReportList->data(), sizeof(HardwareHitReport) * this->numHits, 0);
  }
  catch (std::exception&) {
    std::cerr << "ERROR: unable to read hit reports from hardware client" << std::endl;
    throw;
  }

  return hitReportList;
}
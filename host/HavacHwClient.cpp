#include "HavacHwClient.hpp"
#include "HavacSSV.hpp"
#include "../device/PublicDefines.h"
#include <string>
#include <cassert>
extern "C" {
  #include <FastaVector.h>
  #include <p7HmmReader.h>
}

#define HAVAC_SEQUENCE_BUFFER_GROUP_ID 0
#define HAVAC_SEQUENCE_LENGTH_GROUP_ID 1
#define HAVAC_PHMM_BUFFER_GROUP_ID 2
#define HAVAC_PHMM_LENGTH_GROUP_ID 3
#define HAVAC_HIT_REPORT_BUFFER_GROUP_ID 4
#define HAVAC_HIT_REPORT_LENGTH_GROUP_ID 5

//constructor
HavacHwClient::HavacHwClient(const uint32_t deviceIndex) {
  //identify the device and load the xclbin file.
  havacDevice = xrt::device(deviceIndex);
}


void HavacHwClient::generateKernel(const std::string& xclbinFileSrc, const std::string& havacKernelName) {
  if (!havacDevice) {
    throw std::logic_error("The havac device was not initialized prior to generating the kernel.");
  }

  try {
    havacUuid = havacDevice->load_xclbin(xclbinFileSrc);
  }
  catch (const std::exception& e) {
    std::cerr << "ERROR: HavacHwClient failed to load xclbin file " << xclbinFileSrc << ". does this file exist?\n"<< e.what() << std::endl;
    throw;
  }

  //now load the kernel with the xclbin we just setup
  try {
    havacKernel = xrt::kernel(*havacDevice, *havacUuid, havacKernelName);
  }
  catch (const std::exception& e) {
    std::cerr << "ERROR: HavacHwClient was unable to build the kernel. perhaps the kernel name is wrong?\n" << e.what() << std::endl;
  }
}

void HavacHwClient::allocateBuffers(const uint32_t sequenceBufferLength, const uint32_t phmmBufferLength,
  const uint32_t hitReportBufferLength) {
  assert((uint64_t)sequenceBufferLength <= 4294967295ULL && "Assert fail: sequence buffer length must be less than 4GiB.");
  assert((uint64_t)phmmBufferLength <= 4294967295ULL && "Assert fail: phmm buffer length must be less than 4GiB.");

  int sequenceBankGroup, sequenceLenBankGroup, phmmBankGroup, phmmLenBankGroup, hitReportBankGroup, hitReportLenBankGroup;
  try {
    auto sequenceBankGroup = havacKernel->group_id(HAVAC_SEQUENCE_BUFFER_GROUP_ID);
    auto phmmBankGroup = havacKernel->group_id(HAVAC_PHMM_BUFFER_GROUP_ID);
    auto hitReportBankGroup = havacKernel->group_id(HAVAC_HIT_REPORT_BUFFER_GROUP_ID);
    auto hitReportLenBankGroup = havacKernel->group_id(HAVAC_HIT_REPORT_LENGTH_GROUP_ID);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: could not generate bank group index for one of the buffers\n" << e.what() << std::endl;
    throw;
  }
  try {
    sequenceBuffer = xrt::bo(*havacDevice, sequenceBufferLength, sequenceBankGroup);
    phmmBuffer = xrt::bo(*havacDevice, phmmBufferLength, sequenceBankGroup);
    hitReportBuffer = xrt::bo(*havacDevice, hitReportBufferLength, hitReportBankGroup);
    numHitReportsBuffer = xrt::bo(*havacDevice, sizeof(uint32_t), hitReportLenBankGroup);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: could not allocate memory for one or more hardware-side buffers\n" << e.what() << std::endl;
    throw;
  }
}
void HavacHwClient::writeSequence(const uint8_t* sequenceAsEncodedBytes, const uint32_t sequenceLengthInBytes){
  try {
    phmmBuffer->write(sequenceAsEncodedBytes, sequenceLengthInBytes, 0);
    phmmBuffer->sync(XCL_BO_SYNC_BO_TO_DEVICE);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR: failed to write sequence data to fpga client.\n" << e.what() << std::endl;
    throw;
  }
}
void HavacHwClient::writePhmm(const uint8_t* phmmAsFlattenedArray, const uint32_t phmmLengthInBytes){
  try {
    phmmBuffer->write(phmmAsFlattenedArray, phmmLengthInBytes, 0);
    phmmBuffer->sync(XCL_BO_SYNC_BO_TO_DEVICE);
  }
  catch(std::exception &e){
    std::cerr << "ERROR: failed to write phmm to fpga client.\n" << e.what()<< std::endl;
    throw;
  }
}


void HavacHwClient::invokeHavacSsvAsync() {
  
  this->havacRunObject = (*this->havacKernel)(*this->sequenceBuffer, this->sequenceLengthInSegments, *this->phmmBuffer,
    this->phmmLengthInVectors, *this->hitReportBuffer, *this->numHitReportsBuffer);
}

void HavacHwClient::abort(){
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

boost::optional<std::vector<struct HavacHostHitReport>>& HavacHwClient::getHitReportList() {
  boost::optional<std::vector<struct HavacHostHitReport>> hitReportVector;
  //get the numHitReports
  uint32_t numHitReports = 0;
  try {
    numHitReportsBuffer->read(&numHitReports, sizeof(uint32_t), 0);
  }
  catch (std::exception&) {
    std::cerr << "ERROR: unable to read the number of hit reports " << std::endl;
    throw;
  }

  hitReportVector = std::vector<struct HavacHostHitReport>(numHitReports);
  try {
    hitReportBuffer->read(hitReportVector->data(), sizeof(struct HavacHostHitReport) * numHitReports, 0);
  }
  catch (std::exception&) {
    std::cerr << "ERROR: unable to read hit reports from hardware client" << std::endl;
    throw;
  }

  return hitReportVector;
}
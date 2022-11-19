#include "HavacSSV.hpp"
#include "../device/PublicDefines.h"
#include <iostream>
#include <cstring>
#include <cassert>
#include <cstdint>
#include <vector>
#include <boost/optional.hpp>

// XRT includes
#include "xrt/xrt_bo.h"
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"
#include <experimental/xrt_xclbin.h>
#include <experimental/xrt_ip.h>

#define HAVAC_KERNEL_NAME "HavacKernel"
#define NUM_NUCS_PADDING_BETWEEN_SEQUENCES 8


std::vector<struct HavacHostHitReport>& invokeHavac(const std::vector<uint8_t>& fullPaddedSequence,
  const std::vector<uint8_t>& phmmMatchScores, const std::string& xclbinSrc);
xrt::run& invokeHavacAsync(const std::vector<uint8_t>& fullPaddedSequence,
  const std::vector<uint8_t>& phmmMatchScores, const std::string& xclbinSrc);

boost::optional<std::vector<struct HavacHostHitReport>> havacRunSsv(const struct FastaVector& fastaVector, const struct P7HmmList& phmmList,
  float pValueThreshold, const std::string& xclbinSrc) {


  uint32_t deviceIndex = 0; //oh god how do you get this device index if it's not zero?
  //none of the documentation (of what little exists) uses anything but 0, and even reading the source code
  //doesn't elucidate the deviceIndex, it's just a publically visible member data that doesn't get used in the class




  xrt::run havacRuntimeObject = havacRunSsvAsync(fastaVector, phmmList, pValueThreshold, xclbinSrc);
  ert_cmd_state finishedRunState = havacRuntimeObject.wait();

  if (finishedRunState == ERT_CMD_STATE_TIMEOUT) {
    std::cout << "error upon invoking Havac SSV, timeout occurred." << std::endl;
    return boost::none;
  }
  else if (finishedRunState == ERT_CMD_STATE_NORESPONSE) {
    std::cout << "error upon invoking Havac SSV, no response from HAVAC coprocessor." << std::endl;
    return boost::none;
  }
  else if (finishedRunState == ERT_CMD_STATE_ABORT) {
    std::cout << "error upon invoking Havac SSV, HAVAC coprocessor task aborted." << std::endl;
    return boost::none;
  }
  else if (finishedRunState == ERT_CMD_STATE_ERROR) {
    //fall-through error case to catch anything unexpected.
    std::cout << "error upon invoking Havac SSV, HAVAC coprocessor reported generic error (cmd state ERT_CMD_STATE_ERROR." << std::endl;
    return boost::none;
  }
  else if (finishedRunState != ERT_CMD_STATE_COMPLETED) {
    std::cout << "error upon invoking Havac SSV, state returned code " << finishedRunState << "." << std::endl;
    return boost::none;
  }
  else {

    auto numHitsbuffer = xrt::bo(device, buffer_size_in_bytes, bank_grp_idx_0);
    uint32_t numHitsReported;
    numHitsOutputBuffer.sync();
    numHitsOutputBuffer.read(&numHitsReported);

    allocateBufferForHitReports(...);
    check for good buffer
      bo::sync();
    use bo::read(dst, numHitsReported * sizeof(struct HitReport), 0);
  }

}


xrt::run& havacRunSsvAsync(const struct FastaVector& fastaVector, const struct P7HmmList& phmmList,
  float pValueThreshold, const std::string& xclbinSrc) {
  //puts random padding characters between sequences to help zero out the scores between sequences
    //also makes sure that the total sequence length is a multiple of the number of cell processors 
  std::vector<uint8_t> sequenceAsPaddedString = makePaddedSequenceFromFastaVector(fastaVector, NUM_NUCS_PADDING_BETWEEN_SEQUENCES);
  //aggregates the individual phmms in the list, adds a [-127,-127,-127,-127] vector between phmms as a seperator,
  //and reprojects all the phmms such that they pass the threshold value at exactly 256.
  std::vector<uint8_t> phmmVectorList = createProjectedPhmmFromPhmmList(phmmList, pValueThreshold);

  //setup the structs needed to invoke the hardware

  auto xclbinFile = xrt::xclbin(xclbinSrc);
  int device_index = 0;
  auto device = xrt::device(device_index);
  auto uuid = device.load_xclbin(xclbinFile);

  uint32_t sequenceLengthInBytes = sequenceAsPaddedString.size();
  uint32_t sequenceLengthInSegments = sequenceLengthInBytes / NUM_CELL_PROCESSORS;
  if (sequenceLengthInBytes % NUM_CELL_PROCESSORS != 0) {
    std::cout << "ERROR: sequence byte length " << sequenceLengthInBytes << "was not divisible by the number of cell processors " <<
      NUM_CELL_PROCESSORS << "(remainder " << (sequenceLengthInBytes % NUM_CELL_PROCESSORS) << ")" << std::endl;
  }

  uint32_t phmmLengthInBytes = phmmVectorList.size();

  auto kernel = xrt::kernel(device, uuid, HAVAC_KERNEL_NAME);
  auto bank_grp_arg0 = kernel.group_id(0); // Memory bank index for kernel argument 0
  kernel()
}

std::vector<struct HitReport>& invokeHavac(const std::vector<uint8_t>& fullPaddedSequence,
  const std::vector<uint8_t>& phmmMatchScores, const std::string& xclbinSrc) {

  auto xclbinFile = xrt::xclbin(xclbinSrc);
  int device_index = 0;
  auto device = xrt::device(device_index);
  auto uuid = device.load_xclbin(xclbinFile);

  uint32_t sequenceLengthInBytes = fullPaddedSequence.size();
  uint32_t phmmLengthInBytes = phmmMatchScores.size();

  auto kernel = xrt::kernel(device, uuid, HAVAC_KERNEL_NAME);
  auto bank_grp_arg0 = kernel.group_id(0); // Memory bank index for kernel argument 0
  kernel()
}


xrt::run& invokeHavacAsync(const std::vector<uint8_t>& fullPaddedSequence,
  const std::vector<uint8_t>& phmmMatchScores, const std::string& xclbinSrc) {

  auto xclbinFile = xrt::xclbin(xclbinSrc);
  int device_index = 0;
  auto device = xrt::device(device_index);
  auto uuid = device.load_xclbin(xclbinFile);

  uint32_t sequenceLengthInBytes = fullPaddedSequence.size();
  uint32_t phmmLengthInBytes = phmmMatchScores.size();

  auto kernel = xrt::kernel(device, uuid, HAVAC_KERNEL_NAME);
  auto bank_grp_arg0 = kernel.group_id(0); // Memory bank index for kernel argument 0
}

struct HavacHostBuffers& havacTryAllocateBuffers(xrt::device& device, xrt::kernel& kernel) {
  const uint64_t maxSequenceBufferLengthInBytes = 6 * 1024 * 1024 * 1024; //6GiB
  const uint64_t maxPhmmBufferLengthInBytes = 4 * HAVAC_MAX_SUPPORTED_PHMM_LENGTH;
  const uint64_t maxHitReportMemorySizeInBytes = 512 * 1024 * 1024;  //.5GiB
  const uint8_t sequenceKernelGroupId = 0;
  const uint8_t phmmKernelGroupId = 2;
  const uint8_t hitReportKernelGroupId = 4;

  struct HavacHostBuffers buffers;
  try {
    buffers.sequenceBuffer = xrt::bo(device, maxSequenceBufferLengthInBytes, kernel.group_id(sequenceKernelGroupId));
  }
  catch (...) {
    std::cout << "error: exception generated when attempting to allocate host buffer for sequence data" << std::endl;
    buffers.sequenceBuffer = boost::none;
  }

  try {
    buffers.phmmBuffer = xrt::bo(device, maxPhmmBufferLengthInBytes, kernel.group_id(phmmKernelGroupId));
  }
  catch (...) {
    std::cout << "error: exception generated when attempting to allocate host buffer for phmm data" << std::endl;
    buffers.phmmBuffer = boost::none;
  }

  try {
    buffers.hitReportBuffer = xrt::bo(device, maxHitReportMemorySizeInBytes, kernel.group_id(hitReportKernelGroupId));
  }
  catch (...) {
    std::cout << "error: exception generated when attempting to allocate host buffer for hit report data" << std::endl;
    buffers.hitReportBuffer = boost::none;
  }

  return buffers;
}

ert_cmd_state havacKernelRun(xrt::kernel& kernel, struct HavacHostBuffers& hostBuffers,
  uint32_t sequenceLengthInBytes, uint32_t phmmLengthInBytes) {
  auto run = kernel(hostBuffers.sequenceBuffer, sequenceLengthInBytes, hostBuffers.phmmBuffer,
    phmmLengthInBytes, hostBuffers.hitReportBuffer);
  run.wait();
  return run.state();
}

xrt::run& havacKernelRunAsync(xrt::kernel& kernel, struct HavacHostBuffers& hostBuffers,
  uint32_t sequenceLengthInBytes, uint32_t phmmLengthInBytes) {
  auto run = kernel(hostBuffers.sequenceBuffer, sequenceLengthInBytes, hostBuffers.phmmBuffer,
    phmmLengthInBytes, hostBuffers.hitReportBuffer);
  return run;
}

void havacKernelWaitForCompletion(xrt::run& kernelRunObject) {
  kernelRunObject.wait();
}


/**
 * ERT command state
 *
 * @ERT_CMD_STATE_NEW:         Set by host before submitting a command to
 *                             scheduler
 * @ERT_CMD_STATE_QUEUED:      Internal scheduler state
 * @ERT_CMD_STATE_SUBMITTED:   Internal scheduler state
 * @ERT_CMD_STATE_RUNNING:     Internal scheduler state
 * @ERT_CMD_STATE_COMPLETED:   Set by scheduler when command completes
 * @ERT_CMD_STATE_ERROR:       Set by scheduler if command failed
 * @ERT_CMD_STATE_ABORT:       Set by scheduler if command abort
 * @ERT_CMD_STATE_TIMEOUT:     Set by scheduler if command timeout and reset
 * @ERT_CMD_STATE_NORESPONSE:  Set by scheduler if command timeout and fail to
 *                             reset
 */
ert_cmd_state havacKernelAbortRun(xrt::run& kernelRunObject) {
  return kernelRunObject.abort();
}

ert_cmd_state havacKernelState(xrt::run& kernelRunObject) {
  return kernelRunObject.state();
}
bool havacKernelIsRunning(xrt::run& kernelRunObject) {
  return havacKernelState(kernelRunObject) == ERT_CMD_STATE_RUNNING;
}

void havacWriteDpDataToBuffers(uint8_t* bitCompressedSequenceData, size_t sequenceLengthInBytes,
  uint8_t* projectedPhmmMatchScores, size_t phmmLengthInBytes, struct HavacHostBuffers& hostBuffers) {
  if (hostBuffers.sequenceBuffer) {
    havacWriteToBuffer(*(hostBuffers.sequenceBuffer), bitCompressedSequenceData, sequenceLengthInBytes);
  }
  else {
    std::cout << "error: could not write to sequence buffer; buffer had optional none value" << std::endl;
  }
  if (hostBuffers.phmmBuffer) {
    havacWriteToBuffer(*(hostBuffers.phmmBuffer), projectedPhmmMatchScores, phmmLengthInBytes);
  }
  else {
    std::cout << "error: could not write to phmm buffer; buffer had optional none value" << std::endl;
  }
}

void havacWriteToBuffer(xrt::bo& buffer, void* data, uint64_t dataLengthInBytes) {
  size_t bufferOffset = 0;
  buffer.write(data, dataLengthInBytes, 0);
  buffer.sync(XCL_BO_SYNC_BO_TO_DEVICE);
}

void havacReadHitReportsToArray(xrt::bo& hitReportBuffer, struct HavacHostHitReport* hitReportArray,
  uint32_t numHitReports) {
  size_t bufferOffset = 0;
  hitReportBuffer.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
  hitReportBuffer.read(hitReportArray, numHitReports * sizeof(struct HavacHostHitReport), bufferOffset);
}
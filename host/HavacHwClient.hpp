#ifndef HAVAC_HW_CLIENT_HPP
#define HAVAC_HW_CLIENT_HPP

#include "types/HardwareHitReport.hpp"
#include "HitVerifier.hpp"
#include <cstdint>
#include <string>
#include <memory>
#include <tuple>
#include <vector>

// XRT includes
#include "xrt/xrt_bo.h"
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"
#include <experimental/xrt_xclbin.h>
#include <experimental/xrt_ip.h>
#include <boost/optional.hpp>

using std::shared_ptr;
using std::vector;

class HavacHwClient {
public:
  HavacHwClient(const std::string& xclbinFileSrc, const std::string& havacKernelName, const uint32_t deviceIndex = 0);
  HavacHwClient(HavacHwClient&& hc) = delete;
  HavacHwClient(HavacHwClient& hc) = delete;

  void writeSequence(shared_ptr<vector<uint8_t>> compressedSequence);
  void writePhmm(shared_ptr<vector<int8_t>> phmmAsFlattenedArray);
  void invokeHavacSsvAsync();

  ert_cmd_state getHwState();
  ert_cmd_state waitForHavacSsvAsync(const std::chrono::milliseconds& timeout = std::chrono::milliseconds{ 0 });
  ert_cmd_state abort();

  std::shared_ptr<vector<HardwareHitReport>> getHitReportList();

protected:
  boost::optional<xrt::device> havacDevice;
  boost::optional<xrt::uuid> havacUuid;
  boost::optional<xrt::kernel> havacKernel;
  boost::optional<xrt::bo> sequenceBuffer;
  boost::optional<xrt::bo> phmmBuffer;
  boost::optional<xrt::bo> hitReportBuffer;
  boost::optional<xrt::bo> numHitReportsBuffer;
  boost::optional<xrt::run> havacRunObject;

  uint32_t sequenceLengthInSegments;
  uint32_t phmmLengthInVectors;
  uint32_t numHits;

  static const uint64_t sequenceAllocationSizeInBytes = 512UL * 1024UL * 1024UL * 8UL;  //4.0GiB
  static const uint32_t phmmAllocationSizeInBytes = 512 * 1024 * 1024;         //0.5GiB
  static const uint32_t hitReportAllocationSizeInBytes = 512 * 1024 * 1024;     //0.5GiB

  void generateKernel(const std::string& xclbinFileSrc, const std::string& havacKernelName);
  void allocateBuffers();
};

#endif
#ifndef HAVAC_HW_CLIENT_HPP
#define HAVAC_HW_CLIENT_HPP

// XRT includes
#include <cstdint>
#include <string>

#include "xrt/xrt_bo.h"
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"
#include <experimental/xrt_xclbin.h>
#include <experimental/xrt_ip.h>
#include <boost/optional.hpp>

class HavacHwClient {
public:
  HavacHwClient(const uint32_t deviceIndex);
  HavacHwClient(HavacHwClient&& hc) = delete;
  HavacHwClient(HavacHwClient& hc) = delete;

  void generateKernel(const std::string& xclbinFileSrc, const std::string& havacKernelName);

  void allocateBuffers(const uint32_t sequenceLengthInBytes, const uint32_t phmmLengthInBytes, const uint32_t hitReportBufferLengthInBytes);

  void writeSequence(const uint8_t *sequenceAsEncodedBytes, const uint32_t sequenceLengthInBytes);
  void writePhmm(const uint8_t *phmmAsFlattenedArray, const uint32_t phmmLengthInBytes);
  void invokeHavacSsvAsync();

  ert_cmd_state getHwState();

  void waitForHavacSsvAsync();
  void abort();

  boost::optional<std::vector<struct HavacHostHitReport>>& getHitReportList();


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

};

#endif
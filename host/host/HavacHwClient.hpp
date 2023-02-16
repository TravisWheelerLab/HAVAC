#ifndef HAVAC_HW_CLIENT_HPP
#define HAVAC_HW_CLIENT_HPP

// XRT includes
#include <cstdint>
#include <string>
#include <memory>
#include <tuple>

#include "xrt/xrt_bo.h"
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"
#include <experimental/xrt_xclbin.h>
#include <experimental/xrt_ip.h>
#include <boost/optional.hpp>

#include "HardwareHitList.hpp"


class HavacHwClient {
public:
  HavacHwClient(const std::string& xclbinFileSrc, const std::string& havacKernelName, const uint32_t deviceIndex = 0);
  HavacHwClient(HavacHwClient&& hc) = delete;
  HavacHwClient(HavacHwClient& hc) = delete;

  void writeSequence(const uint8_t* sequenceAsEncodedBytes, const uint32_t sequenceLengthInBytes);
  void writePhmm(int8_t* phmmAsFlattenedArray, const uint32_t phmmLengthInBytes);
  void invokeHavacSsvAsync();

  ert_cmd_state getHwState();

  void waitForHavacSsvAsync();
  void abort();

  std::shared_ptr<std::vector<HavacHardwareHitReport>> getHitReportList();


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
  std::optional<uint32_t> numHits;


  static const uint32_t sequenceAllocationSizeInBytes = 512 * 1024 * 1024 * 8;  //4.0GiB
  static const  uint32_t phmmAllocationSizeInBytes = 512 * 1024 * 1024;         //0.5GiB
  static const  uint32_t hitReportAllocationSizeInBytes = 512 * 1024 * 1024;     //0.5GiB


  void generateKernel(const std::string& xclbinFileSrc, const std::string& havacKernelName);
  void allocateBuffers();


};

#endif
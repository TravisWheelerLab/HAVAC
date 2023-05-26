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
  /// @brief Constructor for HavacHwClient class
  /// @param xclbinFileSrc file source where the hardware design's .xclbin file is found
  /// @param havacKernelName "name" of the kernel. XRT does not explain what this should be 
  /// @param deviceIndex index of the device. this is probably 0, but once again, 
  ///         XRT does not explain what else it could be.
  HavacHwClient(const std::string& xclbinFileSrc, const std::string& havacKernelName, const uint32_t deviceIndex = 0);
  HavacHwClient(HavacHwClient&& hc) = delete;
  HavacHwClient(HavacHwClient& hc) = delete;

  /// @brief writes the given sequence to the DDR memory on the FPGA board
  /// @param compressedSequence sequence to write, as expressed as bit-compressed nucleotides
  ///         each uint8_t should contain 4 contiguous nucleotides as 2-bit encodings. 
  void writeSequence(const vector<uint8_t>& compressedSequence);

  /// @brief writes the given phmm to the DDR memory on the FPGA board 
  /// @param phmmAsFlattenedArray concatenated phmm list, in flattened signed 8-bit scores.
  ///         These scores should have been reprojected using the PhmmReprojection code
  ///         such that a threshold hit is a score of exactly 256
  void writePhmm(shared_ptr<vector<int8_t>> phmmAsFlattenedArray);

  /// @brief requests that the HAVAC hardware begin a computation run on the data that was previously given.
  void invokeHavacSsvAsync();

  /// @brief gets the current state of the FPGA hardware.
  /// @return state of the hardware. possible returns are:
  ///   ERT_CMD_STATE_NEW = 1, 
  ///   ERT_CMD_STATE_QUEUED = 2,
  ///   ERT_CMD_STATE_RUNNING = 3,
  ///   ERT_CMD_STATE_COMPLETED = 4,
  ///   ERT_CMD_STATE_ERROR = 5,
  ///   ERT_CMD_STATE_ABORT = 6,
  ///   ERT_CMD_STATE_SUBMITTED = 7,
  ///   ERT_CMD_STATEIMEOUT = 8,
  ///   ERT_CMD_STATE_NORESPONSE = 9
  ert_cmd_state getHwState();

  /// @brief wait for the current HAVAC hardware run to complete, or for timeout. 
  ///         The XRT docs don't specify what happens if you attempt this while it isn't running,
  //          so just treat calling this while not running as undefined behavior
  /// @param timeout Time, in milliseconds to wait until it gives up and errors out.
  ///         returns ERT_CMD_STATE_TIMEOUT if timeout was triggered
  /// @return state of the hardware upon returning, see getHwState()
  ert_cmd_state waitForHavacSsvAsync(const std::chrono::milliseconds& timeout = std::chrono::milliseconds{ 0 });

  /// @brief aborts the current hardware run. The docs don't specify what happens if you call this 
  ///         while the FPGA isn't running, so it might be undefined behavior.
  /// @return state of the aborted command, which I assume should be ERT_CMD_STATE_ABORT.
  ert_cmd_state abort();

  /// @brief Once a HAVAC run has finished, this retrieves the list of hit reports generated by the hardware
  ///         Note, these are hardware hits, and need to be verified via the HitVerifier.
  /// @return shared_ptr to the list of hardware hits reported by HAVAC.
  std::shared_ptr<vector<HardwareHitReport>> getHitReportList();

protected:
  boost::optional<xrt::device> havacDevice;
  boost::optional<xrt::uuid> havacUuid;
  boost::optional<xrt::kernel> havacKernel;
  boost::optional<xrt::bo> sequenceBuffer;
  boost::optional<xrt::bo> phmmBuffer;
  boost::optional<xrt::bo> hitReportBuffer;
  boost::optional<xrt::bo> hitReportCountBuffer;
  boost::optional<xrt::run> havacRunObject;

  uint32_t sequenceLengthInSegments;
  uint32_t phmmLengthInVectors;

  // static const uint64_t sequenceAllocationSizeInBytes = 512UL * 1024UL * 1024UL * 8UL;  //4.0GiB
  // static const uint64_t sequenceAllocationSizeInBytes = 3UL * 1024UL * 1024UL * 1024UL;  //3GiB
  // static const uint32_t phmmAllocationSizeInBytes = 1UL * 1024UL * 1024UL;         //1MiB
  //HAVAC supports 1024*1024 positions in the phmm, so the allocated size should be 4x that (1byte/symbol/position)
  static const uint32_t hitReportAllocationSizeInBytes = 4UL * 1024UL * 1024UL;     //4MiB

  /// @brief attempts to generate the kernel for the given .xclbin file and kernel name
  /// @param xclbinFileSrc source where the .xclbin file is located
  /// @param havacKernelName "name" of the kernel, XRT docs don't specify further than this.
  void generateKernel(const std::string& xclbinFileSrc, const std::string& havacKernelName);

  /// @brief allocates buffers on the FPGA card for the phmm/sequence/hitReports.
  ///         the allocated buffers should be large enough for any possible inputs/outputs.
  void allocateOutputBuffers();

  /// @breif allocates an individual buffer on the FPGA.
  /// @param argumentIndex index of the argument to define a buffer for on the kernel's top-level function
  /// @param buffer buffer object that will represent the allocated buffer on the hardware.
  /// @param sizeInBytes size to allocate on the hardware. Must be less than 4GiB. 
  void allocateBuffer(const int argumentIndex, boost::optional<xrt::bo>& buffer, const uint64_t sizeInBytes);

  /// @brief reads the number of hits generated from the client. 
  /// Calling this is the first step to reading the unverified hits generated by the client.
  /// @return number of hits the client generated. This will be the number of hit reports available to be read
  /// from the hitReportBuffer.
  uint32_t getNumHits();
};

#endif
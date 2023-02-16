#ifndef HAVAC_HOST_HPP
#define HAVAC_HOST_HPP

#include "HavacHwClient.hpp"
#include <string>
#include <memory>

extern "C"{
  #include <FastaVector.h>
  #include <p7HmmReader.h>
}

using std::shared_ptr;

struct HavacHitLocation {
  uint32_t sequenceIndex;
  uint32_t hitLocationInSequence;
  uint32_t phmmIndex;
  uint32_t hitLocationInPhmm;
};

class Havac {
public:
  /// @brief creates a HAVAC control object, which gives a simple interface to the hardware coprocessor
  /// @param deviceIndex index of the FPGA device, this should be 0 unless multiple cards are installed.
  ///   If this doesn't work, try iterating different deviceIndexes. Unfortunately, the XRT docs doesn't
  ///   speficy anything further, sorry.
  /// @param pValue 
  Havac(const uint32_t deviceIndex = 0, const float pValue = 0.05f);
  Havac(Havac&& havac) = delete;
  Havac(Havac& havac) = delete;
  ~Havac();

  /// @brief load the sequence data from the given fasta file source to the hardware client
  /// @param fastaSrc file location to load
  void loadSequence(const std::string fastaSrc);

  /// @brief load the phmm from the given file src to the hardware client
  /// @param phmmSrc location of file to read
  /// @param desiredPValue p value that would cause a hit to be registered
  void loadPhmm(const std::string phmmSrc, const float desiredPValue));

  /// @brief run the HAVAC hardware client synchronously, and return when finished
  void runHardwareClient();

  /// @brief run the HAVAC hardware client, and return immediately. 
  void runHardwareClientAsync();

  /// @brief while HAVAC client is running, wait until it has finished
  void waitHardwareClientAsync();

  /// @brief abort the current run of HAVAC
  void abortHardwareClient();

  /// @brief retrieve the hits from a finished HAVAC run
  /// @return the list of hits, after having been verified via bounded reference SSV checks.
  std::vector<HavacHitLocation>& getHitsFromFinishedRun();

  /// @brief gets the current hardware state
  /// @return state of the hardware, possible returns are:
  ///   ERT_CMD_STATE_NEW = 1, 
  ///   ERT_CMD_STATE_QUEUED = 2,
  ///   ERT_CMD_STATE_RUNNING = 3,
  ///   ERT_CMD_STATE_COMPLETED = 4,
  ///   ERT_CMD_STATE_ERROR = 5,
  ///   ERT_CMD_STATE_ABORT = 6,
  ///   ERT_CMD_STATE_SUBMITTED = 7,
  ///   ERT_CMD_STATEIMEOUT = 8,
  ///   ERT_CMD_STATE_NORESPONSE = 9
  ert_cmd_state currentHardwareState();

private:
  shared_ptr<HavacHwClient> hwClient;
  FastaVector *fastaVector;
  P7HmmList *p7HmmList;
  shared_ptr<int8_t[]> projectedPhmmMatchScores;
  uint32_t deviceIndex;
  bool phmmLoadedToDevice = false;
  bool sequenceLoadedToDevice = false;

  static const std::string havacXclbinFileSrc = "impl/havac.xclbin";
  static const std::string havacKernelName = "havacKernel";

  bool isRunning();
  void generateProjectedPhmmMatchScores(const float pValue = 0.05f);
};
#endif
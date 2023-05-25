#ifndef HAVAC_HOST_HPP
#define HAVAC_HOST_HPP

#include "HavacHwClient.hpp"
#include "HitVerifier.hpp"
extern "C"{
  #include <FastaVector.h>
  #include <p7HmmReader.h>
}
#include <string>
#include <memory>

using std::shared_ptr;
using std::vector;

enum havac_cmd_state{
  HAVAC_CMD_STATE_NEW = 1, 
  HAVAC_CMD_STATE_QUEUED = 2,
  HAVAC_CMD_STATE_RUNNING = 3,
  HAVAC_CMD_STATE_COMPLETED = 4,
  HAVAC_CMD_STATE_ERROR = 5,
  HAVAC_CMD_STATE_ABORT = 6,
  HAVAC_CMD_STATE_SUBMITTED = 7,
  HAVAC_CMD_STATEIMEOUT = 8,
  HAVAC_CMD_STATE_NORESPONSE = 9
};


const std::string xclbinSrcDefault = "impl/havac.xclbin";

class Havac {
public:
  /// @brief creates a HAVAC control object, which gives a simple interface to the hardware coprocessor
  /// @param deviceIndex index of the FPGA device, this should be 0 unless multiple cards are installed.
  ///   If this doesn't work, try iterating different deviceIndexes. Unfortunately, the XRT docs doesn't
  ///   speficy anything further, sorry.
  /// @param requiredPValue p value that would cause a hit to be registered. If the p value needs to be changed for some reason,
  /// just make another Havac object
  /// @param xclbinSrc 
  Havac(const uint32_t deviceIndex = 0, const float requiredPValue = 0.05f, const std::string xclbinSrc = xclbinSrcDefault);
  Havac(Havac&& havac) = delete;
  Havac(Havac& havac) = delete;
  ~Havac();

  /// @brief load the sequence data from the given fasta file source to the hardware client
  /// @param fastaSrc file location to load
  void loadSequence(const std::string fastaSrc);

  /// @brief load the phmm from the given file src to the hardware client
  /// @param phmmSrc location of file to read
  void loadPhmm(const std::string phmmSrc);

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
  shared_ptr<vector<VerifiedHit>> getHitsFromFinishedRun();

  /// @brief gets the current hardware state
  /// @return state of the hardware, possible returns are:
  ///   HAVAC_CMD_STATE_NEW = 1, 
  ///   HAVAC_CMD_STATE_QUEUED = 2,
  ///   HAVAC_CMD_STATE_RUNNING = 3,
  ///   HAVAC_CMD_STATE_COMPLETED = 4,
  ///   HAVAC_CMD_STATE_ERROR = 5,
  ///   HAVAC_CMD_STATE_ABORT = 6,
  ///   HAVAC_CMD_STATE_SUBMITTED = 7,
  ///   HAVAC_CMD_STATEIMEOUT = 8,
  ///   HAVAC_CMD_STATE_NORESPONSE = 9
  enum havac_cmd_state currentHardwareState();


private:
  shared_ptr<HavacHwClient> hwClient;
  FastaVector *fastaVector;
  P7HmmList *p7HmmList;
  shared_ptr<int8_t[]> projectedPhmmMatchScores;
  uint32_t deviceIndex;
  float requiredPValue;
  bool phmmLoadedToDevice = false;
  bool sequenceLoadedToDevice = false;
  uint32_t numCellGroupsConfigValue;

  std::string havacXclbinFileSrc;
  const std::string havacKernelName = "HavacKernel";
};
#endif
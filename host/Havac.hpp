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
using std::unique_ptr;

struct HavacHitLocation {
  uint32_t sequenceIndex;
  uint32_t hitLocationInSequence;
  uint32_t phmmIndex;
  uint32_t hitLocationInPhmm;
};

class Havac {
public:
  Havac(const uint32_t deviceIndex = 0, const float pValue = 0.05f);
  Havac(Havac&& havac) = delete;
  Havac(Havac& havac) = delete;
  ~Havac();

  void loadSequence(const std::string fastaSrc);
  void loadPhmm(const std::string phmmSrc, const float desiredPValue));
  void runHardwareClient();
  void runHardwareClientAsync();
  void waitHardwareClientAsync();
  void abortHardwareClient();
  std::vector<HavacHitLocation>& getHitsFromFinishedRun();

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
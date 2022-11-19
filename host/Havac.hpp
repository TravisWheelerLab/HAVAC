#ifndef HAVAC_HOST_HPP
#define HAVAC_HOST_HPP

#include "HavacHwClient.hpp"
#include <string>

struct HavacHitLocation{
  uint32_t sequenceIndex;
  uint32_t hitLocationInSequence;
  uint32_t phmmIndex; 
  uint32_t hitLocationInPhmm;
};

class Havac{
  public:
  Havac();
  Havac(Havac && havac) = delete;
  Havac(Havac &havac) = delete;

  void initHardwareClient(const uint32_t deviceIndex);
  void loadPhmmToHardware(const std::string &phmmFileSrc);
  void loadSequenceToHardware(const std::string& fastaFileSrc);
  void runHardwareClient();
  void runHardwareClientAsync();
  void waitHardwareClientAsync();
  std::vector<HavacHitLocation> &getHitsFromFinishedRun();

  private:
  HavacHwClient hwClient;
};
#endif
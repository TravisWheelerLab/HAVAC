#include "../../Havac.hpp"
#include "../../types/HavacHit.hpp"
#include <stdexcept>
#include <memory>
#include <vector>
#include<chrono>


using std::shared_ptr;
using std::vector;
using std::make_shared;

const float desiredPValue = 0.02f;
const int64_t MAX_ALLOWABLE_HIT_POSITIONAL_DIFFERENCE = 16;

void makeSequencesAndPhmms(uint32_t sequenceLength, uint32_t numSequences, const float subProbability);
void fileAppend(const char* fileReadFrom, const char* fileWriteTo);

const float SUB_PROBABILITY = 0.2f;
//args:
//0: program name
//1: xclbin src
//2: fasta src
//3: hmm src
int main(int argc, char** argv) {
  //initialize the RNG seed
  srand(time(NULL));

  //check to see if abort was requested.
  if (argc < 4) {
    std::cout << "error: program requires an argument for the source for the xclbin src, fasta file src, "
      "and hmm src" << std::endl;
    exit(2);
  }

  std::string xclbinSrc = argv[1];
  std::string fastaSrc = argv[2];
  std::string phmmSrc = argv[3];

  std::cout << "fasta src: " << fastaSrc << std::endl;
  std::cout << "phmmSrc: " << phmmSrc << std::endl;
  std::cout << "build start" << std::endl;
  auto havacBuildStart = std::chrono::high_resolution_clock::now();
  shared_ptr<Havac> havac = make_shared<Havac>(0, desiredPValue, xclbinSrc);
  auto havacBuildEnd = std::chrono::high_resolution_clock::now();

  std::cout << "data load start" << std::endl;
  auto dataLoadStart = std::chrono::high_resolution_clock::now();
  std::cout << "\t load phmm" << std::endl;
  havac->loadPhmm(phmmSrc);
  std::cout << "\t load sequence" << std::endl;
  havac->loadSequence(fastaSrc);
  auto dataLoadEnd = std::chrono::high_resolution_clock::now();

  std::cout << "beginning hw run" << std::endl;
  auto havacRunStart = std::chrono::high_resolution_clock::now();
  havac->runHardwareClient();
  auto havacRunEnd = std::chrono::high_resolution_clock::now();

  auto hitVerificationStart = std::chrono::high_resolution_clock::now();
  vector<HavacHit> hardwareHits = havac->getHitsFromFinishedRun();
  auto hitVerificationEnd = std::chrono::high_resolution_clock::now();

  std::cout << "hw generated "<< hardwareHits.size()<< " verified hits."<<std::endl;


  float buildTime = float(std::chrono::duration_cast <std::chrono::microseconds> (havacBuildEnd - havacBuildStart).count());
  float loadTime = float(std::chrono::duration_cast <std::chrono::microseconds> (dataLoadEnd - dataLoadStart).count());
  float runTime = float(std::chrono::duration_cast <std::chrono::microseconds> (havacRunEnd - havacRunStart).count());
  float verifyTime = float(std::chrono::duration_cast <std::chrono::microseconds> (hitVerificationEnd - hitVerificationStart).count());
  float fullTime = float(std::chrono::duration_cast <std::chrono::microseconds> (hitVerificationEnd - havacBuildStart).count());

  std::cout << "timing information:" << std::endl;
  std::cout << "havac build time " << buildTime << " microseconds (" << buildTime / 1000000.0f<< " seconds)." << std::endl;
  std::cout << "havac load time " << loadTime << " microseconds (" << loadTime / 1000000.0f << " seconds)." << std::endl;
  std::cout << "havac run time " << runTime << " microseconds (" << runTime / 1000000.0f << " seconds)." << std::endl;
  std::cout << "havac verify time " << verifyTime << " microseconds (" << verifyTime / 1000000.0f << " seconds)." << std::endl;
  std::cout << "total time taken " << fullTime << " microseconds (" << fullTime / 1000000.0f << " seconds)." << std::endl;


}

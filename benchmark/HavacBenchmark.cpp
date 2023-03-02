#include "HavacBenchmark.hpp"

#include <chrono>
#include <memory>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <getopt.h>
#include "../host/Havac.hpp"


using std::shared_ptr;
using std::vector;
using namespace std::chrono;

float pValueThreshold = 0.05f;
std::string phmmSrc;
std::string fastaSrc;
bool listHits = false;
bool showTimings = false;
int option;


int main(int argc, char** argv) {
  parseOptions(argc, argv);
  constexpr uint32_t deviceId = 0;
  shared_ptr<Havac> havac = std::make_shared<Havac>(deviceId);

  //load the phmm and sequence data to the fpga
  auto dataLoadStartTime = high_resolution_clock::now();
  havac->loadSequence(fastaSrc);
  havac->loadPhmm(phmmSrc, pValueThreshold);
  auto dataLoadEndTime = high_resolution_clock::now();
  auto dataLoadDurationTime = duration_cast<microseconds>(dataLoadEndTime - dataLoadStartTime);

  //execute the SSV run
  auto ssvRunStartTime = high_resolution_clock::now();
  havac->runHardwareClient();
  auto ssvRunEndTime = high_resolution_clock::now();
  auto ssvRunDurationTime = duration_cast<microseconds>(ssvRunEndTime - ssvRunStartTime);


  //get the hits from the hardware and verify them
  auto hitVerificationStartTime = high_resolution_clock::now();
  shared_ptr<vector<VerifiedHit>> verifiedHits = havac->getHitsFromFinishedRun();
  auto hitVerificationEndTimeTime = high_resolution_clock::now();
  auto hitVerificationDurationTime = duration_cast<microseconds>(hitVerificationEndTimeTime - hitVerificationStartTime);

  auto totalDurationTime = duration_cast<microseconds>(hitVerificationEndTimeTime - dataLoadStartTime);

  uint32_t numHitsReturned = verifiedHits->size();

  std::cout << numHitsReturned <<" hits returned." << std::endl;

  if(listHits){
    for (auto &hit: *verifiedHits){
      std::cout << "sequence #"<< hit.sequenceIndex<< ", pos: " << hit.sequencePosition << "; phmm #"<< 
        hit.phmmIndex << ", pos: "<< hit.phmmPosition << std::endl;
    }
  }

  if (showTimings) {
    std::cout << "data load time: " << dataLoadDurationTime.count() << "\nssv run time: " <<
      ssvRunDurationTime.count() << "\nhit verification time: " << hitVerificationDurationTime.count() <<
      "\n total time: " << totalDurationTime.count() << "\n all times in milliseconds" << std::endl;
  }
}


void parseOptions(int argc, char** argv) {
  // Define the options.
  struct option long_options[] = {
      {"f", required_argument, nullptr, 'f'},
      {"m", required_argument, nullptr, 'm'},
      {"p", optional_argument, nullptr, 'p'},
      {"p", no_argument, nullptr, 'l'},
      {"p", no_argument, nullptr, 't'},
      {nullptr, 0, nullptr, 0}
  };
  // Parse the options.
  while ((option = getopt_long(argc, argv, "p::f:m:lt", long_options, nullptr)) != -1) {
    switch (option) {
    case 'p':
      pValueThreshold = std::stof(optarg);
      break;
    case 'f':
      fastaSrc = optarg;
      break;
    case 'm':
      phmmSrc = optarg;
      break;
    case 'l':
      listHits = true;
      break;
    case 't':
      showTimings = true;
      break;
    default:
      std::cerr << "Invalid command line argument" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}
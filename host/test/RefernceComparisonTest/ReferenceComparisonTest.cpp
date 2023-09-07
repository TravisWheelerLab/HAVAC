#include "ReferenceComparisonTest.hpp"
#include "../../Havac.hpp"
#include "../Ssv.hpp"
#include "../../../test/generator/hmmSeqGenerator.h"
#include "../../../device/PublicDefines.h"
#include "../../../PhmmReprojection/PhmmReprojection.h"

extern "C" {
  #include <FastaVector.h>
  #include <p7HmmReader.h>
}


const char* phmmFileSrc = "refTest.hmm";
const char* seqFileSrc = "refTest.fasta";
const float desiredPValue = 0.02f;
// const uint32_t sequenceLength = NUM_CELL_PROCESSORS * 5;


void diagHitCheck(const FastaVector* fastaVector, const P7HmmList* phmmList, const float desiredPValue) {
  std::cout << "beginning diag hit check" << std::endl;
  float scoreMultiplier = findThreshold256ScalingFactor(&phmmList->phmms[0], desiredPValue);

  int16_t accumulatedScore = 0;
  for (size_t position = 0; position < phmmList->phmms[0].header.modelLength;position++) {
    char symbolFromFasta = fastaVector->sequence.charData[position];
    uint8_t symbolEncoding;
    switch (symbolFromFasta) {
    case 'a': case 'A': symbolEncoding = 0; break;
    case 'c': case 'C': symbolEncoding = 1; break;
    case 'g': case 'G': symbolEncoding = 2; break;
    default: symbolEncoding = 3;
    }

    float valFromPhmm = phmmList->phmms[0].model.matchEmissionScores[(4 * position) + symbolEncoding];
    float projectedScore = emissionScoreToProjectedScore(valFromPhmm, scoreMultiplier);
    int8_t projectedScoreAsInt8 = projectedScore;
    accumulatedScore += projectedScoreAsInt8;
    bool hasHit = false;
    if (accumulatedScore < 0) {
      accumulatedScore = 0;
    }
    else if (accumulatedScore >= 256) {
      hasHit = true;
    }
    if (hasHit) {
      accumulatedScore = 0;
    }
  }
}

int main(int argc, char** argv) {
  //initialize the RNG seed
  srand(time(NULL));
  
  bool abortRequested = false;
  std::cout << "beginning" << std::endl;
  uint32_t sequenceLength = 0;

  //run the hardware
  if (argc < 3) {
    std::cout << "error: program requires an argument for the source for the xclbin src, and length of sequence" << std::endl;
    exit(2);
  }

  abortRequested = std::string(argv[2]) == std::string("abort");
  if (abortRequested) {
    std::cout << "abort requested" << std::endl;
  }
  else {
    sequenceLength = atoi(argv[2]);
    if (sequenceLength == 0) {
      std::cout << "error: could not parse sequence length, or 0 given. either way, very illegal!" << std::endl;
      exit(3);
    }
  }
  std::string xclbinSrc = argv[1];
  std::cout << "generating havac object" << std::endl;
  shared_ptr<Havac> havac = make_shared<Havac>(0, desiredPValue, xclbinSrc);

  if (abortRequested) {
    std::cout << "abort requested, current hw state: " << havac->currentHardwareState() << "." << std::endl;
    havac->abortHardwareClient();
    std::cout << "aborted. hw state is: " << havac->currentHardwareState() << "." << std::endl;
    exit(5);
  }

  std::cout << "havac obj constructed, generating hmm seq pair" << std::endl;
  generateRandomHmmSeqPairToFiles(sequenceLength, seqFileSrc, phmmFileSrc);

  std::cout << "loading phmm to device" << std::endl;
  havac->loadPhmm(phmmFileSrc);

  std::cout << "loading sequence to device, " << std::endl;
  havac->loadSequence(seqFileSrc);
  std::cout << "beginning hardware run..." << std::endl;
  havac->runHardwareClient();

  std::cout << "hardware run finished, verifing hits" << std::endl;
  vector<HavacHit> hardwareHits = havac->getHitsFromFinishedRun();
  std::cout << "verification generated " << hardwareHits.size() << " hits." << std::endl;
  std::cout << "hit verification process finished, now checking validity" << std::endl;


  //run the software impl
  //load the fasta and phmm from the files
  FastaVector fastaVector;
  P7HmmList phmmList;
  fastaVectorInit(&fastaVector);
  FastaVectorReturnCode fastaRc = fastaVectorReadFasta(seqFileSrc, &fastaVector);
  P7HmmReturnCode phmmRc = readP7Hmm(phmmFileSrc, &phmmList);

  diagHitCheck(&fastaVector, &phmmList, 0.02f);


  std::cout << "fastavector and p7hmm generated, running soft ssv" << std::endl;
  shared_ptr<vector<ReferenceSsvHit>> softwareHits = HitsFromSsv(&fastaVector,
    &phmmList, desiredPValue);

  std::cout << "soft ssv finished, found " << softwareHits->size() << "hits. comparing hits..." << std::endl;

  compareHardwareSoftwareHits(hardwareHits, softwareHits);

  std::cout << "finished comparison" << std::endl;

  p7HmmListDealloc(&phmmList);
  fastaVectorDealloc(&fastaVector);
}



void compareHardwareSoftwareHits(vector<HavacHit> hardwareHits, shared_ptr<vector<ReferenceSsvHit>> softwareHits) {
  {
    bool foundAllHwHits = true;
    for (auto & hwHit : hardwareHits) {
      bool hwHitFoundinSwHits = false;

      for (ReferenceSsvHit& swHit : *softwareHits) {
        if (hwHit.phmmIndex == swHit.phmmNumber &&
          hwHit.phmmPosition == swHit.phmmPosition &&
          hwHit.sequenceIndex == swHit.sequenceNumber &&
          hwHit.sequencePosition == swHit.sequencePosition) {
          hwHitFoundinSwHits = true;
        }
      }
      if (!hwHitFoundinSwHits) {
        foundAllHwHits = false;
        std::cout << "error: hw hit not found in software hits: " << hwHit.toString() << std::endl;
      }
    }
    if (foundAllHwHits) {
      std::cout << "found all " << hardwareHits.size() << "hardware hits" << std::endl;
    }
  }
  {
    bool foundAllSwHits;
    for (ReferenceSsvHit& swHit : *softwareHits) {
      bool swHitFoundInHwHits = false;
      for (auto & hwHit : hardwareHits) {
        if (hwHit.phmmIndex == swHit.phmmNumber &&
          hwHit.phmmPosition == swHit.phmmPosition &&
          hwHit.sequenceIndex == swHit.sequenceNumber &&
          hwHit.sequencePosition == swHit.sequencePosition) {
          swHitFoundInHwHits = true;
        }
      }
      if (!swHitFoundInHwHits) {
        foundAllSwHits = false;
        std::cout << "error: sw hit not found in hardware hits: " << swHit.toString() << std::endl;
      }
    }
    if (foundAllSwHits) {
      std::cout << "found all " << hardwareHits.size() << "hardware hits" << std::endl;
    }
  }
}
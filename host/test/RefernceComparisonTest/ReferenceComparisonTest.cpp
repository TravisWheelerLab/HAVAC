#include "ReferenceComparisonTest.hpp"
#include "../../Havac.hpp"
#include "Ssv.hpp"
#include "../../../test/generator/hmmSeqGenerator.h"
#include "../../../device/PublicDefines.h"

extern "C"{
  #include <FastaVector.h>
  #include <p7HmmReader.h>
}


const char* phmmFileSrc = "refTest.hmm";
const char* seqFileSrc = "refTest.fasta";
const float desiredPValue = 0.05f;
const uint32_t sequenceLength = NUM_CELL_PROCESSORS * 5;

int main(int argc, char** argv) {
  //run the hardware
  shared_ptr<Havac> havac = make_shared<Havac>(0, desiredPValue);
  generateRandomHmmSeqPairToFiles(sequenceLength, seqFileSrc, phmmFileSrc);
  struct HmmSeqPair generateRandomHmmSeqPair(const uint32_t seqLength);
  havac->loadPhmm(phmmFileSrc);
  havac->loadSequence(seqFileSrc);
  havac->runHardwareClient();
  shared_ptr<vector<VerifiedHit>> hardwareHits = havac->getHitsFromFinishedRun();


  //run the software impl
  //load the fasta and phmm from the files
  FastaVector fastaVector;
  P7HmmList phmmList;
  FastaVectorReturnCode fastaRc = fastaVectorReadFasta(seqFileSrc, &fastaVector);
  P7HmmReturnCode phmmRc = readP7Hmm(phmmFileSrc, &phmmList);

  shared_ptr<vector<ReferenceSsvHit>> softwareHits = HitsFromSsv(&fastaVector,
    &phmmList, desiredPValue);


  compareHardwareSoftwareHits(hardwareHits, softwareHits);
}



void compareHardwareSoftwareHits(shared_ptr<vector<VerifiedHit>> hardwareHits, shared_ptr<vector<ReferenceSsvHit>> softwareHits) {
  {
    bool foundAllHwHits = true;
    for (VerifiedHit& hwHit : *hardwareHits) {
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
        std::cout << "error: hit not found in software hits: " << hwHit.toString() << std::endl;
      }
    }
    if (foundAllHwHits) {
      std::cout << "found all " << hardwareHits->size() << "hardware hits" << std::endl;
    }
  }
  {
    bool foundAllSwHits;
    for (ReferenceSsvHit& swHit : *softwareHits) {
      bool swHitFoundInHwHits = false;
      for (VerifiedHit& hwHit : *hardwareHits) {
        if (hwHit.phmmIndex == swHit.phmmNumber &&
          hwHit.phmmPosition == swHit.phmmPosition &&
          hwHit.sequenceIndex == swHit.sequenceNumber &&
          hwHit.sequencePosition == swHit.sequencePosition) {
          swHitFoundInHwHits = true;
        }
      }
      if (!swHitFoundInHwHits) {
        foundAllSwHits = false;
        std::cout << "error: hit not found in software hits: " << swHit.toString() << std::endl;
      }
    }
    if (foundAllSwHits) {
      std::cout << "found all " << hardwareHits->size() << "hardware hits" << std::endl;
    }
  }
}
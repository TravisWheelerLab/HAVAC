#include "multiInputTest.hpp"
#include "../../Havac.hpp"
#include "../../../test/generator/hmmSeqGenerator.h"
#include "../../../device/PublicDefines.h"
#include "../../../PhmmReprojection/PhmmReprojection.h"
#include <stdexcept>
#include <chrono>

extern "C" {
  #include <FastaVector.h>
  #include <p7HmmReader.h>
}

const char* phmmFileSrc = "refTest.hmm";
const char* seqFileSrc = "refTest.fasta";
const float desiredPValue = 0.05f;
const int64_t MAX_ALLOWABLE_HIT_POSITIONAL_DIFFERENCE = 16;

void makeSequencesAndPhmms(uint32_t sequenceLength, uint32_t numSequences, const float subProbability);
void fileAppend(const char* fileReadFrom, const char* fileWriteTo);

const float SUB_PROBABILITY = 0.2f;
//args:
//0: program name
//1: xclbin src
//2: num sequences
//3: sequence length
int main(int argc, char** argv) {
  //initialize the RNG seed
  srand(time(NULL));

  uint32_t numSequences;
  uint32_t sequenceLength = 0;
  std::string xclbinSrc;

  //check to see if abort was requested.
  if (argc == 2) {
    std::cout << "1 argument given, assuming phmm/fasta files already exist, this tool will not generate them." << std::endl;
    xclbinSrc = argv[1];
  }
  else if (argc < 4) {
    std::cout << "error: program requires an argument for the source for the xclbin src, the number of sequences, "
      "and the length of those sequences" << std::endl;
    exit(2);
  }
  else {
    xclbinSrc = argv[1];
    try {
      numSequences = std::stoi(argv[2]);
    }
    catch (std::invalid_argument& e) {
      std::stringstream ss;
      ss << "was not able to parse integer value from the num sequences argument (value given = '" <<
        argv[2] << "')" << std::endl;
      throw std::runtime_error(ss.str());
    }
    try {
      sequenceLength = std::stoi(argv[3]);
    }
    catch (std::invalid_argument& e) {
      std::stringstream ss;
      ss << "was not able to parse integer value from the sequence length argument (value given = '" <<
        argv[3] << "')" << std::endl;
      throw std::runtime_error(ss.str());
    }

    std::cout << "generating " << numSequences << " sequences, and their corresponding phmms" << std::endl;
    makeSequencesAndPhmms(sequenceLength, numSequences, SUB_PROBABILITY);
  }
  //run the hardware


  std::cout << "generating havac object" << std::endl;
  auto havacCtorStartTime = std::chrono::high_resolution_clock::now();
  shared_ptr<Havac> havac = make_shared<Havac>(0, desiredPValue, xclbinSrc);
  auto havacCtorEndTime = std::chrono::high_resolution_clock::now();

  auto dataLoadStartTime = std::chrono::high_resolution_clock::now();

  havac->loadPhmm(phmmFileSrc);
  havac->loadSequence(seqFileSrc);

  auto dataLoadEndTime = std::chrono::high_resolution_clock::now();
  std::cout << "beginning hardware run..." << std::endl;

  auto clientRunStartTime = std::chrono::high_resolution_clock::now();
  havac->runHardwareClient();
  auto clientRunEndTime = std::chrono::high_resolution_clock::now();


  auto hitValidationStartTime = std::chrono::high_resolution_clock::now();
  shared_ptr<vector<VerifiedHit>> hardwareHits = havac->getHitsFromFinishedRun();
  auto hitValidationEndTime = std::chrono::high_resolution_clock::now();

  auto havacCtorDuration = std::chrono::duration_cast<std::chrono::microseconds>(havacCtorEndTime - havacCtorStartTime);
  auto dataLoadDuration = std::chrono::duration_cast<std::chrono::microseconds>(dataLoadEndTime - dataLoadStartTime);
  auto clientRunDuration = std::chrono::duration_cast<std::chrono::microseconds>(clientRunEndTime - clientRunStartTime);
  auto hitValidationDuration = std::chrono::duration_cast<std::chrono::microseconds>(hitValidationEndTime - hitValidationStartTime);
  auto totalDuration = std::chrono::duration_cast<std::chrono::microseconds>(hitValidationEndTime - havacCtorStartTime);
  std::cout << "havac object construction took " << havacCtorDuration.count() << " microseconds" << std::endl;
  std::cout << "data load took " << dataLoadDuration.count() << " microseconds." << std::endl;
  std::cout << "hardware client run took " << clientRunDuration.count() << " microseconds" << std::endl;
  std::cout << "hit validation took " << hitValidationDuration.count() << " microseconds" << std::endl;
  std::cout << "total runtime: " << totalDuration.count() << " microseconds." << std::endl;

  std::cout << "hardware run finished, verifing hits" << std::endl;
  // std::cout << "verification generated " << hardwareHits->size() << " hits." << std::endl;


  // //run the software impl
  // //load the fasta and phmm from the files
  // FastaVector fastaVector;
  // P7HmmList phmmList;
  // fastaVectorInit(&fastaVector);
  // FastaVectorReturnCode fastaRc = fastaVectorReadFasta(seqFileSrc, &fastaVector);
  // P7HmmReturnCode phmmRc = readP7Hmm(phmmFileSrc, &phmmList);


  // std::cout << "fastavector and p7hmm generated, running soft ssv" << std::endl;
  // shared_ptr<vector<ReferenceSsvHit>> softwareHits = HitsFromSsv(&fastaVector,
  //   &phmmList, desiredPValue);


  // if(softwareHits->size() != hardwareHits->size()){
  //   std::cout << "Error: software hits size "<< softwareHits->size()<< " did not match hardware hits size "<< hardwareHits->size()<< std::endl;
  // }
  // std::cout << "soft ssv finished, found " << softwareHits->size() << "hits. comparing hits..." << std::endl;
  // for (size_t i = 0; i < softwareHits->size();i++) {
  //   std::cout << softwareHits->at(i).toString() << std::endl;
  // }


  // compareHardwareSoftwareHits(hardwareHits, softwareHits, &phmmList, &fastaVector);

  // std::cout << "finished comparison" << std::endl;

  // p7HmmListDealloc(&phmmList);
  // fastaVectorDealloc(&fastaVector);
}



void compareHardwareSoftwareHits(shared_ptr<vector<VerifiedHit>> hardwareHits, shared_ptr<vector<ReferenceSsvHit>> softwareHits,
  P7HmmList* phmmList, FastaVector* fastaVector) {
    {
      bool foundAllHwHits = true;
      for (VerifiedHit& hwHit : *hardwareHits) {
        bool hwHitFoundinSwHits = false;
        bool foundApproximateHitMatch = false;

        for (ReferenceSsvHit& swHit : *softwareHits) {
          if (hwHit.phmmIndex == swHit.phmmNumber &&
            hwHit.phmmPosition == swHit.phmmPosition &&
            hwHit.sequenceIndex == swHit.sequenceNumber &&
            hwHit.sequencePosition == swHit.sequencePosition) {
            hwHitFoundinSwHits = true;
          }
          else if (hwHit.phmmIndex == swHit.phmmNumber &&
            hwHit.sequenceIndex == swHit.sequenceNumber) {
            int64_t hwHitOnDiagonal = (int64_t)hwHit.phmmPosition - (int64_t)hwHit.sequencePosition;
            int64_t swHitOnDiagonal = (int64_t)swHit.phmmPosition - (int64_t)swHit.sequencePosition;
            int64_t phmmPositionDifference = (int64_t)hwHit.phmmPosition - swHit.phmmPosition;
            if (hwHitOnDiagonal == swHitOnDiagonal && abs(phmmPositionDifference) <= MAX_ALLOWABLE_HIT_POSITIONAL_DIFFERENCE) {
              foundApproximateHitMatch = true;
            }
          }
        }
        if (!hwHitFoundinSwHits && !foundApproximateHitMatch) {
          bool isExplainable = hardwareHitIsExplainable(hwHit, phmmList, fastaVector);
          if (isExplainable) {
            std::cout << "hw hit was not found in sw hit, but recreating the hit using the software ssv was successful." << std::endl;
          }
          else {
            foundAllHwHits = false;
            std::cout << "\tERROR: hw hit not found in software hits: " << hwHit.toString() << std::endl;
          }
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
        bool foundApproximateHitMatch = false;
        for (VerifiedHit& hwHit : *hardwareHits) {
          if (hwHit.phmmIndex == swHit.phmmNumber &&
            hwHit.phmmPosition == swHit.phmmPosition &&
            hwHit.sequenceIndex == swHit.sequenceNumber &&
            hwHit.sequencePosition == swHit.sequencePosition) {
            swHitFoundInHwHits = true;
          }
          else if (hwHit.phmmIndex == swHit.phmmNumber &&
            hwHit.sequenceIndex == swHit.sequenceNumber) {
            int64_t hwHitOnDiagonal = (int64_t)hwHit.phmmPosition - (int64_t)hwHit.sequencePosition;
            int64_t swHitOnDiagonal = (int64_t)swHit.phmmPosition - (int64_t)swHit.sequencePosition;
            int64_t phmmPositionDifference = (int64_t)hwHit.phmmPosition - swHit.phmmPosition;
            if (hwHitOnDiagonal == swHitOnDiagonal && abs(phmmPositionDifference) <= MAX_ALLOWABLE_HIT_POSITIONAL_DIFFERENCE) {
              foundApproximateHitMatch = true;
            }
          }
        }
        if (!swHitFoundInHwHits && !foundApproximateHitMatch) {
          foundAllSwHits = false;
          std::cout << "error: sw hit not found in hardware hits: " << swHit.toString() << std::endl;
        }
      }
      if (foundAllSwHits) {
        std::cout << "found all " << hardwareHits->size() << "hardware hits" << std::endl;
      }
    }
}


void fileAppend(const char* fileReadFrom, const char* fileWriteTo) {
  FILE* inputFile = fopen(fileReadFrom, "r");
  if (!inputFile) {
    std::stringstream ss;
    ss << "could not open file " << fileReadFrom << "for reading";
    throw std::runtime_error(ss.str());
  }
  FILE* outputFile = fopen(fileWriteTo, "a+");
  if (!outputFile) {
    fclose(inputFile);
    std::stringstream ss;
    ss << "error opening file " << fileWriteTo << "for appending";
    throw std::runtime_error(ss.str());
  }
  //write a newline to make sure the entries are seperated
  while (true) {
    char c = fgetc(inputFile);
    if (feof(inputFile)) {
      break;
    }
    fputc(c, outputFile);
  }

  fputc('\n', outputFile);
  fclose(inputFile);
  fclose(outputFile);
}

void makeSequencesAndPhmms(uint32_t sequenceLength, uint32_t numSequences, const float subProbability) {
  std::remove(seqFileSrc);
  std::remove(phmmFileSrc);
  std::cout << "generating " << numSequences << " sequences, and their corresponding phmms" << std::endl;
  const char* tmpSeqFileSrc = "_tmp.fasta";
  const char* tmpPhmmFileSrc = "_tmp.phmm";
  std::remove(seqFileSrc);
  std::remove(tmpPhmmFileSrc);

  for (uint32_t i = 0; i < numSequences; i++) {
    generateRandomHmmSeqPairToFiles(sequenceLength, tmpSeqFileSrc, tmpPhmmFileSrc, false, subProbability);


    fileAppend(tmpSeqFileSrc, seqFileSrc);
    fileAppend(tmpPhmmFileSrc, phmmFileSrc);

    std::remove(tmpSeqFileSrc);
    std::remove(tmpPhmmFileSrc);
  }
}


bool hardwareHitIsExplainable(VerifiedHit& hwHit, P7HmmList* phmmList, FastaVector* fastaVector) {
  std::cout << "attempting to explain hw hit:" << hwHit.toString() << std::endl;
  uint32_t walkbackMaxLength = std::min((uint64_t)hwHit.phmmPosition, hwHit.sequencePosition);

  uint32_t sequenceStartingPosition = hwHit.sequenceIndex == 0 ? 0 : fastaVector->metadata.data[hwHit.sequenceIndex - 1].sequenceEndPosition;

  const float  scoreMultiplier = findThreshold256ScalingFactor(&phmmList->phmms[hwHit.phmmIndex], desiredPValue);
  int16_t scoreAccumulator = 0;
  for (uint32_t walkbackStep = 0; walkbackStep <= walkbackMaxLength; walkbackStep++) {
    char symbolFromSequence = fastaVector->sequence.charData[sequenceStartingPosition + hwHit.sequencePosition - walkbackStep];
    uint8_t encodedSymbol;
    switch (symbolFromSequence) {
    case 'a':case 'A': encodedSymbol = 0; break;
    case 'c':case 'C': encodedSymbol = 1; break;
    case 'g':case 'G': encodedSymbol = 2; break;
    default: encodedSymbol = 3;
    }

    float floatMatchScore = phmmList->phmms[hwHit.phmmIndex].model.matchEmissionScores[(4 * (hwHit.phmmPosition - walkbackStep)) + encodedSymbol];
    float projectedMatchScore = emissionScoreToProjectedScore(floatMatchScore, scoreMultiplier);
    int8_t projectedMatchScoreAsInt = projectedMatchScore;

    // std::cout <<symbolFromSequence<< "~ "<<(int)scoreAccumulator<< ", "<< (int)projectedMatchScoreAsInt<< " = "<< scoreAccumulator+projectedMatchScoreAsInt<<
    // "\ts="<<sequenceStartingPosition+hwHit.sequencePosition-walkbackStep<< ", p="<<hwHit.phmmPosition-walkbackStep<<std::endl;

    scoreAccumulator += projectedMatchScoreAsInt;
    if (scoreAccumulator >= 256) {
      return true;
    }
    else if (scoreAccumulator <= 0) {
      return false;
    }

  }
  return false;

}
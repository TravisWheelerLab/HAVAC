#include "../../host/Havac.hpp"
#include "../../PhmmReprojection/PhmmReprojection.h"
#include "hmmerHit.hpp"
#include <iostream>
#include <vector>
#include <memory>
#include <sstream>
#include <string>

using std::shared_ptr;
using std::vector;
using std::make_shared;
using std::string;

struct HmmerWindow {
  uint32_t sequenceIndex;
  uint32_t windowStart;
  uint32_t windowEnd;
  uint32_t windowLength;
  uint32_t phmmPosition;
  char accession[1024];
};


using std::string;
using std::vector;


vector<HmmerWindow> makeHmmerWindowListFromFile(string hmmerWindowFileSrc);
void compareHitsToWindows(shared_ptr<vector<VerifiedHit>> hardwareHits, vector<HmmerWindow> hmmerWindowList, P7HmmList* phmmList);
void compareHavacHitsToHmmerHits(shared_ptr<vector<VerifiedHit>>& hardwareHits, vector<HmmerHitLine>& hmmerHits, P7HmmList *phmmList);

//args:
//xclbin file src
//fasta src
//phmm src
//hit output file from hmmer
int main(int argc, char** argv) {
  if (argc != 5) {
    std::cout << "ERROR: program requires 4 program arguments.\n\txclbin file src\n\tfasta src\n\tphmm src\n\thmmer hit output file src" << std::endl;
    exit(-1);
  }

  const float desiredPValue = 0.02f;
  string xclbinSrc = argv[1];
  string fastaFileSrc = argv[2];
  string phmmFileSrc = argv[3];
  string hmmerWindowFileSrc = argv[4];

  // vector<HmmerWindow> hmmerWindowList = makeHmmerWindowListFromFile(hmmerWindowFileSrc);

  shared_ptr<Havac> havac = make_shared<Havac>(0, desiredPValue, xclbinSrc);
  havac->loadPhmm(phmmFileSrc);
  havac->loadSequence(fastaFileSrc);
  havac->runHardwareClient();
  P7HmmList phmmList;
  readP7Hmm(phmmFileSrc.c_str(), &phmmList);

  uint32_t totalPhmmLength = 0;
for(uint32_t i = 0; i < phmmList.count;i++){
  totalPhmmLength += phmmList.phmms[i].header.modelLength;
}

std::cout << "total concat model lengths: "<< totalPhmmLength<<std::endl;

  shared_ptr<vector<VerifiedHit>> hardwareHits = havac->getHitsFromFinishedRun();
  // std::cout << "# windows (" << hmmerWindowList.size() << ")" << std::endl;
  std::cout << "# hardware hits (" << hardwareHits->size() << "):" << std::endl;

  // compareHitsToWindows(hardwareHits, hmmerWindowList, &phmmList);

  vector<HmmerHitLine> hmmerHits = getHitsFromFile(hmmerWindowFileSrc);
  compareHavacHitsToHmmerHits(hardwareHits, hmmerHits, &phmmList);
}


void compareHavacHitsToHmmerHits(shared_ptr<vector<VerifiedHit>>& hardwareHits, vector<HmmerHitLine>& hmmerHits, P7HmmList* phmmList) {
  uint32_t numHardwareHitsFoundInHmmerHits = 0;
  uint32_t numHmmerHitsFoundInHardwareHits = 0;
  
  for(uint32_t hardwareHitIndex = 0; hardwareHitIndex < hardwareHits->size(); hardwareHitIndex++){
    bool hitFoundForThisHwHit = false;
      VerifiedHit *hardwareHit = &hardwareHits->at(hardwareHitIndex);
      string hardwareAccession = phmmList->phmms[hardwareHit->phmmIndex].header.accessionNumber;

    for(uint32_t hmmerHitIndex = 0; hmmerHitIndex < hmmerHits.size(); hmmerHitIndex++){
      HmmerHitLine *hitLine = &hmmerHits[hmmerHitIndex];
      bool matchesAccession = hitLine->accession.compare(hardwareAccession) == 0;
      bool seqInRange = hardwareHit->sequencePosition >= hitLine->envfrom && hardwareHit->sequencePosition <= hitLine->envto;
      bool hmmInRange = hardwareHit->phmmPosition >= hitLine->hmmfrom && hardwareHit->phmmPosition <= hitLine->hmmto;

      bool hitsMatch = matchesAccession && seqInRange && hmmInRange;

      if(hitsMatch){
        hitFoundForThisHwHit = true;
        numHardwareHitsFoundInHmmerHits++;
        break;
      }
    }
    if (!hitFoundForThisHwHit){
      // std::cout << "ERROR: Verified Hit "<< hardwareHit->toString()<< " could not find a matching hmmer hit"<< std::endl;
    }
  }


  for(uint32_t hmmerHitIndex  = 0; hmmerHitIndex < hmmerHits.size(); hmmerHitIndex++){
    bool hitFoundForThisHmmerHit = false;
    HmmerHitLine* hitLine = &hmmerHits[hmmerHitIndex];

    for (uint32_t hardwareHitIndex = 0; hardwareHitIndex < hardwareHits->size(); hardwareHitIndex++) {
      VerifiedHit* hardwareHit = &hardwareHits->at(hardwareHitIndex);
      string hardwareAccession = phmmList->phmms[hardwareHit->phmmIndex].header.accessionNumber;

      bool matchesAccession = hitLine->accession.compare(hardwareAccession) == 0;
      bool seqInRange = hardwareHit->sequencePosition >= hitLine->envfrom && hardwareHit->sequencePosition <= hitLine->envto;
      bool hmmInRange = hardwareHit->phmmPosition >= hitLine->hmmfrom && hardwareHit->phmmPosition <= hitLine->hmmto;

      bool hitsMatch = matchesAccession && seqInRange && hmmInRange;

      if(hitsMatch){
        hitFoundForThisHmmerHit = true;
        numHmmerHitsFoundInHardwareHits++;
        break;
      }
    }
    if(!hitFoundForThisHmmerHit){
      std::cout << "ERROR: Hmmer hit "<< hitLine->toString() << " could not find a matching hardware hit"<<std::endl;
    }
  }
  std::cout << "final result: " << numHardwareHitsFoundInHmmerHits<<"/"<< hardwareHits->size()<< "hardware hits located, "<<
  numHmmerHitsFoundInHardwareHits<< "/"<< hmmerHits.size()<< "hmmer hits located"<< std::endl;
}



struct RefSsvHit {
  uint32_t seqPos;
  uint32_t phmmPos;
};


uint8_t encodeSequenceSymbol(char symbol) {
  switch (symbol) {
  case 'a': case 'A': return 0;
  case 'c': case 'C': return 1;
  case 'g': case 'G': return 2;
  case 't': case 'T': return 3;
  default: return rand() % 4;
  }
}
void checkWindowsAgainstReferenceSsv(vector<HmmerWindow>& windowList, FastaVector* fastaVector, P7HmmList* phmmList, const float pValue) {
  vector<RefSsvHit> ssvHits;
  P7Hmm* phmm = &phmmList->phmms[0];
  const uint32_t sequenceLength = fastaVector->sequence.count;
  const uint32_t modelLength = phmm->header.modelLength;
  // const float scoreMultiplier = findThreshold256ScalingFactor(&phmmList->phmms[0],pValue);
  vector<uint8_t> cellScoreList(modelLength * 4);

  for (uint32_t i = 0; i < modelLength * 4; i++) {
    cellScoreList[i] = 0;
  }

  vector<int8_t> projectedScores(modelLength);
  p7HmmProjectForThreshold256(phmm, pValue, projectedScores.data());

  for (uint32_t seqPosition = 0; seqPosition < sequenceLength; seqPosition++) {
    for (uint32_t phmmPosition = modelLength - 1; phmmPosition > 0; phmmPosition++) {
      char symbol = fastaVector->sequence.charData[seqPosition];
      uint8_t encodedSymbol = encodeSequenceSymbol(symbol);
      int16_t matchScore = cellScoreList[seqPosition + encodedSymbol];
      int16_t result = (int16_t)cellScoreList[phmmPosition - 1] + matchScore;
      if (result >= 256) {
        RefSsvHit hit;
        hit.phmmPos = phmmPosition;
        hit.seqPos = seqPosition;
        ssvHits.push_back(hit);
      }
      if (result >= 256 || result < 0) {
        result = 0;
      }

      cellScoreList[phmmPosition] = result;
    }
  }
  char symbol = fastaVector->sequence.charData[0];
  uint8_t encodedSymbol = encodeSequenceSymbol(symbol);
  cellScoreList[0] = cellScoreList[encodedSymbol];
}


vector<HmmerWindow> makeHmmerWindowListFromFile(string hmmerWindowFileSrc) {
  vector<HmmerWindow> hmmerWindowList;
  FILE* file = fopen(hmmerWindowFileSrc.c_str(), "r");
  if (!file) {
    std::cout << "ERROR: could not open hmmer window file at source " << hmmerWindowFileSrc << std::endl;
    exit(-2);
  }
  //read the header
  char buffer[1024];
  fgets(buffer, 1024, file);

  while (!feof(file)) {
    HmmerWindow thisWindow;
    fgets(buffer, 2048, file);
    int numScanned = sscanf(buffer, "%i\t%i\t%i\t%i\t%i\t%s\t\n", &thisWindow.sequenceIndex, &thisWindow.windowStart, &thisWindow.windowEnd,
      &thisWindow.windowLength, &thisWindow.phmmPosition, thisWindow.accession);
    if (numScanned < 6) {
      std::cout << "could not parse line, found " << numScanned << "entries on line: \"" << buffer << std::endl;
    }
    else {
      hmmerWindowList.push_back(thisWindow);
    }
  }
  fclose(file);
  return hmmerWindowList;
}


bool hitMatchesWindow(VerifiedHit hit, char* phmmAccession, HmmerWindow window) {
  // uint32_t sequenceWindowStartPosition = window.windowStart;
  // uint32_t sequenceWindowEndPosition = window.windowEnd;
  bool matchesSequenceIndex = window.sequenceIndex == hit.sequenceIndex;
  bool sequencePositionInsideWindow = hit.sequencePosition >= window.windowStart &&
    hit.sequencePosition <= window.windowEnd;
  // bool phmmPositionInsideWindow = hit.phmmPosition <= window.phmmPosition &&
  //   hit.phmmPosition >= (window.phmmPosition - window.windowLength);
  bool matchesAccession = strcmp(phmmAccession, window.accession) == 0;
  bool hitInsideWindow = matchesSequenceIndex && sequencePositionInsideWindow;
  if (hitInsideWindow && !matchesAccession) {
    printf("hit @ hw pos %zu occurs withing window [%u,%u], but accessions dont match (hw acc %s, hmmer %s)\n",
      hit.sequencePosition, window.windowStart, window.windowEnd, phmmAccession, window.accession);
  }

  bool fullMatch = hitInsideWindow && matchesAccession;
  return fullMatch;
}

void compareHitsToWindows(shared_ptr<vector<VerifiedHit>> hardwareHits, vector<HmmerWindow> hmmerWindowList, P7HmmList* phmmList) {
  uint32_t numHwHitsFoundInHmmerHits = 0;
  uint32_t numHmmerHitsFoundInHwHits = 0;
  for (size_t hardwareHitIndex = 0; hardwareHitIndex < hardwareHits->size();hardwareHitIndex++) {
    VerifiedHit hardwareHit = hardwareHits->at(hardwareHitIndex);
    bool windowFoundForHit = false;
    char* phmmAccession = phmmList->phmms[hardwareHit.phmmIndex].header.accessionNumber;

    for (size_t windowIndex = 0; windowIndex < hmmerWindowList.size(); windowIndex++) {
      HmmerWindow window = hmmerWindowList.at(windowIndex);
      if (hitMatchesWindow(hardwareHit, phmmAccession, window)) {
        windowFoundForHit = true;
        // std::cout << "hit index " << hardwareHitIndex << "matches window index " << windowIndex << std::endl;
        break;
      }
    }
    if (windowFoundForHit) {
      numHwHitsFoundInHmmerHits++;
    }
    else {
      std::cout << "FAIL: hw hit #" << hardwareHitIndex << " was not found to be inside a window. hw hit: " << hardwareHit.toString() << std::endl;
    }
  }

  for (size_t windowIndex = 0; windowIndex < hmmerWindowList.size(); windowIndex++) {
    bool hitFoundInsideWindow = false;
    HmmerWindow window = hmmerWindowList.at(windowIndex);

    for (size_t hardwareHitIndex = 0; hardwareHitIndex < hardwareHits->size();hardwareHitIndex++) {
      VerifiedHit hardwareHit = hardwareHits->at(hardwareHitIndex);
      char* phmmAccession = phmmList->phmms[hardwareHit.phmmIndex].header.accessionNumber;

      if (hitMatchesWindow(hardwareHit, phmmAccession, window)) {
        hitFoundInsideWindow = true;
        // std::cout << "window index " << windowIndex << "matches hit index " << hardwareHitIndex << std::endl;
        break;
      }
    }

    if (hitFoundInsideWindow) {
      numHmmerHitsFoundInHwHits++;
    }
    else {
      std::cout << "FAIL: window did not seem to contain a hit: window #" << windowIndex << std::endl;
    }
  }

  std::cout << "final results: \n\t "<< numHwHitsFoundInHmmerHits<<"/"<< hardwareHits->size()<< " hardware hits found matching windows"<<std::endl;
  std::cout << numHmmerHitsFoundInHwHits<<"/"<< hmmerWindowList.size()<< " hmmer windows found matching hardware hits."<<std::endl;
}
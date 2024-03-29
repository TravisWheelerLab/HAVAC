#include<vector>
#include<memory>
#include<cstdint>
#include<iostream>
#include<sstream>
#include<string.h>
#include<fstream>
#include<cstring>
#include<cmath>
#include "../../PhmmReprojection/PhmmReprojection.h"
extern "C" {
  #include <FastaVector.h>
  #include <p7HmmReader.h>
}


/*
  This file is used to test the hits generated by hmmer, against either a reference SSV, or the hits generated by havac.

*/

using std::vector;
using std::shared_ptr;
using std::make_shared;
using std::string;


struct TracebackData {
  uint32_t tracebackSteps;
  int16_t maximumScore;
};

struct HmmerWindow {
  uint32_t sequenceId;
  uint32_t windowStart;
  uint32_t windowEnd;
  uint32_t windowLength;
  uint32_t hmmPosition;
  string accession;
};



vector<HmmerWindow> readHmmerHits(string hmmerHitsFileSrc);
uint8_t encodeSymbol(const char symbol);
vector<float> getEmissionsFromHmmFile(string hmmFileSrc);
string getSequenceSectionWithEasel(uint32_t seqStartingPosition, uint32_t seqEndingPosition, string fastaFileSrc,
  string sequenceHeader, string tmpOutputFileSrc);
void printHitData(uint32_t sequencePosition, uint32_t phmmPosition, TracebackData tracebackData, int16_t cellResult);

bool symbolIsAmbiguous(char symbol);
void compareHmmHitSequencesWithFastaVector(vector<HmmerWindow>& hmmerWindows, char* fastaFileSrc, char* tmpOutputFileSrc);
// void compareHmmerHitsWithRefSsv(vector<HmmerWindow>& hmmerHits, FastaVector* fastaVector, P7Hmm* phmm, int8_t* const emissionScores);
void checkHmmerHitsAgainstRefSsv(vector<HmmerWindow>& hmmerHits, vector<int8_t>& emissionScoresVector,
  vector<float>& emissionsAsFloats, vector<float>& emissionsAsScaledFloats, char* fastaFileSrc,
  char* seqHeader, char* tmpOutputFileSrc, uint8_t hmmerThreshold);
uint16_t refSsv(string sequence, vector<int8_t>& emissionScoresVector);
float refSsvFloat(string sequence, vector<float>& emissionScoresVector);
uint16_t refSsvDiagonal(string sequence, vector<int8_t> emissionScoresVector, HmmerWindow& window);
float refSsvDiagonalFloat(string sequence, vector<float> emissionScoresVector, HmmerWindow& window);
uint16_t attemptTraceback(uint32_t windowStart, uint32_t windowEnd, uint32_t phmmPosition, vector<int8_t> emissionScores, FastaVector* fastaVector);
uint16_t attemptTracebackFloat(uint32_t windowStart, uint32_t windowEnd, uint32_t phmmPosition, vector<float> emissionScores, FastaVector* fastaVector);
vector<int8_t> deprojectAndRoundEmissions(vector<float>& emissionScores, const uint8_t thresholdScore);



uint8_t hmmerThresholdScore = 1;
int main(int argc, char** argv) {
  if (argc != 6) {
    printf("program requires 3 inputs.\n1. fasta file src\n2. hmm file src\n3. hmmerHit txt file\n4.p value for search\n");
    exit(1);
  }
  const uint8_t HMMER_THRESHOLD = 187;
  char* fastaFileSrc = argv[1];
  char* hmmFileSrc = argv[2];
  char* hmmerHitFileSrc = argv[3];
  float pValue;
  try {
    pValue = std::stof(argv[4]);
  }
  catch (std::invalid_argument& e) {
    std::cerr << "could not parse float for p-value string " << argv[4] << " (arg 4)." << std::endl;
    exit(2);
  }
  try {
    hmmerThresholdScore = std::stoi(argv[5]);
  }
  catch (std::invalid_argument& e) {
    std::cerr << "could not parse int value for hmmer threshold score " << argv[5] << " (arg 5)." << std::endl;
    exit(3);
  }

  P7HmmList phmmList;

  enum P7HmmReturnCode rc = readP7Hmm(hmmFileSrc, &phmmList);
  if (rc != p7HmmSuccess) {
    throw std::logic_error("reading the phmm file was not sucessful");
  }
  P7Hmm* phmm = &phmmList.phmms[0];

  const float emissionScalingFactor = findThreshold256ScalingFactor(phmm, pValue);

  std::cout << "generating logo emissions..." << std::endl;
  vector<float> hmmLogoEmissions = getEmissionsFromHmmFile(hmmFileSrc);
  vector<float> projectedLogoEmissions(hmmLogoEmissions.size());
  for (size_t i = 0; i < hmmLogoEmissions.size();i++) {
    projectedLogoEmissions[i] = hmmLogoEmissions[i] * emissionScalingFactor;
  }
  std::cout << "reading hmmer hits..." << std::endl;
  vector<HmmerWindow> hmmerHits = readHmmerHits(hmmerHitFileSrc);

  std::cout << "projecting scores..." << std::endl;
  vector<int8_t> emissionScoresVector(phmm->header.modelLength * 4);
  int8_t* emissionScoresPtr = emissionScoresVector.data();
  p7HmmProjectForThreshold256(phmm, pValue, emissionScoresPtr);
  std::cout << "projected scores" << std::endl;


  for (size_t i = 0; i < emissionScoresVector.size(); i++) {
    std::cout << (int)emissionScoresVector[i] << ", ";
  }
  std::cout << std::endl;

  char seqHeader[64];
  char tmpOutputFileSrc[64];
  strcpy(seqHeader, "chr22");
  strcpy(tmpOutputFileSrc, ".tmpSeqSegment.txt");
  checkHmmerHitsAgainstRefSsv(hmmerHits, emissionScoresVector, hmmLogoEmissions, projectedLogoEmissions,
    fastaFileSrc, seqHeader, tmpOutputFileSrc, HMMER_THRESHOLD);
  // compareHmmHitSequencesWithFastaVector(hmmerHits, fastaFileSrc, tmpOutputFileSrc);



}


void compareHmmHitSequencesWithFastaVector(vector<HmmerWindow>& hmmerWindows, char* fastaFileSrc, char* tmpOutputFileSrc) {
  std::cout << "in compare func" << std::endl;
  FastaVector fastaVector;
  std::cout << "initing fv" << std::endl;
  fastaVectorInit(&fastaVector);
  std::cout << "reading fasta at " << fastaFileSrc << std::endl;
  fastaVectorReadFasta(fastaFileSrc, &fastaVector);

  char seqHeader[32];
  strcpy(seqHeader, "chr22");

  for (HmmerWindow& window : hmmerWindows) {
    string sequenceSegment = getSequenceSectionWithEasel(window.windowStart, window.windowEnd, fastaFileSrc,
      seqHeader, tmpOutputFileSrc);

    char* fastaVectorSeqPtr;
    size_t seqLen;
    fastaVectorFastaGetSequence(&fastaVector, 0, &fastaVectorSeqPtr, &seqLen);
    string fastaVectorWindowString;
    fastaVectorWindowString.assign(fastaVectorSeqPtr + window.windowStart, window.windowLength);
    if (sequenceSegment.compare(fastaVectorWindowString)) {
      std::cout << "ERROR: substrings do not match: hmmer: " << sequenceSegment << ", havac: " << fastaVectorWindowString << std::endl;
    }
  }
  fastaVectorDealloc(&fastaVector);
}



void checkHmmerHitsAgainstRefSsv(vector<HmmerWindow>& hmmerHits, vector<int8_t>& emissionScoresVector,
  vector<float>& emissionsAsFloats, vector<float>& emissionsAsScaledFloats, char* fastaFileSrc,
  char* seqHeader, char* tmpOutputFileSrc, uint8_t hmmerThreshold) {
  uint16_t lowestU16Result = 256;
  float lowestFloatResult = 256;
  uint32_t numHits = hmmerHits.size();
  uint32_t numHitsThatPassReferenceSsv = 0;
  uint32_t numHitsThatPassFloatReference = 0;
  uint32_t numHitsThatPass250 = 0;
  uint32_t numHitsThatPass250Float = 0;

  FastaVector fastaVector;
  fastaVectorInit(&fastaVector);
  fastaVectorReadFasta(fastaFileSrc, &fastaVector);

  vector<int8_t> emissionsForHmmerThreshold = deprojectAndRoundEmissions(emissionsAsScaledFloats, hmmerThresholdScore);
  for (HmmerWindow& window : hmmerHits) {

    string sequenceSegment = getSequenceSectionWithEasel(window.windowStart, window.windowEnd, fastaFileSrc,
      seqHeader, tmpOutputFileSrc);

    // uint16_t largestScoreForWindow = refSsv(sequenceSegment, emissionScoresVector);
    // float scoreAsFloat = refSsvFloat(sequenceSegment, emissionsAsScaledFloats);
    // float scoreUsingLogo = refSsvFloat(sequenceSegment, emissionsAsFloats);

    uint16_t largestScoreForWindow = attemptTraceback(window.windowStart - 1, window.windowEnd - 1, window.hmmPosition + 1, emissionScoresVector, &fastaVector);
    // uint16_t largestScoreForWindow = attemptTraceback(window.windowStart - 1, window.windowEnd - 1,
    //   window.hmmPosition + 1, emissionsForHmmerThreshold, &fastaVector);
    float scoreAsFloat = attemptTracebackFloat(window.windowStart - 1, window.windowEnd - 1, window.hmmPosition + 1, emissionsAsScaledFloats, &fastaVector);
    float scoreUsingLogo = attemptTracebackFloat(window.windowStart - 1, window.windowEnd - 1, window.hmmPosition + 1, emissionsAsFloats, &fastaVector);


    // uint16_t largestScoreForWindow = refSsvDiagonal(sequenceSegment, emissionScoresVector, window);
    // float scoreAsFloat = refSsvDiagonalFloat(sequenceSegment, emissionsAsScaledFloats, window);
    // float scoreUsingLogo = refSsvDiagonalFloat(sequenceSegment, emissionsAsFloats, window);

    if (largestScoreForWindow > 250) {
      numHitsThatPass250++;
    }
    if (scoreAsFloat > 250) {
      numHitsThatPass250Float++;
    }
    if (scoreAsFloat >= 256) {
      numHitsThatPassFloatReference++;
    }
    else {
      lowestFloatResult = std::min(lowestFloatResult, scoreAsFloat);
    }
    if (largestScoreForWindow >= 256) {
      std::cout << "Hmmer window passes with threshold " << largestScoreForWindow << "." << std::endl;
      numHitsThatPassReferenceSsv++;
    }
    else {
      lowestU16Result = std::min(lowestU16Result, largestScoreForWindow);
      std::cout << "\033[31mHmmer window did not pass threshold, only got thresh\033[0m=" << largestScoreForWindow << std::endl;

      // float scoreAsFloat = refSsvFloat(sequenceSegment, emissionsAsScaledFloats);
      std::cout << "\t\x1B[32mfloat highest score " << scoreAsFloat << "\x1B[0m (using unprojected hmm logo values, " << scoreUsingLogo <<
        ")" << std::endl;
    }
  }
  std::cout << "could verify " << numHitsThatPassReferenceSsv << "/" << numHits << " hmmer windows." << std::endl;
  std::cout << "using uncompressed floats " << numHitsThatPassFloatReference << "/" << numHits << "hmmer windows." << std::endl;
  std::cout << "lowest 16-bit int result: " << lowestU16Result << ", lowest float " << lowestFloatResult << std::endl;
  std::cout << "pass 250 (int): " << numHitsThatPass250 << ", pass 250 (float): " << numHitsThatPass250Float << std::endl;
}


uint16_t refSsvDiagonal(string sequence, vector<int8_t> emissionScoresVector, HmmerWindow& window) {
  int16_t accumulatedScore = 0;
  int16_t highestScoreSeen = 0;
  uint32_t phmmStartIndex = window.hmmPosition - sequence.size();
  for (uint32_t i = 0; i < sequence.size();i++) {
    char symbol = sequence[i];
    if (symbolIsAmbiguous(symbol)) {
      std::cout << "SYMBOL @ window pos " << i << "(" << symbol << ") is ambiguous!" << std::endl;
    }
    uint8_t encodedSymbol = encodeSymbol(symbol);
    int8_t matchScore = emissionScoresVector.data()[(phmmStartIndex + i + 2) * 4 + encodedSymbol];
    accumulatedScore += matchScore;
    if (accumulatedScore > highestScoreSeen) {
      highestScoreSeen = accumulatedScore;
    }
    if (accumulatedScore < 0 || accumulatedScore > 256) {
      accumulatedScore = 0;
    }
  }
  return highestScoreSeen;
}

float refSsvDiagonalFloat(string sequence, vector<float> emissionScoresVector, HmmerWindow& window) {
  float accumulatedScore = 0;
  float highestScoreSeen = 0;
  uint32_t phmmStartIndex = window.hmmPosition - sequence.size();
  for (uint32_t i = 0; i < sequence.size();i++) {
    char symbol = sequence[i];
    if (symbolIsAmbiguous(symbol)) {
      std::cout << "SYMBOL @ window pos " << i << "(" << symbol << ") is ambiguous!" << std::endl;
    }
    uint8_t encodedSymbol = encodeSymbol(symbol);
    float matchScore = emissionScoresVector.data()[(phmmStartIndex + i + 2) * 4 + encodedSymbol];
    accumulatedScore += matchScore;
    if (accumulatedScore > highestScoreSeen) {
      highestScoreSeen = accumulatedScore;
    }
    if (accumulatedScore < 0 || accumulatedScore > 256) {
      accumulatedScore = 0;
    }
  }
  return highestScoreSeen;
}

uint16_t refSsv(string sequence, vector<int8_t>& emissionScoresVector) {
  uint16_t largestScoreSeen = 0;
  vector<uint8_t> cellScores(sequence.size());
  for (uint32_t i = 0; i < sequence.size();i++) {
    cellScores[i] = 0;
  }
  for (uint32_t vectorIndex = 0; vectorIndex < emissionScoresVector.size(); vectorIndex += 4) {
    int8_t* phmmVector = &emissionScoresVector.data()[vectorIndex];

    for (uint32_t seqIndex = sequence.size() - 1; seqIndex > 0; seqIndex--) {
      uint8_t encodedSymbol = encodeSymbol(sequence[seqIndex]);
      int8_t matchScore = phmmVector[encodedSymbol];

      int16_t result = (int16_t)cellScores[seqIndex - 1] + (int16_t)matchScore;

      if (result > largestScoreSeen) {
        largestScoreSeen = result;
      }

      if (result >= 256) {
        // std::cout << "vecindex " << vectorIndex / 4 << ", seqIndex " << seqIndex << ", seqlen " << sequence.size() << std::endl;
        // std::cout << (int16_t)cellScores[seqIndex - 1] << " + " << (int16_t)matchScore << " = " << result << std::endl;
        // for (int8_t i = 0; i < 4; i++) {
        //   std::cout << (int)phmmVector[i] << ", ";
        // }
        // std::cout << std::endl;

        result = 0;
      }
      else if (result < 0) {
        result = 0;
      }
      cellScores[seqIndex] = result;
    }
    int8_t firstCellMatchScore = phmmVector[encodeSymbol(sequence[0])];
    if (firstCellMatchScore > 0) {
      cellScores[0] = firstCellMatchScore;
    }
    else {
      cellScores[0] = 0;
    }
    if (cellScores[0] > largestScoreSeen) {
      largestScoreSeen = cellScores[0];
    }
  }

  return largestScoreSeen;
}


float refSsvFloat(string sequence, vector<float>& emissionScoresVector) {
  float largestScoreSeen = 0;
  vector<float> cellScores(sequence.size());
  for (uint32_t i = 0; i < sequence.size();i++) {
    cellScores[i] = 0;
  }
  for (uint32_t vectorIndex = 0; vectorIndex < emissionScoresVector.size(); vectorIndex += 4) {
    float* phmmVector = &emissionScoresVector.data()[vectorIndex];

    for (uint32_t seqIndex = sequence.size() - 1; seqIndex > 0; seqIndex--) {
      uint8_t encodedSymbol = encodeSymbol(sequence[seqIndex]);
      float matchScore = phmmVector[encodedSymbol];

      float result = cellScores[seqIndex - 1] + matchScore;

      if (result > largestScoreSeen) {
        largestScoreSeen = result;
      }

      if (result >= 256) {
        // std::cout << "vecindex " << vectorIndex / 4 << ", seqIndex " << seqIndex << ", seqlen " << sequence.size() << std::endl;
        // std::cout << (int16_t)cellScores[seqIndex - 1] << " + " << (int16_t)matchScore << " = " << result << std::endl;
        // for (int8_t i = 0; i < 4; i++) {
        //   std::cout << (int)phmmVector[i] << ", ";
        // }
        // std::cout << std::endl;

        result = 0;
      }
      else if (result < 0) {
        result = 0;
      }
      cellScores[seqIndex] = result;
    }
    float firstCellMatchScore = phmmVector[encodeSymbol(sequence[0])];
    if (firstCellMatchScore > 0) {
      cellScores[0] = firstCellMatchScore;
    }
    else {
      cellScores[0] = 0;
    }
    if (cellScores[0] > largestScoreSeen) {
      largestScoreSeen = cellScores[0];
    }
  }

  return largestScoreSeen;
}



vector<HmmerWindow> readHmmerHits(string hmmerHitsFileSrc) {
  vector<HmmerWindow> hmmerWindowList;
  FILE* file = fopen(hmmerHitsFileSrc.c_str(), "r");
  if (!file) {
    std::cout << "ERROR: could not open hmmer window file at source " << hmmerHitsFileSrc << std::endl;
    exit(-2);
  }
  //read the header
  char buffer[2048];
  fgets(buffer, 2048, file);
  char accession[128];

  while (!feof(file)) {
    uint32_t sIndex, winStart, winEnd, winLen, phmmPos;
    fgets(buffer, 2048, file);
    int numScanned = sscanf(buffer, "%i\t%i\t%i\t%i\t%i\t%s\t\n", &sIndex, &winStart,
      &winEnd, &winLen, &phmmPos, accession);
    if (numScanned < 6) {
      std::cout << "could not parse line, found " << numScanned << "entries on line: \"" << buffer << std::endl;
    }
    else {
      HmmerWindow window;
      window.sequenceId = sIndex;
      window.windowStart = winStart;
      window.windowEnd = winEnd;
      window.hmmPosition = phmmPos;
      window.windowLength = winLen;
      window.accession = accession;

      hmmerWindowList.push_back(window);
    }
  }
  fclose(file);
  return hmmerWindowList;
}


uint8_t encodeSymbol(const char symbol) {
  switch (symbol) {
  case 'a': case 'A': return 0;
  case 'c': case 'C': return 1;
  case 'g': case 'G': return 2;
  case 't': case 'T': return 3;
  default: return rand() % 4;
  }
}

void printHitData(uint32_t seqPosition, uint32_t phmmPosition, TracebackData tracebackData, int16_t cellResult) {
  std::cout << "hmmer hit @ i=" << seqPosition << ", j=" << phmmPosition << ", cellScore= " << cellResult <<
    ", maxTracebackScore = " << tracebackData.maximumScore << ", steps=" << tracebackData.tracebackSteps << "." << std::endl;
}



string getSequenceSectionWithEasel(uint32_t seqStartingPosition, uint32_t seqEndingPosition, string fastaFileSrc,
  string sequenceHeader, string tmpOutputFileSrc) {
  const uint32_t eachFlankExtraSequence = 1;
  std::stringstream ss;
  ss << "esl-sfetch -c " << seqStartingPosition - eachFlankExtraSequence << "." << seqEndingPosition << " ";
  ss << fastaFileSrc << " " << sequenceHeader;
  ss << " > " << tmpOutputFileSrc;
  system(ss.str().c_str());

  std::string sequenceSegment;
  std::ifstream fileInputStream;
  {
    std::string buffer;
    fileInputStream.open(tmpOutputFileSrc);
    //first getline gets the header
    std::getline(fileInputStream, buffer);

    while (buffer.size() != 0) {
      std::getline(fileInputStream, buffer);
      sequenceSegment.append(buffer);
    }
    std::getline(fileInputStream, sequenceSegment);
  }
  std::remove(tmpOutputFileSrc.c_str());

  return sequenceSegment;
}


vector<float> getEmissionsFromHmmFile(string hmmFileSrc) {
  vector<float> emissionScoresVector;
  const string hmmLogoTmpFileSrc = ".hmmLogoTmpFile.txt";
  std::stringstream ss;

  ss << "hmmlogo --height_score " << hmmFileSrc << " > " << hmmLogoTmpFileSrc;
  std::cout << "running " << ss.str() << std::endl;

  system(ss.str().c_str());
  std::cout << "hmm logo finished" << std::endl;
  std::string sequenceSegment;
  std::ifstream fileInputStream;
  fileInputStream.open(hmmLogoTmpFileSrc);

  string lineBuffer;
  //first, read the header "Residue Heights"
  std::getline(fileInputStream, lineBuffer);
  std::getline(fileInputStream, lineBuffer);
  do {
    std::cout << "comparing line buffer... " << lineBuffer << std::endl;
    float a_value, c_value, g_value, t_value;
    int modelPosition;
    int numScanned = std::sscanf(lineBuffer.c_str(), "%i:  %f  %f  %f  %f", &modelPosition, &a_value, &c_value, &g_value, &t_value);
    if (numScanned != 5) {
      std::cout << "scanned " << numScanned << " items" << std::endl;
      std::cout << lineBuffer << std::endl;
      throw std::logic_error("not 5 things to scan!");
    }
    emissionScoresVector.push_back(a_value);
    emissionScoresVector.push_back(c_value);
    emissionScoresVector.push_back(g_value);
    emissionScoresVector.push_back(t_value);


    std::getline(fileInputStream, lineBuffer);
  } while (lineBuffer[0] != 'I');

  return emissionScoresVector;
}



uint16_t attemptTraceback(uint32_t windowStart, uint32_t windowEnd, uint32_t phmmPosition, vector<int8_t> emissionScores, FastaVector* fastaVector) {

  uint32_t phmmLength = emissionScores.size() / 4;
  int16_t currentAccumulatedScore = 0;
  int16_t highestAccumulatedScore = 0;
  uint32_t windowLength = windowEnd - windowStart;
  uint32_t phmmStartPosition = phmmPosition - windowLength;

  for (uint32_t diagonalStep = 0; diagonalStep <= windowLength; diagonalStep++) {
    uint32_t currentSequencePosition = windowStart + diagonalStep;
    uint32_t currentPhmmPosition = phmmStartPosition + diagonalStep;
    char symbol = fastaVector->sequence.charData[currentSequencePosition];
    uint8_t encodedSymbol = encodeSymbol(symbol);
    int8_t* phmmVector = &emissionScores.data()[currentPhmmPosition * 4];
    int8_t matchScore = phmmVector[encodedSymbol];
    currentAccumulatedScore += matchScore;

    if (currentAccumulatedScore > highestAccumulatedScore) {
      highestAccumulatedScore = currentAccumulatedScore;
    }
    if (currentAccumulatedScore < 0) {
      currentAccumulatedScore = 0;
    }
    if (currentAccumulatedScore >= 256) {
      return currentAccumulatedScore;
    }
  }

  std::cout << "traceback required!" << std::endl;
  //now, if we still haven't reached the threshold, walkback from the start of the window
  currentAccumulatedScore = highestAccumulatedScore;
  uint32_t maxTraceback = std::min(windowStart, phmmStartPosition);

  for (uint32_t tracebackStep = 1; tracebackStep <= maxTraceback; tracebackStep++) {
    uint32_t currentSeqPosition = windowStart - tracebackStep;
    uint32_t currentPhmmPosition = phmmStartPosition - tracebackStep;
    char symbol = fastaVector->sequence.charData[currentSeqPosition];
    uint8_t encodedSymbol = encodeSymbol(symbol);
    int8_t* phmmVector = &emissionScores.data()[currentPhmmPosition * 4];
    int8_t matchScore = phmmVector[encodedSymbol];

    currentAccumulatedScore += matchScore;

    if (currentAccumulatedScore > highestAccumulatedScore) {
      highestAccumulatedScore = currentAccumulatedScore;
    }
    if (currentAccumulatedScore < 0) {
      break;
    }
    if (currentAccumulatedScore >= 256) {
      return currentAccumulatedScore;
    }
  }
  std::cout << "traceforward required!" << std::endl;
  //if we STILL haven't reached the threshold, check with a walkforward
  currentAccumulatedScore = highestAccumulatedScore;
  uint32_t seqRemaining = fastaVector->sequence.count - windowEnd;
  uint32_t phmmRemaining = phmmLength - phmmPosition;
  uint32_t maxTraceforward = std::min(seqRemaining, phmmRemaining);
  for (uint32_t traceforwardStep = 1; traceforwardStep < maxTraceforward; traceforwardStep++) {
    uint32_t currentSeqPosition = windowEnd + traceforwardStep;
    uint32_t currentPhmmPosition = phmmPosition + traceforwardStep;

    char symbol = fastaVector->sequence.charData[currentSeqPosition];
    uint8_t encodedSymbol = encodeSymbol(symbol);
    int8_t* phmmVector = &emissionScores.data()[currentPhmmPosition * 4];
    int8_t matchScore = phmmVector[encodedSymbol];
    currentAccumulatedScore += matchScore;


    if (currentAccumulatedScore > highestAccumulatedScore) {
      highestAccumulatedScore = currentAccumulatedScore;
    }
    if (currentAccumulatedScore < 0) {
      break;
    }
    if (currentAccumulatedScore >= 256) {
      return currentAccumulatedScore;
    }
  }

  return highestAccumulatedScore;
}


uint16_t attemptTracebackFloat(uint32_t windowStart, uint32_t windowEnd, uint32_t phmmPosition, vector<float> emissionScores, FastaVector* fastaVector) {
  uint32_t phmmLength = emissionScores.size() / 4;
  float currentAccumulatedScore = 0;
  float highestAccumulatedScore = 0;
  uint32_t windowLength = windowEnd - windowStart;
  uint32_t phmmStartPosition = phmmPosition - windowLength;

  for (uint32_t diagonalStep = 0; diagonalStep < windowLength; diagonalStep++) {
    uint32_t currentSequencePosition = windowStart + diagonalStep;
    uint32_t currentPhmmPosition = phmmStartPosition + diagonalStep;
    char symbol = fastaVector->sequence.charData[currentSequencePosition];
    uint8_t encodedSymbol = encodeSymbol(symbol);
    float* phmmVector = &emissionScores.data()[currentPhmmPosition * 4];
    float matchScore = phmmVector[encodedSymbol];
    currentAccumulatedScore += matchScore;

    if (currentAccumulatedScore > highestAccumulatedScore) {
      highestAccumulatedScore = currentAccumulatedScore;
    }
    if (currentAccumulatedScore < 0) {
      currentAccumulatedScore = 0;
    }
    if (currentAccumulatedScore >= 256) {
      return currentAccumulatedScore;
    }
  }

  //now, if we still haven't reached the threshold, walkback from the start of the window
  currentAccumulatedScore = highestAccumulatedScore;
  uint32_t maxTraceback = std::min(windowStart, phmmStartPosition);

  for (uint32_t tracebackStep = 1; tracebackStep <= maxTraceback; tracebackStep++) {
    uint32_t currentSeqPosition = windowStart - tracebackStep;
    uint32_t currentPhmmPosition = phmmStartPosition - tracebackStep;
    char symbol = fastaVector->sequence.charData[currentSeqPosition];
    uint8_t encodedSymbol = encodeSymbol(symbol);
    float* phmmVector = &emissionScores.data()[currentPhmmPosition * 4];
    float matchScore = phmmVector[encodedSymbol];

    currentAccumulatedScore += matchScore;

    if (currentAccumulatedScore > highestAccumulatedScore) {
      highestAccumulatedScore = currentAccumulatedScore;
    }
    if (currentAccumulatedScore < 0) {
      break;
    }
    if (currentAccumulatedScore >= 256) {
      return currentAccumulatedScore;
    }
  }

  //if we STILL haven't reached the threshold, check with a walkforward
  currentAccumulatedScore = highestAccumulatedScore;
  uint32_t seqRemaining = fastaVector->sequence.count - windowEnd;
  uint32_t phmmRemaining = phmmLength - phmmPosition;
  uint32_t maxTraceforward = std::min(seqRemaining, phmmRemaining);
  for (uint32_t traceforwardStep = 0; traceforwardStep < maxTraceforward; traceforwardStep++) {
    uint32_t currentSeqPosition = windowEnd + traceforwardStep;
    uint32_t currentPhmmPosition = phmmPosition + traceforwardStep;

    char symbol = fastaVector->sequence.charData[currentSeqPosition];
    uint8_t encodedSymbol = encodeSymbol(symbol);
    float* phmmVector = &emissionScores.data()[currentPhmmPosition * 4];
    float matchScore = phmmVector[encodedSymbol];
    currentAccumulatedScore += matchScore;

    if (currentAccumulatedScore > highestAccumulatedScore) {
      highestAccumulatedScore = std::max(highestAccumulatedScore, currentAccumulatedScore);
      highestAccumulatedScore = currentAccumulatedScore;
    }
    if (currentAccumulatedScore < 0) {
      break;
    }
    if (currentAccumulatedScore >= 256) {
      return currentAccumulatedScore;
    }
  }

  return highestAccumulatedScore;
}


bool symbolIsAmbiguous(char symbol) {
  switch (symbol) {
  case 'a': case 'A': case 'c': case 'C': case 'g': case 'G': case 't': case 'T': return false;
  default: return true;
  }
  return true;
}

vector<int8_t> deprojectAndRoundEmissions(vector<float>& emissionScores, const uint8_t thresholdScore) {
  vector<int8_t> roundedEmissions(emissionScores.size());

  for (uint32_t i = 0; i < emissionScores.size();i++) {
    float scoreAdjustedForThreshold = (emissionScores[i] / 256) * thresholdScore;
    roundedEmissions[i] = std::round(scoreAdjustedForThreshold);
  }

  return roundedEmissions;
}
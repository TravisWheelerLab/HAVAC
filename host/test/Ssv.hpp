#pragma once

#include <memory>
#include <vector>
#include <sstream>
extern "C" {
  #include <FastaVector.h>
  #include <p7HmmReader.h>
}



using std::shared_ptr;
using std::vector;
using std::make_shared;

class ReferenceSsvHit {
  public:
  uint32_t sequenceNumber;
  uint32_t phmmNumber;
  uint32_t sequencePosition;
  uint32_t phmmPosition;

  std::string toString(){
    std::stringstream stringStream;
    stringStream << "sequence $" << this->sequenceNumber << ", position " << this->sequencePosition <<
      "; phmm #" << this->phmmNumber << " position " << this->phmmPosition;
    return stringStream.str();
  }
};

shared_ptr<vector<ReferenceSsvHit>> HitsFromSsv(FastaVector* fastaVector, 
  P7HmmList* phmmList, const float desiredPValue);


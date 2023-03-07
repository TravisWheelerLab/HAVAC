#include "VerifiedHit.hpp"
#include <sstream>

VerifiedHit::VerifiedHit()
  :sequencePosition(0), sequenceIndex(0), phmmPosition(0), phmmIndex(0) {
}

VerifiedHit::VerifiedHit(uint64_t sequencePosition, uint32_t sequenceIndex,
  uint32_t phmmPosition, uint32_t phmmIndex)
  :sequencePosition(sequencePosition), sequenceIndex(sequenceIndex),
  phmmPosition(phmmPosition), phmmIndex(phmmIndex) {
}

std::string VerifiedHit::toString() {
  std::stringstream stringStream;
  stringStream << "sequence $" << this->sequenceIndex << ", position "<< this->sequencePosition <<
    "; phmm #"<< this->phmmIndex << " position "<< this->phmmPosition;
  return stringStream.str();
}
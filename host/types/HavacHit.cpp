#include "HavacHit.hpp"
#include <sstream>
#include <stdexcept>

HavacHit::HavacHit(const uint64_t sequencePosition, const uint32_t sequenceIndex, const uint32_t phmmPosition, const uint32_t phmmIndex)
:sequencePosition(sequencePosition),
sequenceIndex(sequenceIndex),
phmmPosition(phmmPosition),
phmmIndex(phmmIndex){

}

std::string HavacHit::toString() {
  std::stringstream stringStream;
  stringStream << "sequence $" << this->sequenceIndex << ", position " << this->sequencePosition <<
    "; phmm #" << this->phmmIndex << " position " << this->phmmPosition;
  return stringStream.str();
}
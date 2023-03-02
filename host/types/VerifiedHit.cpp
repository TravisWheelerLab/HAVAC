#include "VerifiedHit.hpp"
#include <format>

VerifiedHit::VerifiedHit()
  :sequencePosition(0), sequenceIndex(0), phmmPosition(0), phmmIndex(0) {
}

VerifiedHit::VerifiedHit(uint64_t sequencePosition, uint32_t sequenceIndex,
  uint32_t phmmPosition, uint32_t phmmIndex)
  :sequencePosition(sequencePosition), sequenceIndex(sequenceIndex),
  phmmPosition(phmmPosition), phmmIndex(phmmIndex) {
}

std::string VerifiedHit::toString() {
  return std::format("sequence #{}, position {}; phmm #{}, position {}", this->sequenceIndex,
    this->sequencePosition, this->phmmIndex, this->phmmPosition);
}
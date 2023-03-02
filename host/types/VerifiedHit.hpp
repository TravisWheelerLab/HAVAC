#ifndef HAVAC_VERIFIED_HIT_HPP
#define HAVAC_VERIFIED_HIT_HPP

#include <cstdint>
#include <string>

class VerifiedHit {
public:
  VerifiedHit();
  VerifiedHit(uint64_t sequencePosition, uint32_t sequenceIndex, uint32_t phmmPosition, uint32_t phmmIndex);
  uint64_t sequencePosition;
  uint32_t sequenceIndex;
  uint32_t phmmPosition;
  uint32_t phmmIndex;

  std::string toString();
};

#endif
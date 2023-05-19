#ifndef HAVAC_HARDWARE_HIT_REPORT_HPP
#define HAVAC_HARDWARE_HIT_REPORT_HPP

#include <cstdint>
#include <string>
#include <device/PublicDefines.h>

#define GROUP_HIT_BITS_NUM_BYTES (NUM_CELL_GROUPS /8)

class HardwareHitReport {
public:
  HardwareHitReport();
  HardwareHitReport(uint32_t sequencePassIndex, uint8_t groupHitBits[GROUP_HIT_BITS_NUM_BYTES],
    uint32_t phmmPosition);

  static_assert(NUM_CELL_GROUPS % 8 == 0, "num groups is required to be a multiple of 8, but was not.");

  uint32_t phmmPosition;
  uint32_t sequencePassIndex;
  uint8_t groupHitBits[GROUP_HIT_BITS_NUM_BYTES];
  std::string toString();
};

#endif
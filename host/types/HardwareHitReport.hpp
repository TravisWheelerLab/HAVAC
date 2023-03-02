#ifndef HAVAC_HARDWARE_HIT_REPORT_HPP
#define HAVAC_HARDWARE_HIT_REPORT_HPP

#include <cstdint>

class HardwareHitReport{
  public:
    HardwareHitReport();
    HardwareHitReport(uint32_t sequencePassIndex, uint32_t sequenceGroupIndex,
      uint32_t phmmPosition);
    uint32_t sequencePassIndex;
    uint32_t sequenceGroupIndex;
    uint32_t phmmPosition;
};

#endif
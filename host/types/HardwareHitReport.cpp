#include "HardwareHitReport.hpp"


HardwareHitReport::HardwareHitReport()
  :sequencePassIndex(0),
  sequenceGroupIndex(0),
  phmmPosition(0) {
}


HardwareHitReport::HardwareHitReport(uint32_t sequencePassIndex, uint32_t sequenceGroupIndex,
  uint32_t phmmPosition)
  :sequencePassIndex(sequencePassIndex),
  sequenceGroupIndex(sequenceGroupIndex),
  phmmPosition(phmmPosition) {

}
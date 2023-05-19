#include "HardwareHitReport.hpp"
#include <sstream>
#include <bitset>
#include <string.h>
#include "../../device/PublicDefines.h"

HardwareHitReport::HardwareHitReport()
  :sequencePassIndex(0),
  groupHitBits({ 0 }),
  phmmPosition(0) {
}


HardwareHitReport::HardwareHitReport(uint32_t sequencePassIndex, uint8_t groupHitBits[GROUP_HIT_BITS_NUM_BYTES],
  uint32_t phmmPosition)
  :sequencePassIndex(sequencePassIndex),
  phmmPosition(phmmPosition) {
  memcpy(this->groupHitBits, groupHitBits, GROUP_HIT_BITS_NUM_BYTES);

}


std::string HardwareHitReport::toString() {
  std::stringstream ss;
  ss << "phmm position: " << this->phmmPosition << ", sequence segment #" << this->sequencePassIndex << ", groups: ";
  for (uint32_t i = GROUP_HIT_BITS_NUM_BYTES; i > 0;i--) {
    ss << std::bitset<8>(this->groupHitBits[i]).to_string() << "_";
  }
  ss << std::bitset<8>(this->groupHitBits[0]).to_string();

  return ss.str();
}
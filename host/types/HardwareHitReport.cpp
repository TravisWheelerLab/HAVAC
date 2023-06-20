#include "HardwareHitReport.hpp"
#include <immintrin.h>
#include <x86intrin.h>
#include <sstream>
#include <bitset>
#include <string.h>

HardwareHitReport::HardwareHitReport()
  :phmmPosition(),
  sequencePassIndex(),
  groupHitBits() {
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
  bool foundGroupReportingHit;
  for (uint32_t i = 0; i < GROUP_HIT_BITS_NUM_BYTES;i++) {
    for (uint32_t bitInByte = 0; bitInByte < 8; bitInByte++) {
      if (this->groupHitBits[i] & (1 << bitInByte)) {
        if (foundGroupReportingHit) {
          ss << ", ";
        }
        else {
          ss << "[";
        }
        foundGroupReportingHit = true;
        ss << i * 8 + bitInByte;
      }
    }
  }
  if (foundGroupReportingHit) {
    ss << "]";
  }
  else {
    ss << "NoneFound";
  }

  return ss.str();
}


std::string HardwareHitReport::toGroupListString() {
  std::stringstream ss;

  ss << "phmm position: " << this->phmmPosition << ", sequence segment #" <<
    this->sequencePassIndex << ", groups: [";
  for (size_t i = 0; i < GROUP_HIT_BITS_NUM_BYTES;i++) {
    uint8_t thisGroupBits = this->groupHitBits[i];
    while (thisGroupBits) {

      uint32_t groupIndex = _bit_scan_forward(thisGroupBits);
      //remove the group we're about to validaten from the sequenceGroupIndex data
      thisGroupBits ^= 1 << groupIndex;

      groupIndex += (i * 8);

      ss << groupIndex << " ";
    }
  }
  ss << "]";
  return ss.str();
}
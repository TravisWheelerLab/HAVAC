#ifndef HAVAC_HOST_MAIN_HPP
#define HAVAC_HOST_MAIN_HPP

#include "xrt/xrt_bo.h"
#include <boost/optional.hpp>

struct HavacHostBuffers{
  boost::optional<xrt::bo> sequenceBuffer;
  boost::optional<xrt::bo> phmmBuffer;
  boost::optional<xrt::bo> hitReportBuffer;
};

struct HavacHostHitReport{
  uint32_t phmmIndex;
  uint32_t sequenceIndex;
  uint64_t groupsPassingThreshold[2];
};


#endif

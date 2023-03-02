#ifndef HAVAC_PHMM_REPROJECTION_H
#define HAVAC_PHMM_REPROJECTION_H

extern "C" {
  #include "p7HmmReader.h"
}
void p7HmmProjectForThreshold256(struct P7Hmm* phmm, float desiredPValue, int8_t* outputArray);

#endif
#ifndef HAVAC_PHMM_REPROJECTION_H
#define HAVAC_PHMM_REPROJECTION_H

extern "C" {
  #include "../P7HmmReader/src/p7HmmReader.h"
}
void p7HmmProjectForThreshold256(struct P7Hmm* phmm, float desiredPValue, int8_t* projectedScoresArray);

#endif
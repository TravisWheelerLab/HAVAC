#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "PhmmReprojection.h"

#define eslSMALLX1    5e-9
#define eslCONST_LOG2 0.69314718055994529
#define HAVAC_LOG_ONE_FOURTH -1.3862943611198906

//this function is straight out of HMMER
double esl_gumbel_invsurv(double p, double mu, double lambda) {
  /* The real calculation is mu - ( log(-1. * log(1-p)) / lambda).
  *  But there's a problem with small p:
  *     for p<1e-15, 1-p will be viewed as 1, so
  *     log ( -log(1-p) ) == log (0) -> inf
  *  Instead, use two approximations;
  *    (1) log( 1-p) ~= -p   for small p (e.g. p<0.001)
  *      so log(-1. * log(1-p)) ~= log(p)
  *    (2) log (p) ~= (p^p - 1) / p
  *
  *    See notes Mar 1, 2010.
  */
  double log_part;
  if (p < eslSMALLX1) {
    log_part = (pow(p, p) - 1) / p;
  }
  else {
    log_part = log(-1. * log(1 - p));
  }

  //test 2
  return mu - (log_part / lambda);
}


//mu and lambda come from the STATS MSV line on the profile hmm
//max length also from profile Hmm, same with modelLength
float findThreshold256ScalingFactor(double p, double mu, double lambda, double maxLength, double modelLength) {
  double gumbleInverse = esl_gumbel_invsurv(p, mu, lambda);
  float tloop = logf((float)maxLength / (float)(maxLength + 3));
  float tloop_total = tloop * maxLength;
  float tmove = logf(3.0f / (float)(maxLength + 3));
  float tbmk = logf(2.0f / ((float)modelLength * (float)(modelLength + 1)));
  float transitionEToC = logf(1.0f / 2);
  float bg_P = (float)maxLength / (float)(maxLength + 1);
  float nullSc = (float)maxLength * log(bg_P) + log(1.0 - bg_P);

  float scThresh =
    (gumbleInverse * eslCONST_LOG2) + nullSc - (tmove + tloop_total + tmove + tbmk + transitionEToC);
  // that's in nats; let's shift to bits
  scThresh /= eslCONST_LOG2;

  float scaleFactor = 255 / scThresh;

  return scaleFactor;
}

void p7HmmProjectForThreshold256(struct P7Hmm* phmm, float desiredPValue) {
  //there's some log/exp tricks here, the naive approach to rescaling these values is to use the formula:
  // scaledScore = eslCONST_LOG2 *log(exp(-1* baseScore)/backgroundProbability) * scaleFactor;
  // where backgroundProbability = .25, therefore
  // scaledScore = eslCONST_LOG2 *log(exp(-1* baseScore)/.25) * scaleFactor;
  //
  //or more mathy, y = c1 * log(exp(-1*x)/.25)* c2
  //  = c1 * log(exp(x)^-1/.25 * c2 = c1 * log(1/exp(x)/.25) * c2
  //  = c1 * log(1/(exp(x)*.25)) *c2 =  c1 * (log(1)-log(exp(x)*.25))) *c2
  //  = c1 * (-log(exp(x)*.25))) *c2 = c1 * (log(exp(x)*.25))) *c2 * -1
  //  = c1 * (log(exp(x))+log(.25) ) *c2 * -1
  //  = c1 * (x+log(.25) ) *c2 * -1
  //  = (x + log(.25)) * (c1 * c2 * -1);
  //assuming c1, c2, and x are real (they're floats).
  //therefore, we can scale the scores simply with c1 * c2 * -1 * (x + log(.25))
  // where c1 is the log(2), c2 is the scaleFactor for the desired p-value, and log(.25) can be a precomputed constant.
  float mu = phmm->stats.msvGumbelMu;
  float lambda = phmm->stats.msvGumbleLambda;
  float modelLength = phmm->header.modelLength;
  float maxLength = phmm->header.maxLength;
  float scaleFactor = findThreshold256ScalingFactor(desiredPValue, mu, lambda, maxLength, modelLength);

  const float scoreMultiplier = eslCONST_LOG2 * scaleFactor * -1.0f;

  uint32_t alphabetCardinality = p7HmmGetAlphabetCardinality(phmm);

  for (size_t i = 0; i < alphabetCardinality * modelLength; i++) {
    //phmm->model.matchEmissionScores[i] = eslCONST_LOG2 * log(exp(-1* phmm->model.matchEmissionScores[i])/backgroundProbability) * scaleFactor;
    phmm->model.matchEmissionScores[i] = scoreMultiplier * (HAVAC_LOG_ONE_FOURTH + phmm->model.matchEmissionScores[i]);
  }

}

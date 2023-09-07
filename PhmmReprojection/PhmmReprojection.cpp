#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "PhmmReprojection.h"

#define HAVAC_GUMBEL_EPSILON 5e-9
#define HAVAC_NAT_LOG_2 0.69314718055994529
#define HAVAC_NAT_LOG_ONE_FOURTH -1.3862943611198906

//this function is adapted from the easel function esl_gumbel_invsurv (esl_gumbel.h)
//defines the inverse gumbel distribution, that is,
// for a given p-value, what score is required to hit that value?
// note that this is for the full model, not the single-hit reduced model we're working with,
// so we'll need to add some penalties later to make it come out to the right values.
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
  double log_part = (p < HAVAC_GUMBEL_EPSILON) ?
    (pow(p, p) - 1) / p :
    log(-1. * log(1 - p));
  return mu - (log_part / lambda);
}


//mu and lambda come from the STATS MSV line on the profile hmm
//max length also from profile Hmm, same with modelLength
float findThreshold256ScalingFactor(const struct P7Hmm* const phmm, const float pValue) {
  const float mu = phmm->stats.msvGumbelMu;
  const float lambda = phmm->stats.msvGumbelLambda;
  const float maxLength = phmm->header.maxLength;
  const float modelLength = phmm->header.modelLength;
    //given the probability of survival, give me the score. what score gives the p value?
    //what score do we need to hit the required p value? this value is for the full model. 
    double scoreRequiredForFullModelPvalue = esl_gumbel_invsurv(pValue, mu, lambda); //inverse survival (gumbel distribution)
    //we're only working with the single-hit model, so we need to add some penalties to make the probabilities match.
    float nStateLoopPenalty = logf(maxLength / (maxLength + 3));  //probablility of staying in the n state (or c state), looping.
    float nStateLoopPenaltyTotal = nStateLoopPenalty * maxLength;  //total length penalty for the max sequence length
    float nStateEscapePenalty = logf(3.0f / (maxLength + 3));  //penalty for escaping the n or c states (increases by 1 every time seq len doubles)
    float bStateToAnyMStatePenalty = logf(2.0f / (modelLength * (modelLength + 1)));  //penalty for transition from b to kth 'm' option (msubk)
    float transitionEToC = logf(1.0f / 2);
    float coreModelAdjustment = (nStateEscapePenalty + nStateLoopPenaltyTotal +
      nStateEscapePenalty + bStateToAnyMStatePenalty + transitionEToC);

    float backgroundLoopProbability = maxLength / (maxLength + 1); //background loop probability
    float backgroundLoopPenaltyTotal = maxLength * log(backgroundLoopProbability);  //total length penalty for the max sequence length background, change to bg_loop total or something
    float backgroundMovePenalty = log(1.0 - backgroundLoopProbability);
    float backgroundScore = backgroundLoopPenaltyTotal + backgroundMovePenalty;

    float thresholdScoreInNats =
      (scoreRequiredForFullModelPvalue * HAVAC_NAT_LOG_2) + backgroundScore - coreModelAdjustment;
    // that's in nats; let's shift to bits
    float thresholdScoreInBits = thresholdScoreInNats / HAVAC_NAT_LOG_2;

    float scaleFactor = 256.0f / thresholdScoreInBits; //a hit should be triggered if a cell hits 256 (causes an 8-bit int overflow)

    return scaleFactor;


  //original implementation, left for posterity
  // #define eslCONST_LOG2 0.69314718055994529
  // // (these are all in nats)
  // float tloop = logf((float)maxLength / (float)(maxLength + 3));
  // float tloop_total = tloop * maxLength;
  // float tmove = logf(3.0f / (float)(maxLength + 3));
  // float tbmk = logf(2.0f / ((float)modelLength * (float)(modelLength + 1)));
  // float tec = logf(1.0f / 2);
  // float bg_P = (float)maxLength / (float)(maxLength + 1);
  // float null_sc = (float)maxLength * log(bg_P) + log(1.0 - bg_P);

  // float invP = esl_gumbel_invsurv(pValue, mu, lambda);
  // // (5) Get the score threshold
  //   float sc_thresh = (invP * eslCONST_LOG2) + null_sc - (tmove + tloop_total + tmove + tbmk + tec);
  // // that's in nats; let's shift to bits
  // sc_thresh /= eslCONST_LOG2;

  // return 256/sc_thresh; 
}

inline float emissionScoreToProjectedScore(const float emissionScore, const float scoreMultiplier) {
  //this code takes a value from negative natural log likelihood space, and converts it
  //into a bits score scaled such that the threshold score will be 256.
  //this is a simplification of the following function: 
  //log2((exp(-emissionScore) / (1/4))) * multiplier;
  // where the 1/4 is the background distribution (even across 4 nucleotides)
  //let's rewrite this so we can simplify: 
  //y = log2(e^-s / (1/4) ) * m
  //y = log2(4e^-s ) * m
  //y = (log2(4) + log2(e^-s) ) * m
  //y = (log2(4) - s*log2(e)) * m
  //y = (2 - s*log2(e)) * m
  //y = -log2(e) * (s -(2/log2(e)) ) * m

  // return log2((exp(-emissionScore)/.25))*multiplier;  //this is equivalent to the following
  const float LOG2_E = 1.44269504089;
  float projectedScore = -LOG2_E * (emissionScore - (2 / LOG2_E)) * scoreMultiplier;
  return projectedScore;
}

void p7HmmProjectForThreshold256(const struct P7Hmm* const phmm, const float desiredPValue, int8_t* outputArray) {

  float modelLength = phmm->header.modelLength;
  float scoreMultiplier = findThreshold256ScalingFactor(phmm, desiredPValue);
  uint32_t alphabetCardinality = p7HmmGetAlphabetCardinality(phmm);

  for (size_t i = 0; i < alphabetCardinality * modelLength; i++) {
    float valueFromHmmFile = phmm->model.matchEmissionScores[i];
    float projectedScore = emissionScoreToProjectedScore(valueFromHmmFile, scoreMultiplier);

    outputArray[i] = round(projectedScore);
  }
}

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "PhmmReprojection.h"

#define eslSMALLX1    5e-9
#define eslCONST_LOG2 0.69314718055994529
#define HAVAC_LOG_ONE_FOURTH -1.3862943611198906

//this function is straight out of HMMER
//defines the inverse gumble distribution, that is,
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
  //given the probability of survival, give me the score. what score gives the p value?
  //what score do we need to hit the required p value? this value is for the full model. 
  double scoreRequiredForFullModelPvalue = esl_gumbel_invsurv(p, mu, lambda); //inverse survival (gumble distribution)
  //we're only working with the single-hit model, so we need to add some penalties to make the probabilities match.
  float nStateLoopPenalty = logf((float)maxLength / (float)(maxLength + 3));  //probablility of staying in the n state (or c state), looping.
  float nStateLoopPenaltyTotal = nStateLoopPenalty * maxLength;  //total length penalty for the max sequence length
  float nStateEscapePenalty = logf(3.0f / (float)(maxLength + 3));  //penalty for escaping the n or c states (increases by 1 every time seq len doubles)
  float bStateToKthMStatePenalty = logf(2.0f / ((float)modelLength * (float)(modelLength + 1)));  //penalty for transition from b to kth 'm' option (msubk)
  float transitionEToC = logf(1.0f / 2);
  float coreModelAdjustment = (nStateEscapePenalty + nStateLoopPenaltyTotal +
    nStateEscapePenalty + bStateToKthMStatePenalty + transitionEToC);

  float backgroundLoopProbability = (float)maxLength / (float)(maxLength + 1); //background loop probability
  float backgroundLoopPenaltyTotal = (float)maxLength * log(backgroundLoopProbability);  //total length penalty for the max sequence length background, change to bg_loop total or something
  float backgroundMovePenalty = log(1.0 - backgroundLoopProbability);
  float backgroundScore = backgroundLoopPenaltyTotal + backgroundMovePenalty;

  float baseScoreThreshold =
    (scoreRequiredForFullModelPvalue * eslCONST_LOG2) - coreModelAdjustment + backgroundScore;
  // that's in nats; let's shift to bits
  baseScoreThreshold /= eslCONST_LOG2;

  float scaleFactor = 256.0f / baseScoreThreshold; //a hit should be triggered if a cell hits 256 (causes an 8-bit int overflow)

  return scaleFactor;
}


float generateScoreMultiplierForPhmmScore(const struct P7Hmm* const phmm, const float desiredPValue) {
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

  return scoreMultiplier;
}

 float projectPhmmScoreWithMultiplier(const float phmmScore, const float scoreMultiplier) {
  return scoreMultiplier * (HAVAC_LOG_ONE_FOURTH + phmmScore);
}

void p7HmmProjectForThreshold256(const struct P7Hmm* const phmm, const float desiredPValue, int8_t* outputArray) {

  float modelLength = phmm->header.modelLength;
  float scoreMultiplier = generateScoreMultiplierForPhmmScore(phmm, desiredPValue);
  uint32_t alphabetCardinality = p7HmmGetAlphabetCardinality(phmm);

  for (size_t i = 0; i < alphabetCardinality * modelLength; i++) {
    //phmm->model.matchEmissionScores[i] = eslCONST_LOG2 * log(exp(-1* phmm->model.matchEmissionScores[i])/backgroundProbability) * scaleFactor;
    outputArray[i] = scoreMultiplier * (HAVAC_LOG_ONE_FOURTH + phmm->model.matchEmissionScores[i]);
  }
}

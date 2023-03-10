#ifndef HAVAC_PHMM_REPROJECTION_H
#define HAVAC_PHMM_REPROJECTION_H

extern "C" {
  #include "p7HmmReader.h"
}

/// @brief Finds the multiplier that would need to be applied to each score in the phmm to make it
///   so that a score of exactly 256 would count as a threshold hit for the given p-value.
///   Please note that the .phmm file stores the scores in negative log-odds, so this multiplier
///   will multiply by -1, thus creating positive scores for good matches, and negative scores for bad matches.
/// @param phmm pointer to the phmm to project. these score multipliers are on a per-model basis
/// @param desiredPValue The p value required to generate a score of 256, thus generating a hit.
/// @return  score multipler needed for each score in the given phmm.
float generateScoreMultiplierForPhmmScore(const struct P7Hmm* const phmm, const float desiredPValue);


/// @brief projects a score as given in a phmm file into a score for use with a constant 256 hit threshold.
/// @param phmmScore score from the phmm file. this will be a negative log-odds score.
/// @param scoreMultiplier multiplier found via the generateScoreMultiplierForPhmmScore() function.
/// @return projected score in int8_t space.
float projectPhmmScoreWithMultiplier(const float phmmScore, const float scoreMultiplier);

/// @brief projects the scores inside the given phmm such that an accumulated score of 256 represents a hit at
///     the given p-value. the phmm is not modified in this function. Instead, the projected score are used to
///     populate the outputArray 
/// @param phmm model to project the scores into int8_t space.
/// @param desiredPValue The p value required to generate a score of 256, thus generating a hit.
/// @param outputArray flattened array that will contain the projected scores from the phmm.
void p7HmmProjectForThreshold256(const struct P7Hmm* const phmm, const float desiredPValue, int8_t* outputArray);

#endif
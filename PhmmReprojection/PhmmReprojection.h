#ifndef HAVAC_PHMM_REPROJECTION_H
#define HAVAC_PHMM_REPROJECTION_H

extern "C" {
  #include "p7HmmReader.h"
}

/// @brief Finds the score multiplier needed to project the phmm scores so that
/// a hit will be generated on an accumulated score of 256.
/// @param phmm phmm to project. since the scaling factor depends on the maxLength attribute,
///   the scaling factor will be unique to a particular hmm.
/// @param pValue p-value of a diagonal that should generate a hit.
float findThreshold256ScalingFactor(const struct P7Hmm* const phmm, const float pValue);


/// @brief projects a given score taken from the phmm file (negative natural log likelihood) into
///  a score havac can use (in bits) with the score scaling to make hits occur at a score of 256.
/// @param emissionScore score from the phmm file. this value is stored as a negative log likelihood value.
/// @param scoreMultiplier multiplier to apply to the score (in bits) to make hits occur at a score of 256.
/// @return 
float emissionScoreToProjectedScore(const float emissionScore, const float scoreMultiplier);


/// @brief projects the scores inside the given phmm such that an accumulated score of 256 represents a hit at
///     the given p-value. the phmm is not modified in this function. Instead, the projected score are used to
///     populate the outputArray 
/// @param phmm model to project the scores into int8_t space.
/// @param desiredPValue The p value required to generate a score of 256, thus generating a hit.
/// @param outputArray flattened array that will contain the projected scores from the phmm.
void p7HmmProjectForThreshold256(const struct P7Hmm* const phmm, const float desiredPValue, int8_t* outputArray);

#endif
#ifndef HAVAC_HIT_VERIFIER_HPP
#define HAVAC_HIT_VERIFIER_HPP


#include "types/HardwareHitReport.hpp"
#include "types/VerifiedHit.hpp"
#include <boost/optional.hpp>
#include "HavacSSV.hpp"
#include <tuple>
#include <memory>

extern "C" {
  #include <FastaVector.h>
  #include <p7HmmReader.h>
}

using std::shared_ptr;
using std::vector;

#define SSV_VERIFICATION_PRIOR_FLANK_RANGE 64
#define SSV_VERIFICATION_POST_FLANK_RANGE 16
#define SSV_PHMM_VERIFICATION_RANGE (SSV_VERIFICATION_PRIOR_FLANK_RANGE + SSV_VERIFICATION_POST_FLANK_RANGE)

struct SsvReferenceBounds {
  uint64_t sequenceStartPosition;
  uint64_t sequenceEndPosition;
  uint32_t phmmStartPosition;
  uint32_t phmmEndPosition;
  uint32_t phmmIndexInList;
};

/// @brief Class to verify any hits generated by HAVAC hardware, and to identify the exact locations of the hits described.
///   Objects of this class should only be used once, and a new HitVerifier should be generated if another HAVAC run needs to be verified.
class HitVerifier {
public:
  HitVerifier(FastaVector *fastaVector, P7HmmList *phmmList);
  HitVerifier(HitVerifier& hv) = delete;
  HitVerifier(HitVerifier&& hv) = delete;
  ~HitVerifier() = default;

  /// @brief verifies the hits reported from the hardware. This throws out hits that only
  ///     happened because they strided a phmm or sequence boundary, and resolves
  ///     the position of any hits in a hit group to an exact cell location.
  ///     Note: these hits are not guaranteed to be unique, it's unlikely but possible that
  ///     hits may show up multiple times if located near each other. 
  /// @param hits  shared_ptr to the list of hits returned by the hardware
  /// @return  shared_ptr to a list of hits that could be verified via a software SSV implementation.
  /// @param desiredPValue p-value required to verify a hit. this should be the same value given to HAVAC at hardware runtime.
  shared_ptr<vector<VerifiedHit>> verify(shared_ptr<vector<HardwareHitReport>> hits, const float desiredPValue);


protected:
  FastaVector *fastaVector;
  P7HmmList *phmmList;
  shared_ptr<vector<float>> ssvCellScores;

  /// @brief attempts to verify a single hardware hit report. Any hits found via the reference SSV
  ///     are added to the verifiedHitList
  /// @param hardwareHitReport hit report generated by the hardware.
  /// @param verifiedHitList list of verified hits. any hits verified in this function are appended to this list
  /// @param desiredPValue p-value required to verify a hit. this should be the same value given to HAVAC at hardware runtime.
  void verifyHit(const HardwareHitReport& hardwareHitReport, shared_ptr<vector<VerifiedHit>> verifiedHitList, const float desiredPValue);
  
  /// @brief given a hardware hit report, find the index of the phmm the hit occurred in, and the position
  ///   in that phmm the hit corresponds to. This function will throw an exception if the position is
  ///   outside the possible range for the given phmmlist
  /// @param hardwareHitReport report to query for the phmm index and position
  /// @return tuple of the phmm index, and the hit position in said index.
  std::tuple<uint32_t, uint32_t> getPhmmIndexFromPosition(const HardwareHitReport& hardwareHitReport);

  /// @brief Invokes a reference SSV around the hit area to determine if the hit is genuine, and if so,
  ///   finds the exact cell where the hit occurrs. Instead of running an entire SSV matrix, this just explores
  ///   the diagonals that are on the cells that could've caused the hit.
  /// @param hitLocatedInPhmmNumber index of the phmm inside the phmmList where the hit occurred.
  /// @param sequenceNumber index of the sequence inside the fasta to explore with SSV.
  /// @param hitOccurredAtPhmmIndex the local position in the phmm where the hit occurred.
  /// @param localPhmmSsvLength length of the ssv matrix along the phmm axis.
  /// @param possibleSequenceIndexStart the first position in the sequence that could've caused the hit, relative to the sequence start
  /// @param possibleSequenceIndexEnd the last position in the sequence, exclusive, that could've caused the hit, relateive to the sequence start.
  ///   In other words, this is the location of the cell AFTER the final one that could be the hit!
  /// @param desiredPValue p-value that will generate a hit. this will likely be the same as what was used.
  ///     when projecting the PHMM for use in the HAVAC hardware.
  /// @param verifiedHitList  list of verified hits. if any hits are found in this function, they will be appended to
  ///   this vector.
  void verifyWithReferenceSsv(const uint32_t hitLocatedInPhmmNumber, const uint32_t sequenceNumber,
    const uint32_t hitOccurredAtPhmmIndex, const uint32_t possibleSequenceIndexStart, const uint32_t possibleSequenceIndexEnd,
    const float desiredPValue, shared_ptr<vector<VerifiedHit>> verifiedHitList);

  /// @brief for a given hit group, verifies any and all hits that might have occurred in the report.
  /// @param localPhmmPosition index in the phmm where the hit occurred, local to the given model
  /// @param phmmIndex the index of which model in the phmm file should be checked.
  /// @param sequenceSegmentIndex segment pass for the hit
  /// @param groupIndex index of the group that caused the hit. Since hit reports may have hits in 
  ///   multiple groups, this handles one group at a time. the group should be resolved from bit-position to index.
  /// @param verifiedHitList list to append any found hits to.
  /// @param desiredPValue p value needed to be reached for a hit to be reported.
  void verifyHitForGroup(uint32_t localPhmmPosition, uint32_t phmmIndex, uint32_t sequenceSegmentIndex, uint32_t groupIndex,
    shared_ptr<vector<VerifiedHit>> verifiedHitList, const float desiredPValue);

};
#endif
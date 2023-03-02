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
  HitVerifier(shared_ptr<FastaVector>& fastaVector, shared_ptr<P7HmmList>& phmmList);
  HitVerifier(HitVerifier& hv) = delete;
  HitVerifier(HitVerifier&& hv) = delete;
  ~HitVerifier() = default;
  shared_ptr<vector<VerifiedHit>> verify(shared_ptr<vector<HardwareHitReport>> hits);


protected:
  shared_ptr<FastaVector> fastaVector;
  shared_ptr<P7HmmList> phmmList;
  shared_ptr<vector<float>> ssvCellScores;


  void verifyHit(const HardwareHitReport& hardwareHitReport, shared_ptr<vector<VerifiedHit>> verifiedHitList);
  
  /// @brief given a hardware hit report, find the index of the phmm the hit occurred in, and the position
  ///   in that phmm the hit corresponds to. This function will throw an exception if the position is
  ///   outside the possible range for the given phmmlist
  /// @param hardwareHitReport report to query for the phmm index and position
  /// @return tuple of the phmm index, and the hit position in said index.
  std::tuple<uint32_t, uint32_t> getPhmmIndexFromPosition(const HardwareHitReport& hardwareHitReport);
  void verifyWithReferenceSsv(const uint32_t hitLocatedInPhmmNumber, const uint32_t sequenceNumber,
    const uint32_t localSsvPhmmStartPosition, const uint32_t localPhmmSsvLength,
    const uint32_t ssvSequenceRangeBegin, const uint32_t ssvSequenceRangeEnd);

};
#endif
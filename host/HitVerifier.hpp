#ifndef HAVAC_HIT_VERIFIER_HPP
#define HAVAC_HIT_VERIFIER_HPP


#include <boost/optional.hpp>
#include "HavacSSV.hpp"
#include "Havac.hpp"
#include <tuple>
#include <memory>

extern "C" {
  #include "FastaVector/src/FastaVector.h"
  #include "P7HmmReader/src/p7HmmReader.h"
}

using std::shared_ptr;
using std::unique_ptr;
#define SSV_VERIFICATION_PRIOR_FLANK_RANGE 64
#define SSV_VERIFICATION_POST_FLANK_RANGE 16
#define SSV_VERIFICATION_RANGE (SSV_VERIFICATION_PRIOR_FLANK_RANGE + SSV_VERIFICATION_POST_FLANK_RANGE)

struct SsvVerifiedHit {
  uint64_t sequencePosition;
  uint32_t sequenceIndex;
  uint32_t phmmPosition;
  uint32_t phmmIndex;
};

struct HavacHardwareHitReport {
  uint32_t sequencePassIndex;
  uint32_t sequenceGroupIndex;
  uint32_t phmmPosition;
};

struct SsvReferenceBounds {
  uint64_t sequenceStartPosition;
  uint64_t sequenceEndPosition;
  uint32_t phmmStartPosition;
  uint32_t phmmEndPosition;
  uint32_t phmmIndexInList;
};

class HitVerifier {
public:
  HitVerifier(shared_ptr<FastaVector>& fastaVector, shared_ptr<P7HmmList>& phmmList);
  HitVerifier(HitVerifier& hv) = delete;
  HitVerifier(HitVerifier&& hv) = delete;
  ~HitVerifier() = default;
  void setFastaVector(shared_ptr<FastaVector>& fastaVector);
  void setPhmm(shared_ptr<P7HmmList>& phmmList);
  std::tuple<std::shared_ptr<std::vector<SsvVerifiedHit>>, 
    std::shared_ptr<std::vector<HavacHardwareHit>>> 
    verify(std::shared_ptr<std::vector<HavacHardwareHitReport>> hits);
  void getVerifiedHits();
  void getRejectedHits();
  void clear();


protected:
  shared_ptr<FastaVector> fastaVector;
  shared_ptr<P7HmmList> phmmList;
  unique_ptr<float[]> ssvCellScores;
  shared_ptr<std::vector<HavacHardwareHit> rejectedHits;

  bool verifyHit(const struct HavacHardwareHit& hardwareHitReport);
  int getPhmmIndexFromPosition(const struct HavacHardwareHit& hardwareHitReport);
  void verifyWithReferenceSsv(const uint32_t hitLocatedInPhmmNumber, const uint32_t sequenceNumber,
    const uint32_t localSsvPhmmStartPosition, const uint32_t localPhmmSsvLength,
    const uint32_t ssvSequenceRangeBegin, const uint32_t ssvSequenceRangeEnd);

};
#endif
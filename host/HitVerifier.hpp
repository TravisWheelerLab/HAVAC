#ifndef HAVAC_HIT_VERIFIER_HPP
#define HAVAC_HIT_VERIFIER_HPP


#include <boost/optional.hpp>
#include "HavacSSV.hpp"
#include "Havac.hpp"
extern "C" {
  #include "FastaVector/src/FastaVector.h"
  #include "P7HmmReader/src/p7HmmReader.h"
}

struct SsvReferenceBounds {
  uint32_t sequenceStartIndex;
  uint32_t sequenceEndIndex;
  uint32_t phmmStartIndex;
  uint32_t phmmEndIndex;
  uint32_t phmmNumInList;
};

class HitVerifier {
public:
  HitVerifier(struct FastaVector& fastaVector, struct P7HmmList& phmmList);
  HitVerifier(HitVerifier& hv) = delete;
  HitVerifier(HitVerifier&& hv) = delete;

  void verify(boost::optional<std::vector<struct HavacHostHitReport>>& hits);
  void setReferences(struct FastaVector &fastaVector, struct P7HmmList &phmmList);  //automatically invoked via the constructor.
  //clears the hit locations vector
  void clear();


protected:
  struct FastaVector fastaVector;
  struct P7HmmList phmmList;
  std::vector<struct HavacHitLocation> actualHitLocations;

  struct SsvReferenceBounds &getBoundsForReferenceSsv(const uint32_t sequencePassIndex, const uint32_t sequenceGroupIndex,
    uint32_t phmmIndex);
  bool checkBoundedSsv(const struct SsvReferenceBounds &bounds);
};
#endif
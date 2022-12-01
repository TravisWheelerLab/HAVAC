#include "HitVerifier.hpp"
#include "device/PublicDefines.h"
#include <exception>
#include <stdexcept>

HitVerifier::HitVerifier(struct FastaVector& fastaVector, struct P7HmmList& phmmList) {
  setReferences(fastaVector, phmmList);
}


void HitVerifier::setReferences(struct FastaVector& fastaVector, struct P7HmmList& phmmList) {
  this->fastaVector = fastaVector;
  this->phmmList = phmmList;
}


void HitVerifier::clear() {
  this->actualHitLocations.clear();
}


//can throw std::domain_error
//sequence domain may straddle multiple sequences in fastaVector
struct SsvReferenceBounds& HitVerifier::getBoundsForReferenceSsv(const uint32_t sequencePassIndex, const uint32_t sequenceGroupIndex,
  uint32_t phmmPosition) {
  const uint8_t separatorVectorsPerPhmm = 2;
  const uint32_t maximumPhmmBoundRange = 64;
  struct SsvReferenceBounds refBounds;

  //find the position in the phmm, and which phmm it's in
  uint32_t phmmNumber = 0;
  uint32_t positionInPhmm = 0;
  uint32_t lastPhmmEndPosition = 0;
  bool foundInPhmmRange = false;  //used to verify that the position is actually in range of the phmm list.

  for (uint32_t i = 0; i < this->phmmList.count; i++) {
    uint32_t numSeperators = (i + 1) * separatorVectorsPerPhmm; //there's 2 seperator vectors appended for each phmm
    uint32_t endPositionWithSeperators = this->phmmList.phmms[i].header.modelLength + numSeperators;
    if (phmmPosition < endPositionWithSeperators) {
      phmmNumber = i;
      positionInPhmm = phmmPosition - lastPhmmEndPosition;
      foundInPhmmRange = true;
      break;
    }
    else {
      lastPhmmEndPosition = endPositionWithSeperators;
    }
  }
  if (!foundInPhmmRange) {
    char buffer[1024];
    snprintf(buffer, sizeof(buffer), "phmmLocation %u was past expected maximum range of %u\n", phmmPosition, lastPhmmEndPosition);
    throw new std::domain_error(buffer);
  }

  refBounds.phmmEndIndex = positionInPhmm;
  refBounds.phmmNumInList = phmmNumber;


  //find the beginning bound for the phmm range
  if (positionInPhmm < maximumPhmmBoundRange) {
    refBounds.phmmStartIndex = 0;
  }
  else {
    refBounds.phmmStartIndex = positionInPhmm - maximumPhmmBoundRange;
  }


  //make the sequence start and end positions
  uint64_t sequenceDomainStartPosition = (sequencePassIndex * NUM_CELL_PROCESSORS) + (sequenceGroupIndex * CELLS_PER_GROUP);
  uint64_t sequenceDomainEndPosition = sequenceDomainStartPosition + CELLS_PER_GROUP;
  refBounds.sequenceStartIndex = (sequencePassIndex * NUM_CELL_PROCESSORS) + (sequenceGroupIndex * CELLS_PER_GROUP);
  refBounds.sequenceEndIndex = refBounds.sequenceStartIndex + CELLS_PER_GROUP;

  return refBounds;
}


bool HitVerifier::checkBoundedSsv(const struct SsvReferenceBounds& bounds) {
  //in this function, note that HAVAC should have seperators 'A's at the end of every sequence.
  //this helps drop the score between sequences (in 3/4ths of all cases), and makes the 
  //positions in FastaVector align with the HAVAC positions for sequences,
  //since FastaVector separates sequences with null terminators.
  uint32_t startingSequenceNumber = 0;
  uint32_t endingSequenceNumber = 0;
  uint32_t thisSequenceStartPosition = 0;
  for (uint32_t i = 0; i < this->fastaVector.metadata.count; i++) {
    //use the bounds as the total position inside the fastavector sequence
    //perform ssv against the phmm, with the null terminator ('/0') zeroing out the score 
    //any cells that pass threshold, determine what sequence they fall into 
  }
}
#include "SequencePreprocessor.hpp"
#include "../../../device/PublicDefines.h"

#include <cstdlib>

SequencePreprocessor::SequencePreprocessor(FastaVector* fastaVector) {
  const uint32_t numSequencesInFasta = fastaVector->metadata.count;
  this->originalSequenceLength = fastaVector->sequence.count;
  //round up to the next sequence segment boundary
  this->compressedSequenceLengthInSegments = (originalSequenceLength + (NUM_CELL_PROCESSORS - 1)) / NUM_CELL_PROCESSORS;
  this->compressedSequenceLengthInSymbols = compressedSequenceLengthInSegments * NUM_CELL_PROCESSORS;
  this->compressedSequenceLengthInBytes = compressedSequenceLengthInSymbols / SEQUENCE_PREPROCESSOR_SYMBOLS_PER_BYTE;
  this->compressedSequenceBuffer = std::make_shared<uint8_t[]>(this->compressedSequenceLengthInBytes);

  this->setCompressedSequence(fastaVector);
}

shared_ptr<uint8_t[]> SequencePreprocessor::getCompressedSequenceBuffer() {
  return this->compressedSequenceBuffer;
}

uint32_t SequencePreprocessor::getCompressedSequenceLengthInSegments() {
  return this->compressedSequenceLengthInSegments;
}

uint32_t SequencePreprocessor::getCompressedSequenceLengthInSymbols() {
  return this->compressedSequenceLengthInSymbols;
}

uint32_t SequencePreprocessor::getCompressedSquenceLengthInBytes() {
  return this->compressedSequenceLengthInBytes;
}

void SequencePreprocessor::setCompressedSequence(struct FastaVector* fastaVector) {
  char* sequence = fastaVector->sequence.charData;
  for (uint32_t i = 0; i < this->compressedSequenceLengthInSymbols; i++) {

    //if the symbol is past the end of the original sequence, just write 'a's until the end.
    const char originalSymbol = (i < this->originalSequenceLength) ? sequence[i] : 'a';
    const uint8_t compressedSymbol = this->getCompressedSymbol(originalSymbol);

    const uint32_t byteInCompressedSequence = i / SEQUENCE_PREPROCESSOR_SYMBOLS_PER_BYTE;
    const uint8_t bitIndexInByte = (i % SEQUENCE_PREPROCESSOR_SYMBOLS_PER_BYTE) * 2;

    //mask away whatever was there in the compressed sequence
    const uint8_t bitmask = ~(0x3 << bitIndexInByte);
    (*this->compressedSequenceBuffer)[byteInCompressedSequence] &= bitmask;

    const uint8_t shiftedCompressedSymbol = compressedSymbol << bitIndexInByte;
    (*this->compressedSequenceBuffer)[byteInCompressedSequence] = shiftedCompressedSymbol;
  }
}


uint8_t SequencePreprocessor::getCompressedSymbol(const char originalSymbol) {
  //this function can handle any singe nucleotide, or any ambiguity codes between 2 nucs.
  //Anything else (3 possible nucs or 'N' unknowns) get resolved randomly
  //this first switch should handle all the common symbols
  switch (originalSymbol) {
  case 'a': case 'A': return 0;
  case 'c': case 'C': return 1;
  case 'g': case 'G': return 2;
  case 't': case 'T': return 3;
  default:
    //if it wasn't one of the common symbols, it'll need some amount of randomization
    const uint8_t randMod2 = rand() % 2;
    switch (originalSymbol) {
    case 'r': case 'R': return randMod2 << 1;     //A or G (0 or 2)
    case 'y': case 'Y': return randMod2 << 1 + 1; //C or T (1 or 3)
    case 's': case 'S': return randMod2 + 1;      //C or G (1 or 2)
    case 'w': case 'W': return randMod2 * 3;      //A or T (0 or 3)
    case 'k': case 'K': return randMod2 + 2;      //G or T (2 or 3)
    case 'm': case 'M': return randMod2;          //A or C (0 or 1)
    default: return rand() % 4;
    }
  }
}
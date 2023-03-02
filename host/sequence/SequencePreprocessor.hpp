#ifndef HAVAC_SEQUENCE_PREPROCESSOR_HPP
#define HAVAC_SEQUENCE_PREPROCESSOR_HPP

#include <cstdint>
#include <memory>
#include <vector>

extern "C" {
  #include <FastaVector.h>
}

using std::shared_ptr;
using std::vector;

#define SEQUENCE_PREPROCESSOR_SYMBOLS_PER_BYTE 4

class SequencePreprocessor {
public:
  /// @brief preprocesses the ascii data inside a FastaVector into a bit-compressed
  ///   format that HAVAC can use.  
  /// NOTE: objects from this class should not be reused. If you want to process another fasta, 
  ///     create a new object of this class.
  /// @param fastaVector struct containing the sequences to be processed
  SequencePreprocessor(struct FastaVector* fastaVector);

  /// @brief gets the compressed data representing the sequences inside the fasta from the constructor 
  /// @return shared_ptr to the compressed sequence data
  shared_ptr<vector<uint8_t>> getCompressedSequenceBuffer();

  /// @brief gets the total sequence length in segments. This is a required argument for the HAVAC invocation.
  /// @return number of sequence segments making up the full compressed sequence 
  uint32_t getCompressedSequenceLengthInSegments();

  /// @brief gets the total number of symbols from the entire concatenated sequence 
  /// @return number of symbols represented in the compressed sequence buffer
  uint32_t getCompressedSequenceLengthInSymbols();

  /// @brief gets the length of the compressed sequence, in bytes
  /// @return the length of the compressed sequence, in bytes
  uint32_t getCompressedSquenceLengthInBytes();

private:
  uint32_t originalSequenceLength;
  uint32_t compressedSequenceLengthInSegments;
  uint32_t compressedSequenceLengthInSymbols;
  uint32_t compressedSequenceLengthInBytes;
  shared_ptr<vector<uint8_t>> compressedSequenceBuffer;

  void setCompressedSequence(struct FastaVector* fastaVector);
  uint8_t getCompressedSymbol(const char originalSymbol);
};

#endif
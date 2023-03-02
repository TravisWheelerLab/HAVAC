#ifndef HAVAC_PHMM_PREPROCESSOR_HPP
#define HAVAC_PHMM_PREPROCESSOR_HPP


#include <cstdint>
#include <memory>
#include <vector>
extern "C" {
  #include <p7HmmReader.h>
}

using std::shared_ptr;
using std::vector;

class PhmmPreprocessor {
public:
  /// @brief Class to preprocess a phmm in P7HmmList form into the int8_t* form required by HAVAC 
  /// @param phmmList loaded phmm list to be processed into the phmm data for a single HAVAC call
  /// @param desiredPvalue pvalue used to scale the the scores in the phmm 
  PhmmPreprocessor(std::shared_ptr<P7HmmList> phmmList, const float desiredPvalue = 0.05f);
  
  /// @brief get the processed phmm data that can be loaded onto the HAVAC FPGA system. 
  /// @return shared_ptr to the data array.
  shared_ptr<vector<int8_t>> getProcessedPhmmData();

  /// @brief gets the length of the processed phmm data array
  /// @return length of the ray in bytes
  uint32_t getPhmmLengthInBytes();
  
  /// @brief thesum of the model lengths for all phmms in the phmm list given in the constructor. 
  /// @return the sum of lengths of all models
  uint32_t getPhmmListLengthInVectors();

private:
  shared_ptr<vector<int8_t>> phmmData;
  uint32_t phmmDataLengthInBytes;
  uint32_t phmmDataLengthInVectors;

};

#endif
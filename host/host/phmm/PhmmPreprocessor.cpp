#include "PhmmPreprocessor.hpp"

#include "../../../PhmmReprojection/PhmmReprojection.h"

#define BYTES_PER_PHMM_VECTOR 4

PhmmPreprocessor::PhmmPreprocessor(struct P7HmmList* phmmList, const float desiredPvalue) {
  //concatenates the lengths of all vectors in the phmm file.
  this->phmmDataLengthInVectors = 0;
  for (uint32_t i = 0; i < phmmList->count; i++) {
    this->phmmDataLengthInVectors += phmmList->phmms[i].header.modelLength;
  }

  this->phmmDataLengthInBytes = this->phmmDataLengthInVectors * BYTES_PER_PHMM_VECTOR;
  this->phmmData = std::make_shared<int8_t[]>(this->phmmDataLengthInBytes);

  //individually projects each model such that a hit is generated on a score of 256,
  //and sets the data for the phmmData member variable.
  uint32_t currentIndexIntoProcessedPhmmData = 0;
  for(uint32_t i = 0; i < phmmList->count; i++){
    const uint32_t thisModelLength = phmmList->phmms[i].header.modelLength;
    p7HmmProjectForThreshold256(&phmmList->phmms[i], desiredPvalue, &(*this->phmmData)[currentIndexIntoProcessedPhmmData]);
    currentIndexIntoProcessedPhmmData += thisModelLength * BYTES_PER_PHMM_VECTOR;
  }
}

shared_ptr<int8_t[]> PhmmPreprocessor::getProcessedPhmmData() {
  return this->phmmData;
}

uint32_t PhmmPreprocessor::getPhmmLengthInBytes(){
  return this->phmmDataLengthInBytes;
}

uint32_t PhmmPreprocessor::getPhmmListLengthInVectors(){
  return this->phmmDataLengthInVectors;
}
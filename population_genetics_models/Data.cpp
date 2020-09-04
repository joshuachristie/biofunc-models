#include <vector>
#include "Data.h"

void DataPersistenceInfiniteApprox::set_infinite_approx_persistence(bool persistence){
  _persistence = persistence;
}

void DataPersistenceInfiniteApprox::print_infinite_approx_persistence(){
  // code
}

void DataPersistenceMultipleGens::append_persistence_by_gen(bool persistence){
  _persistence_by_gen.push_back(persistence);
}

void DataPersistenceMultipleGens::print_persistence_by_gen(){
  // code
}

void DataAlleleFreq::append_allele_A_freq_by_gen(double allele_A_freq){
  _allele_A_freq_by_gen.push_back(allele_A_freq);
}

void DataAlleleFreq::print_allele_A_freq_by_gen(){
  // code
}


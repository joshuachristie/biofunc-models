#include <vector>
#include "Data.h"

void DataPersistenceByGen::append_persistence_by_gen(bool persistence){
  _persistence_by_gen.push_back(persistence);
}

void DataAlleleFreq::append_allele_A_freq_by_gen(double allele_A_freq){
  _allele_A_freq_by_gen.push_back(allele_A_freq);
}

void DataAlleleFreqAndPGB::append_allele_A_freq_by_gen(double allele_A_freq){
  _allele_A_freq_by_gen.push_back(allele_A_freq);
}

#include <vector>
#include "Data.h"

void Data::append_persistence_by_gen(bool persistence){
  _persistence_by_gen.push_back(persistence);
}

void Data::append_allele_A_freq_by_gen(double allele_A_freq){
  _allele_A_freq_by_gen.push_back(allele_A_freq);
}

#include <vector>
#include "Data.h"

void Data::append_conditional_existence_by_gen(bool conditional_existence){
  _conditional_existence_by_gen.push_back(conditional_existence);
}

void Data::append_trait_freq_by_gen(double trait_freq){
  _trait_freq_by_gen.push_back(trait_freq);
}

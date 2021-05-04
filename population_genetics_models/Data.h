#ifndef DATA_H
#define DATA_H

#include <vector>

class Data {
public:
  bool _conditional_existence;
  std::vector<bool> _conditional_existence_by_gen;
  std::vector<double> _trait_freq_by_gen;

  void append_conditional_existence_by_gen(bool);
  void append_trait_freq_by_gen(double);
  
  Data(int number_gens_to_record_pp, int reserve_length_af) {
    _conditional_existence_by_gen.reserve(number_gens_to_record_pp);
    _trait_freq_by_gen.reserve(reserve_length_af);
  }
    
};

#endif

#ifndef DATA_CONTAINER_H
#define DATA_CONTAINER_H

#include <vector>
#include "Data.h"

class DataContainer {
public:
  std::vector<Data> _simulation_data;

  const double get_conditional_existence_infinite_approx();
  void set_conditional_existence_outcome_infinite(const int replicate, const bool outcome);
  const std::vector<double> get_conditional_existence_by_gen();
  const void append_conditional_existence(const int replicate, const bool conditional_existence);
  const std::vector<double> get_trait_freqs(const int replicate);
  void append_trait_freq(const int replicate, const double trait_freq);

  DataContainer(int number_replicates, int number_gens_to_record_pp, int reserve_length_af){
    for (int i = 0; i < number_replicates; i++){
      _simulation_data.push_back(Data(number_gens_to_record_pp, reserve_length_af));
    }
  }

};

#endif

#ifndef DATA_CONTAINERS_H
#define DATA_CONTAINERS_H

#include <vector>
#include "Data.h"

class DataContainer {
public:
  std::vector<Data> _simulation_data;

  const double get_persistence_infinite_approx();
  void set_persistence_outcome_infinite(const int replicate, const bool outcome);
  const std::vector<double> get_persistence_by_gen();
  const void append_persistence(const int replicate, const bool persistence);
  const std::vector<double> get_allele_A_freqs(const int replicate);
  void append_allele_A_freq(const int replicate, const double allele_A_freq);

  DataContainer(int number_replicates, int number_gens_to_record_pp, int reserve_length_af){
    for (int i = 0; i < number_replicates; i++){
      _simulation_data.push_back(Data(number_gens_to_record_pp, reserve_length_af));
    }
  }

};

#endif

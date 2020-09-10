#include "DataContainer.h"

const double DataContainer::get_persistence_infinite_approx(){
  int persistence_count = 0;
  for (auto it = _simulation_data.begin(); it != _simulation_data.end(); it++){
    persistence_count += it->_persistence;
  }
  return static_cast<double>(persistence_count) / static_cast<double>(_simulation_data.size());
};

void DataContainer::set_persistence_outcome_infinite(const int replicate, const bool outcome){
  _simulation_data[replicate]._persistence = outcome;
}

const std::vector<double> DataContainer::get_persistence_by_gen(){
  std::vector<double> persistence_probability(_simulation_data[0]._persistence_by_gen.size(), 0.0);
  for (std::size_t i = 0; i < _simulation_data.size(); i++){ // length: number_replicates
    for (std::size_t j = 0; j < persistence_probability.size(); j++){ // length: number_gens_to_output_pp
      persistence_probability[j] += (static_cast<double>(_simulation_data[i]._persistence_by_gen[j]) /
				     static_cast<double>(_simulation_data.size()));
    }
  }
  return persistence_probability;
}

const void DataContainer::append_persistence(const int replicate, const bool persistence){
  _simulation_data[replicate].append_persistence_by_gen(persistence);
}

const std::vector<double> DataContainer::get_allele_A_freqs(const int replicate){
  return _simulation_data[replicate]._allele_A_freq_by_gen;
}
  
void DataContainer::append_allele_A_freq(const int replicate, const double allele_A_freq){
  _simulation_data[replicate].append_allele_A_freq_by_gen(allele_A_freq);
}

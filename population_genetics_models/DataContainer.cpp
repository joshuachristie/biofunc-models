#include "DataContainer.h"

const double DataContainer::get_conditional_existence_infinite_approx(){
  int conditional_existence_count = 0;
  for (auto it = _simulation_data.begin(); it != _simulation_data.end(); it++){
    conditional_existence_count += it->_conditional_existence;
  }
  return static_cast<double>(conditional_existence_count) / static_cast<double>(_simulation_data.size());
};

void DataContainer::set_conditional_existence_outcome_infinite(const int replicate, const bool outcome){
  _simulation_data[replicate]._conditional_existence = outcome;
}

const std::vector<double> DataContainer::get_conditional_existence_by_gen(){
  std::vector<double> conditional_existence_probability(_simulation_data[0]._conditional_existence_by_gen.size(), 0.0);
  for (std::size_t i = 0; i < _simulation_data.size(); i++){ // length: number_replicates
    for (std::size_t j = 0; j < conditional_existence_probability.size(); j++){ // length: number_gens_to_output_pp
      conditional_existence_probability[j] +=
	(static_cast<double>(_simulation_data[i]._conditional_existence_by_gen[j]) /
				     static_cast<double>(_simulation_data.size()));
    }
  }
  return conditional_existence_probability;
}

const void DataContainer::append_conditional_existence(const int replicate, const bool conditional_existence){
  _simulation_data[replicate].append_conditional_existence_by_gen(conditional_existence);
}

const std::vector<double> DataContainer::get_trait_freqs(const int replicate){
  return _simulation_data[replicate]._trait_freq_by_gen;
}
  
void DataContainer::append_trait_freq(const int replicate, const double trait_freq){
  _simulation_data[replicate].append_trait_freq_by_gen(trait_freq);
}

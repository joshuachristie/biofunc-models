#ifndef HTE_H
#define HTE_H

#include <vector>
#include <random>
#include "Parameters.h"

namespace HTE {

  HTE_Model_Parameters parse_parameter_values(int argc, char* argv[]);
  
  std::vector<double> get_fitness_function(const HTE_Model_Parameters &parameters);
  
  void expected_allele_freqs(double &allele_A_freq, const std::vector<double> &haploid_fitnesses,
			     int env_state);
  
  void realised_allele_freqs(double &allele_A_freq, const HTE_Model_Parameters &parameters, std::mt19937 &rng);
  
  void run_simulation(double &allele_A_freq, const std::vector<double> &haploid_fitnesses,
		      const HTE_Model_Parameters &parameters, std::mt19937 &rng);
  
}

#endif

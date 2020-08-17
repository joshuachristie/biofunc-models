#ifndef HSE_H
#define HSE_H

#include <vector>
#include <random>
#include "Parameters.h"

namespace HSE{
  
  HSE_Model_Parameters parse_parameter_values(int argc, char* argv[]);
 
  std::vector<double> get_fitness_function(const HSE_Model_Parameters &parameters);
 
  void expected_allele_freqs(double &allele_A_freq, const std::vector<double> &haploid_fitnesses);
 
  void realised_allele_freqs(double &allele_A_freq, const HSE_Model_Parameters &parameters, std::mt19937 &rng);
 
  void run_simulation(const std::vector<double> &haploid_fitnesses, const HSE_Model_Parameters &parameters,
		      std::mt19937 &rng, std::vector<bool> &final_A_freqs);
  
  void run_model(int argc, char* argv[]);

}

#endif

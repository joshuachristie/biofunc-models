#ifndef HTEOE_H
#define HTEOE_H

#include <vector>
#include <random>
#include "Parameters.h"

namespace HTEOE {
  
  HTEOE_Model_Parameters parse_parameter_values(int argc, char* argv[]);
  
  std::vector<double> get_fitness_function(const HTEOE_Model_Parameters &parameters);
  
  void expected_allele_freqs(double &allele_A_freq, const std::vector<double> &haploid_fitnesses);
  
  void realised_allele_freqs(double &allele_A_freq, const HTEOE_Model_Parameters &parameters, std::mt19937 &rng);
  
  void run_simulation(double &allele_A_freq, const std::vector<double> &haploid_fitnesses,
		      const HTEOE_Model_Parameters &parameters, std::mt19937 &rng);

  void run_model(int argc, char* argv[]);

}

#endif

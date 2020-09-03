#ifndef HTE_H
#define HTE_H

#include <vector>
#include <random>
#include "Parameters.h"

namespace HTE {

  const HTE_Model_Parameters parse_parameter_values(int argc, char* argv[]);
  
  const std::vector<double> get_fitness_function(const HTE_Model_Parameters &parameters);

  void calculate_allele_freqs(double &allele_A_freq, const std::vector<double> &fitnesses,
			      const HTE_Model_Parameters &parameters, std::mt19937 &rng, int &gen);
  
  void run_model(int argc, char* argv[]);

}

#endif

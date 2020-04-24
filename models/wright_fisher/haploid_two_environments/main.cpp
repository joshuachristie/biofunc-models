// Wright-fisher haploid model (one locus with two alleles; one environment with two types)
// At any given time, there is one environment. The environment can change during the simulation however.
// Each allele has two components to its fitness: fitness in env 1 and env 2
// Requires specification of a particular selection regime (i.e. environment at each generation)
// So instead of transition probabilities, there are parameters for the length of env 1 and env 2
// Models the evolved for/maintained by distinction


#include <iostream>
#include <vector>
#include <numeric>
#include "functions.h"
#include "../include/rng.h"
#include "../include/helper_functions.h"
int main(int argc, char* argv[]){
  // parameters
  const double selection_coefficient_A_env_1 = atof(argv[1]);
  const double selection_coefficient_A_env_2 = atof(argv[2]);
  const double selection_coefficient_a_env_1 = atof(argv[3]);
  const double selection_coefficient_a_env_2 = atof(argv[4]);
  const int population_size = atoi(argv[5]);
  const int gen_env_1 = atoi(argv[6]);
  const int gen_env_2 = atoi(argv[7]);
  const int number_replicates = atoi(argv[8]);

  const int number_generations = gen_env_1 + gen_env_2
  // set rng
  std::mt19937 rng = initialiseRNG();
  const std::vector<double> haploid_fitnesses = getFitnessFunction(selection_coefficient_A_env_1,
								   selection_coefficient_A_env_2,
								   selection_coefficient_a_env_1,
								   selection_coefficient_a_env_2);
  std::vector<bool> final_A_freqs;
  // iterate over generations and replicates
  for (int rep = 0; rep < number_replicates; rep++){
    const double tolerance = 0.000001; // for double comparison
    double allele_A_freq = 1.0 / static_cast<double>(population_size); // initial freq is 1/N
    iterateOverGenerations(allele_A_freq, haploid_fitnesses, population_size, number_generations, rng,
			   tolerance);
    // record 0 if extinct, 1 otherwise (indicating the allele exists at a non-zero proportion)
    closeToValue(allele_A_freq, 0.0, tolerance) ? final_A_freqs.push_back(0) : final_A_freqs.push_back(1);
  }
  // print fixation probability
  std::cout << std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) / static_cast<double>(number_replicates) << std::endl;
  return 0;
}

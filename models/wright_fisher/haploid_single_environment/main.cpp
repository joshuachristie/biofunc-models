// This model explores biological function in the simplest case: a single trait with a single effect
// The set-up is a Wright-Fisher haploid model: one locus, two alleles, single environment
// There is no need to apportion the function metric between components using PID because there is a simple
// mapping between trait -> selection_coefficient -> persistence_probability (i.e. the function can be wholly
// assigned to the single effect of the trait)

#include <iostream>
#include <vector>
#include <numeric>
#include "functions.h"
#include "../include/rng.h"
#include "../include/helper_functions.h"
int main(int argc, char* argv[]){
  // parameters
  const double selection_coefficient = atof(argv[1]);
  const int population_size = atoi(argv[2]);
  const int number_generations = atoi(argv[3]);
  const int number_replicates = atoi(argv[4]);
  // set rng
  std::mt19937 rng = initialiseRNG();
  const std::vector<double> haploid_fitnesses = getFitnessFunction(selection_coefficient);
  std::vector<bool> final_A_freqs;
  // iterate over generations and replicates
  for (int rep = 0; rep < number_replicates; rep++){
    const double tolerance = 1.0 / static_cast<double>(population_size * 2); // for double comparison
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

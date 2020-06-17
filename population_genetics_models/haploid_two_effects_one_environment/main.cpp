// This model explores biological function in the context of a trait that causes two phenotypic (fitness) effects
// The step-up is a Wright-Fisher haploid model: two locus (each with two alleles), single environment.
// (But as it is a haploid model, and there is no recombination between the loci, and the two loci are linked.)

// Let's say that the resident population has trait (allele) a which is invaded by our allele of interest, A
// Trait A has two effects, given by selection_coefficient_A1 and selection_coefficient_A2
// The resident trait, a, also has two effects, selection_coefficient_a1 and selection_coefficient_a2
// The fitness of the trait of interest is  wA = 1.0 + selection_coefficient_A1 + selection_coefficient_A2
// The fitness of the resident trait is  wa = 1.0 + selection_coefficient_a1 + selection_coefficient_a2
// For PID purposes, I will fix the fitness of trait a (and thus selection_coefficient_a1 and
// selection_coefficient_a2) systematically varying selection_coefficient_A1 and selection_coefficient_A2.
// PID is used to apportion trait A's function between selection_coefficient_A1 and selection_coefficient_A2
// (and to their redundant/synergistic effects)

#include <iostream>
#include <vector>
#include <numeric>
#include "functions.h"
#include "../include/rng.h"
#include "../include/helper_functions.h"

int main(int argc, char* argv[]){
  // parameters
  const double selection_coefficient_A1 = atof(argv[1]);
  const double selection_coefficient_A2 = atof(argv[2]);
  const double selection_coefficient_a1 = atof(argv[3]);
  const double selection_coefficient_a2 = atof(argv[4]);
  const int population_size = atoi(argv[5]);
  const int number_generations = atoi(argv[6]);
  const int number_replicates = atoi(argv[7]);
  // set rng
  std::mt19937 rng = initialiseRNG();
  const std::vector<double> haploid_fitnesses = getFitnessFunction(selection_coefficient_A1,
								   selection_coefficient_A2,
								   selection_coefficient_a1,
								   selection_coefficient_a2);
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

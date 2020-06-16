// This model explores biological function in the context of a trait that experiences two different environments
// and has a separate fitness in each environment.
// The step-up is a Wright-Fisher haploid model: one locus with two alleles; one environment with two types

// At any given time, there is one environment. The environment can change during the simulation however.
// Let's say that the resident population has trait (allele) a which is invaded by our allele of interest, A
// Each allele (A/a) has two components to its fitness: fitness in env 1 and env 2
// Trait/allele a is the resident trait and in environment 1 its fitness is 1 + selection_coefficient_a_env_1
// while in environment 2 its fitness is 1 + selection_coefficient_a_env_2
// Trait/allele A is our trait of interest and in environment 1 its fitness is 1 + selection_coefficient_A_env_1
// while in environment 2 its fitness is 1 + selection_coefficient_A_env_2
// For PID purposes, I will fix the fitness of trait a (and thus also fix selection_coefficient_a_env_1 and
// selection_coefficient_a_env_2).
// The two "sources" are thus selection_coefficient_A_env_1 and selection_coefficient_A_env_2.
// PID is used to apportion trait A's function between selection_coefficient_A_env_1 and
// selection_coefficient_A_env_2 (and to their redundant/synergistic effects).

// Requires specification of a particular selection regime (i.e. environment at each generation)
// So instead of transition probabilities, there are parameters for the length of env 1 and env 2
// For simplicity, I only consider a selection regime in which environment 1 exists for gen_env_1 generations,
// which is followed by environment 2 existing for gen_env_2 generations
// This gets at the evolved for/maintained by distinction (where evolved for is env 1 and maintained by is env 2)

// Note that here PID is not strictly necessary in order to calculate the overall function of trait/allele A
// (since both sources are ultimately derived from trait/allele A), but it is necessary in order to apportion the
// function of trait/allele A between its effects.
// If a trait evolves in environment 1 but is maintained due to its effects in environment 2, we might want to
// quantify the function of the trait's effect in environment 1 separately from its effect in environment 2


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

  // set rng
  std::mt19937 rng = initialiseRNG();
  const std::vector<double> haploid_fitnesses = getFitnessFunction(selection_coefficient_A_env_1,
								   selection_coefficient_A_env_2,
								   selection_coefficient_a_env_1,
								   selection_coefficient_a_env_2);
  std::vector<bool> final_A_freqs;
  // iterate over generations and replicates
  for (int rep = 0; rep < number_replicates; rep++){
    const double tolerance = 1.0 / static_cast<double>(population_size * 2); // for double comparison
    double allele_A_freq = 1.0 / static_cast<double>(population_size); // initial freq is 1/N
    iterateOverGenerations(allele_A_freq, haploid_fitnesses, population_size, gen_env_1, gen_env_2, rng,
			   tolerance);
    // record 0 if extinct, 1 otherwise (indicating the allele exists at a non-zero proportion)
    closeToValue(allele_A_freq, 0.0, tolerance) ? final_A_freqs.push_back(0) : final_A_freqs.push_back(1);
  }
  // print fixation probability
  std::cout << std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) / static_cast<double>(number_replicates) << std::endl;
  return 0;
}

// This model explores biological function in the context of a trait whose spread is influenced, not only by
// its own effects, but also by the effects of another trait.
// The set-up is a Wright-Fisher diploid model: one locus, two alleles, single environment

// Let's say that the resident population has trait (genotype) aa, which is "invaded" by our allele of interest, A (i.e. an aa individual mutates to Aa).
// Let's say that our trait of interest is the heterozygote Aa.
// Whether A allele persists in the population depends not only on the fitness of Aa but also on the fitness of AA
// Since I use allele frequencies---not genotype frequencies---to calculate the function metric, the function of
// allele A needs to be apportioned over the traits (genotypes) Aa and AA (as well as their redundant/synergistic effects).
// For this I use PID, where the two "sources" are the fitnesses of Aa and AA.
// By combining the function metric and the PID, we can calculate the function of traits Aa and AA.

#include <iostream>
#include <vector>
#include <numeric>
#include "functions.h"
#include "../include/rng.h"
#include "../include/helper_functions.h"

int main(int argc, char* argv[]){
  // parameters
  const double selection_coefficient_homozygote = atof(argv[1]);
  const double selection_coefficient_heterozygote = atof(argv[2]);
  const int population_size = atoi(argv[3]);
  const int number_generations = atoi(argv[4]);
  const int number_replicates = atoi(argv[5]);
  // set rng
  std::mt19937 rng = initialiseRNG();
  const std::vector<double> genotype_fitnesses = getFitnessFunction(selection_coefficient_homozygote,
								    selection_coefficient_heterozygote);
  std::vector<bool> final_A_freqs;
  // iterate over generations and replicates
  for (int rep = 0; rep < number_replicates; rep++){
    const double tolerance = 1.0 / static_cast<double>(population_size * 4); // for double comparison
    double allele_A_freq = 1.0 / static_cast<double>(population_size * 2); // initial freq is 1/2N
    iterateOverGenerations(allele_A_freq, genotype_fitnesses, population_size, number_generations, rng,
			   tolerance);
    // record 0 if extinct, 1 otherwise (indicating the allele exists at a non-zero proportion)
    closeToValue(allele_A_freq, 0.0, tolerance) ? final_A_freqs.push_back(0) : final_A_freqs.push_back(1);
  }
  // print fixation probability
  std::cout << std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) / static_cast<double>(number_replicates) << std::endl;
  return 0;
}

#include <random>
#include <numeric>
#include "../include/helper_functions.h"
#include "functions.h"

std::vector<double> getFitnessFunction(const double selection_coefficient_A_env_1,
				       const double selection_coefficient_A_env_2,
				       const double selection_coefficient_a_env_1,
				       const double selection_coefficient_a_env_2){
  // haploid_fitnesses[w_A_1, w_A_2, w_a_1, w_a_2] gives relative fitness of A/a alleles in envs 1/2
  std::vector<double> haploid_fitnesses {1.0 + selection_coefficient_A_env_1,
					 1.0 + selection_coefficient_A_env_2,
					 1.0 + selection_coefficient_a_env_1,
					 1.0 + selection_coefficient_a_env_2};
  return haploid_fitnesses;
}

void expectedAlleleFreqs(double &allele_A_freq, const std::vector<double> &haploid_fitnesses, int env_state){
  // expected frequency of allele A after selection
  std::vector<double> expected_allele_freq_raw(2);
  expected_allele_freq_raw[0] = allele_A_freq * haploid_fitnesses[env_state];
  expected_allele_freq_raw[1] = (1.0 - allele_A_freq) * haploid_fitnesses[2 + env_state];
  // replace allele_A_freq with normalised expectation
  allele_A_freq = expected_allele_freq_raw[0] / (std::accumulate(expected_allele_freq_raw.begin(),
								 expected_allele_freq_raw.end(), 0.0));
}

void realisedAlleleFreqs(double &allele_A_freq, const int population_size, std::mt19937 &rng){
  // realised frequency of allele A after selection
  std::binomial_distribution<int> surviving_As(population_size, allele_A_freq); // use expectation of allele A
  allele_A_freq = static_cast<double>(surviving_As(rng)) / static_cast<double>(population_size); // proportion
}

void iterateOverGenerations(double &allele_A_freq, const std::vector<double> &haploid_fitnesses,
			    const int population_size, const int gen_env_1, const int gen_env_2,
			    std::mt19937 &rng, const double tolerance){
  // iterate through generations until one allele is fixed or the max number of generations is reached
  int gen = 0;
  int env_state = 0; // env 1 = 0; env 2 = 1
  while (gen < (gen_env_1 + gen_env_2) && !(closeToValue(allele_A_freq, 0.0, tolerance) ||
				       closeToValue(allele_A_freq, 1.0, tolerance))){
    expectedAlleleFreqs(allele_A_freq, haploid_fitnesses, env_state);
    realisedAlleleFreqs(allele_A_freq, population_size, rng);
    ++gen;
    if (gen == gen_env_1 - 1) {env_state = 1;}
  }
}

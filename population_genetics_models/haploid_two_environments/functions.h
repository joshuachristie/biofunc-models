#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <random>

// haploid_fitnesses[w_A_1, w_A_2, w_a_1, w_a_2] gives relative fitness of A/a alleles in envs 1/2
std::vector<double> getFitnessFunction(const double selection_coefficient_A_env_1,
				       const double selection_coefficient_A_env_2,
				       const double selection_coefficient_a_env_1,
				       const double selection_coefficient_a_env_2);
// expected frequency of allele A after selection
void expectedAlleleFreqs(double &allele_A_freq, const std::vector<double> &haploid_fitnesses, int env_state);
// realised frequency of allele A after selection
void realisedAlleleFreqs(double &allele_A_freq, const int population_size, std::mt19937 &rng);
// iterate through generations until one allele is fixed or the max number of generations is reached
void iterateOverGenerations(double &allele_A_freq, const std::vector<double> &haploid_fitnesses,
			    const int population_size, const int gen_env_1, const int gen_env_2, std::mt19937 &rng,
			    const double tolerance);
#endif

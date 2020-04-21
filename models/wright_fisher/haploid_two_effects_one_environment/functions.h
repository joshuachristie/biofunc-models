#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <random>

// haploid_fitnesses[w_A, w_a] gives relative fitness of A allele
std::vector<double> getFitnessFunction(const double selection_coefficient_A1,
				       const double selection_coefficient_A2,
				       const double selection_coefficient_a1,
				       const double selection_coefficient_a2);
// expected frequency of allele A after selection
void expectedAlleleFreqs(double &allele_A_freq, const std::vector<double> &haploid_fitnesses);
// realised frequency of allele A after selection
void realisedAlleleFreqs(double &allele_A_freq, const int population_size, std::mt19937 &rng);
// iterate through generations until one allele is fixed or the max number of generations is reached
void iterateOverGenerations(double &allele_A_freq, const std::vector<double> &haploid_fitnesses,
			    const int population_size, const int number_generations, std::mt19937 &rng,
			    const double tolerance);

#endif

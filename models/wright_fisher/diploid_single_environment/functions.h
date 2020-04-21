#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <random>

// diploid_fitnesses[w_AA, w_Aa, w_aa] gives relative fitness of genotypes
std::vector<double> getFitnessFunction(const double selection_coefficient_homozygote,
				       const double selection_coefficient_heterozygote);
// expected frequency of allele A after selection
void expectedAlleleFreqs(double &allele_A_freq, const std::vector<double> genotype_fitnesses);
// realised frequency of allele A after selection
void realisedAlleleFreqs(double &allele_A_freq, const int population_size, std::mt19937 &rng);
// iterate through generations until one allele is fixed or the max number of generations is reached
void iterateOverGenerations(double &allele_A_freq, const std::vector<double> genotype_fitnesses,
			    const int population_size, const int number_generations, std::mt19937 &rng,
			    const double tolerance);
#endif

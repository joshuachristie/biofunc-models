#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <random>

std::vector<double> getFitnessFunction(const double selection_coefficient);
void expectedAlleleFreqs(double &allele_A_freq, std::vector<double> &haploid_fitnesses);
void realisedAlleleFreqs(double &allele_A_freq, const int population_size, std::mt19937 &rng);
void iterateOverGenerations(double &allele_A_freq, const std::vector<double> haploid_fitnesses,
			    const int population_size, const int number_generations);
bool closeToValue(double allele_freq, double value);

#endif

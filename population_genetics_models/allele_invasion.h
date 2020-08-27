/**
   @file allele_invasion.h
*/

#ifndef ALLELE_INVASION_H
#define ALLELE_INVASION_H

#include <random>
#include <vector>
#include "helper_functions.h"

/**
   @brief Runs a single invasion attempt of an allele
   @param[in] fitnesses Vector containing fitnesses
   @param[in] parameters.shared.population_size Number of individuals in the population
   @param[in] parameters.fixed.tolerance Tolerance for comparing equality of doubles
   @param[in, out] rng Random number generator
   @param[in, out] final_A_freqs Vector of bools storing whether allele A persists (true/false) at census
   @param[in] calculate_allele_method Template for method to calcluate allele frequency (one of HSE::calculate_allele_freqs, HTE::calculate_allele_freqs, DSE::calculate_allele_freqs, or HTEOE::calculate_allele_freqs)
   @return Nothing (but alters \p final_A_freqs)
*/

template <class P, class F>
void allele_invasion(const std::vector<double> &fitnesses, const P &parameters, std::mt19937 &rng,
		     std::vector<bool> &final_A_freqs, double &allele_A_freq, F calculate_allele_method){
  int gen = 0;
  while (help::is_neither_fixed_nor_extinct(gen, allele_A_freq, parameters)){
    calculate_allele_method(allele_A_freq, fitnesses, parameters, rng, gen);
  }
  help::close_to_value(allele_A_freq, 0.0, parameters.fixed.tolerance) ? final_A_freqs.push_back(0) :
    final_A_freqs.push_back(1);
}

#endif

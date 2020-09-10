/**
   @file allele_invasion.h
*/

#ifndef ALLELE_INVASION_H
#define ALLELE_INVASION_H

#include <random>
#include <vector>
#include "helper_functions.h"
#include "DataContainers.h"
/**
   @brief Runs a single invasion attempt of an allele
   @param[in] fitnesses Vector containing fitnesses
   @param[in] parameters.shared.population_size Number of individuals in the population
   @param[in] parameters.fixed.tolerance Tolerance for comparing equality of doubles
   @param[in, out] rng Random number generator
   @param[in] calculate_allele_method Template for method to calcluate allele frequency (one of HSE::calculate_allele_freqs, HTE::calculate_allele_freqs, DSE::calculate_allele_freqs, or HTEOE::calculate_allele_freqs)
   @param[in, out] final_A_freqs Vector storing frequencies of the A allele
   @return Nothing (but alters \p allele_A_freq and can alter final_A_freqs depending on parameters.shared.num_gens_to_output_pp)
*/
template <class P, class F>
void allele_invasion(const std::vector<double> &fitnesses, const P &parameters, std::mt19937 &rng,
		     double &allele_A_freq, F calculate_allele_freqs_function, DataContainer &data,
		     const int replicate){
  int gen = -1;
  while (help::is_neither_fixed_nor_extinct(gen, allele_A_freq, parameters) ||
	 gen < parameters.shared.number_gens_to_output_pp){
    calculate_allele_freqs_function(allele_A_freq, fitnesses, parameters, rng, gen);
    if (gen < parameters.shared.number_gens_to_output_pp){
      help::record_A_allele_presence_by_gen(allele_A_freq, parameters, replicate, data);
    }
  }
  help::record_A_allele_presence_infinite(allele_A_freq, parameters, replicate, data);
}

#endif

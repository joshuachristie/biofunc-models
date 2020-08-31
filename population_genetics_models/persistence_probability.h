/**
   @file persistence.probability.h
*/
#ifndef PERSISTENCE_PROBABILITY_H
#define PERSISTENCE_PROBABILITY_H

#include <vector>
#include <random>
#include <numeric>
#include "allele_invasion.h"
#include "helper_functions.h"

/**
   @brief Template function to run replicates and calculate persistence probability for the pop gen models
   @param[in] params Template for HSE_Model_Parameters, DSE_Model_Parameters, HTE_Model_Parameters, or HTEOE_Model_Parameters
   @param[in] fitnesses Vector of allele or genotype fitnesses
   @param[in, out] rng Random number generator
   @param[in, out] final_A_freqs Vector of bools storing whether allele A persists (true/false) at census
   @param[in] calculate_allele_method Template for method to calcluate allele frequency (one of HSE::calculate_allele_freqs, HTE::calculate_allele_freqs, DSE::calculate_allele_freqs, or HTEOE::calculate_allele_freqs)
   @return Persistence_probability The probability that the trait persists in the population
*/
template <class P, class F>
const double calculate_persistence_probability(const P &params, std::mt19937 &rng, const std::vector<double>
					       &fitnesses, std::vector<bool> &final_A_freqs,
					       F calculate_allele_method){
  for (int i = 0; i < params.fixed.number_replicates; i++){
    double allele_A_freq = params.shared.initial_A_freq;
    // run simulation to see whether allele A invades and either becomes fixed or withstands 1000000 gens
    allele_invasion(fitnesses, params, rng, allele_A_freq, calculate_allele_method);
    int reinvasions = 0;
    // run reinvasion attempts by resident while allele A remains (if number_reinvasions is non-zero)
    while (help::is_not_extinct(allele_A_freq, params) && reinvasions < params.shared.number_reinvasions){
      allele_A_freq -= params.shared.initial_A_freq; // replace single A allele with an a allele
      allele_invasion(fitnesses, params, rng, allele_A_freq, calculate_allele_method);
      reinvasions++;
    }
    help::is_not_extinct(allele_A_freq, params) ? final_A_freqs.push_back(1) : final_A_freqs.push_back(0);
  }
  return std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) /
    static_cast<double>(params.fixed.number_replicates);
}

#endif

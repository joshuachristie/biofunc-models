/**
   @file persistence.probability.h
*/
#ifndef PERSISTENCE_PROBABILITY_H
#define PERSISTENCE_PROBABILITY_H

#include <vector>
#include <random>
#include <numeric>
#include "allele_invasion.h"

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
double calculate_persistence_probability(const P &params, std::mt19937 &rng,
					 const std::vector<double> &fitnesses,
					 std::vector<bool> &final_A_freqs, F calculate_allele_method){
  for (int i = 0; i < params.fixed.number_replicates; i++){
    double allele_A_freq = params.model.initial_A_freq;
    allele_invasion(fitnesses, params, rng, final_A_freqs, allele_A_freq, calculate_allele_method);
  }
  return std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) /
    static_cast<double>(params.fixed.number_replicates);
}

#endif

/**
   @file conditional_existence_probability.h
*/
#ifndef CONDITIONAL_EXISTENCE_PROBABILITY_H
#define CONDITIONAL_EXISTENCE_PROBABILITY_H

#include <vector>
#include <random>
#include <numeric>
#include "trait_invasion.h"
#include "conditional_existence_status.h"
#include "DataContainer.h"
#include "trait_freq.h"

/**
   @brief Template function to run replicates and calculate conditional existence probability for the pop gen models
   @param[in] params Template for HSE_Model_Parameters, DSE_Model_Parameters, HTE_Model_Parameters, or HTEOE_Model_Parameters
   @param[in] fitnesses Vector of allele or genotype fitnesses
   @param[in, out] rng Random number generator
   @param[in] calculate_trait_freqs Template for method to calcluate trait frequency (one of HSE::calculate_trait_freqs, HTE::calculate_trait_freqs, DSE::calculate_trait_freqs, or HTEOE::calculate_trait_freqs)
   @param[in, out] data DataContainer class object
   @return Nothing (but modifies \p data)
*/
template <class P, class F>
void calculate_conditional_existence_probability(const P &params, std::mt19937 &rng, const std::vector<double>
						 &fitnesses, F calculate_trait_freqs, DataContainer &data){

  for (int i = 0; i < params.fixed.number_replicates; i++){
    std::vector<double> trait_freq = trait_freq::initialise_trait_freq(params);
    int reinvasions = -1;
    // run simulation to see whether trait invades and either becomes fixed or withstands 1000000 gens
    trait_invasion(fitnesses, params, rng, trait_freq, calculate_trait_freqs, data, i, reinvasions);
    // run reinvasion attempts by resident while trait remains (if number_reinvasions is non-zero)
    while (!conditional_existence_status::trait_extinct(trait_freq, params) &&
	   reinvasions < params.shared.number_reinvasions - 1){

      reinvasions++;
      // replace single individual carrying trait of interest with single individual carrying resident trait
      trait_freq[ params.shared.trait_info[0] ] -= params.shared.initial_trait_freq;
      trait_invasion(fitnesses, params, rng, trait_freq, calculate_trait_freqs, data, i, reinvasions);
    }
  }
}

#endif

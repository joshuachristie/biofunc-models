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
#include "DataContainers.h"
#include "fixed_parameters.h"

/**
   @brief Template function to run replicates and calculate persistence probability for the pop gen models
   @param[in] params Template for HSE_Model_Parameters, DSE_Model_Parameters, HTE_Model_Parameters, or HTEOE_Model_Parameters
   @param[in] fitnesses Vector of allele or genotype fitnesses
   @param[in, out] rng Random number generator
   @param[in] calculate_allele_freqs_function Template for method to calcluate allele frequency (one of HSE::calculate_allele_freqs, HTE::calculate_allele_freqs, DSE::calculate_allele_freqs, or HTEOE::calculate_allele_freqs)
   @param[in, out] data DataContainer class object
   @return Nothing (but modifies \p data)
*/
template <class P, class F>
void calculate_persistence_probability(const P &params, std::mt19937 &rng, const std::vector<double>
				       &fitnesses, F calculate_allele_freqs_function, DataContainer &data){

  for (int i = 0; i < params.fixed.number_replicates; i++){
    double allele_A_freq = params.shared.initial_A_freq;
    int reinvasions = -1;
    // run simulation to see whether allele A invades and either becomes fixed or withstands 1000000 gens
    allele_invasion(fitnesses, params, rng, allele_A_freq, calculate_allele_freqs_function, data, i, reinvasions);
    // run reinvasion attempts by resident while allele A remains (if number_reinvasions is non-zero)
    while (help::is_not_extinct(allele_A_freq, params) && reinvasions < params.shared.number_reinvasions - 1){
      reinvasions++;
      allele_A_freq -= params.shared.initial_A_freq; // replace single A allele with an a allele
      allele_invasion(fitnesses, params, rng, allele_A_freq, calculate_allele_freqs_function, data, i, reinvasions);
    }
  }
}

#endif

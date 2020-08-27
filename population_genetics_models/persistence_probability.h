/**
   @file persistence.probability.h
*/
#ifndef PERSISTENCE_PROBABILITY_H
#define PERSISTENCE_PROBABILITY_H

#include <vector>
#include <random>
#include <numeric>

/**
   @brief Template function to run replicates and calculate persistence probability for the pop gen models
   @param[in] params Template for HSE_Model_Parameters, DSE_Model_Parameters, HTE_Model_Parameters, or HTEOE_Model_Parameters
   @param[in] run_sim Template for HSE::run_simulation, DSE::run_simulation, HTE::run_simulation or HTEOE::run_simulation
   @param[in] fitnesses Either \p haploid_fitnesses (HSE, HTE, HTEOE) or \p genotype_fitnesses (DSE)
   @param[in, out] rng Random number generator
   @param[in, out] final_A_freqs Vector of bools storing whether allele A persists (true/false) at census
   @return Persistence_probability The probability that the trait persists in the population
*/
template <class P, class F>
double calculate_persistence_probability(const P &params, F run_sim, std::mt19937 &rng, const
					 std::vector<double> &fitnesses, std::vector<bool> &final_A_freqs){
  for (int i = 0; i < params.fixed.number_replicates; i++){
    run_sim(fitnesses, params, rng, final_A_freqs);
    double allele_A_freq = params.model.initial_A_freq;
  }
  return std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) /
    static_cast<double>(params.fixed.number_replicates);
}

#endif

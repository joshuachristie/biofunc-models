/**
   @file persistence.probability.h
*/
#ifndef PERSISTENCE_PROBABILITY_H
#define PERSISTENCE_PROBABILITY_H

#include <vector>
#include <random>
#include <numeric>
#include "allele_invasion.h"
#include "persistence_status.h"
#include "DataContainer.h"

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
    double trait_freq = params.shared.initial_trait_freq;
    int reinvasions = -1;
    // run simulation to see whether trait invades and either becomes fixed or withstands 1000000 gens
    allele_invasion(fitnesses, params, rng, trait_freq, calculate_allele_freqs_function, data, i, reinvasions);
    // run reinvasion attempts by resident while allele A remains (if number_reinvasions is non-zero)
    while (persist_status::is_not_extinct(trait_freq, params) && reinvasions < params.shared.number_reinvasions - 1){
      reinvasions++;
       // I need to decide how to handle this. If I'm doing traits then for the diploid case I have 3 types
      // which means that I can't just decrease the frequency. As I thought how to handle this, I realised
      // that I don't want to do these reinvasions for all of the different types of models (I think the HSE
      // is enough to prove the point). So I might need to rework this entire allele_invasion method
      // (i.e. separate it into a regular invasion and a reinvasion)
      // but not sure whether that solves anything--the main issue is that I need trait_freq to be a vector
      // for the DSE. Actually this might be totally fine? I can just make it a vector either of length 1 or 2
      // if the trait of interest is always in the first position, then I can always access it with
      // trait_freq[0]. I only need trait_freq[1] for DSE and calculate allele freqs.
      // probably makes sense that the trait_freq that isn't included is the initial resident (aa for diploid)

      // one potential problem is how to deal with having AA as the trait of interest vs Aa as the trait of
      // interest. If I always refer to trait_freq[0] then I'd need to swap the position of AA and Aa in
      // trait_freq. I think this should be okay; I just need another parameter that indicates
      // which trait is the one of interest and use this to index into trait_freq (it would need to be set
      // for all the models (0 for the haploid variants) but only provided as a command line argument to DSE)
      // I could actually make it a using declaration so I can continue to refer to trait_freq
      trait_freq -= params.shared.initial_trait_freq; // replace single trait with an a allele
      exit (EXIT_FAILURE);
      allele_invasion(fitnesses, params, rng, trait_freq, calculate_allele_freqs_function, data, i, reinvasions);
    }
  }
}

#endif

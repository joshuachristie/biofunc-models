/**
   @file trait_invasion.h
*/

#ifndef TRAIT_INVASION_H
#define TRAIT_INVASION_H

#include <random>
#include <vector>
#include "record_data.h"
#include "persistence_status.h"
#include "DataContainer.h"
/**
   @brief Runs a single invasion attempt of a trait
   @param[in] fitnesses Vector containing fitnesses
   @param[in] parameters.shared.population_size Number of individuals in the population
   @param[in] parameters.fixed.tolerance Tolerance for comparing equality of doubles
   @param[in, out] rng Random number generator
   @param[in] calculate_trait_freqs Template for method to calcluate traitfrequency (one of HSE::calculate_trait_freqs, HTE::calculate_trait_freqs, DSE::calculate_trait_freqs, or HTEOE::calculate_trait_freqs)
   @param[in, out] data DataContainer class object
   @param[in] replicate Simulation number (for indexing into \p data)
   @param[in] reinvasions Equals -1 when the invasion is the initial one (i.e. trait invading resident)
   @return Nothing (but alters \p trait_freq)
*/
template <class P, class F>
void trait_invasion(const std::vector<double> &fitnesses, const P &parameters, std::mt19937 &rng,
		    std::vector<double> &trait_freq, F calculate_trait_freqs,
		    DataContainer &data, const int replicate, const int reinvasions){
  int gen = -1;
  bool extinct = false; // for formatting printing trait raw data

  bool allele_A_extinct, allele_A_fixed, reached_max_gen, output_pp_by_gen;
  do {
    
    calculate_trait_freqs(trait_freq, fitnesses, parameters, rng, gen);

    allele_A_extinct = persist_status::allele_A_extinct(trait_freq, parameters);
    allele_A_fixed = persist_status::allele_A_fixed(trait_freq, parameters);
    reached_max_gen = persist_status::reached_max_gen(gen, parameters);
    output_pp_by_gen = persist_status::output_pp_by_gen(gen, parameters);
    
    if (output_pp_by_gen && reinvasions == -1){
      record::trait_presence_by_gen(trait_freq, parameters, replicate, data);
    }
    
    if ((parameters.shared.print_trait_raw_data && reinvasions == -1 && !extinct)) {
      record::trait_freq(trait_freq, parameters, replicate, data);
    }
    
    if (allele_A_extinct) {
      extinct = true;
    }

  }
  while ( (!allele_A_extinct && !allele_A_fixed && !reached_max_gen) || output_pp_by_gen );
  
  record::trait_presence_infinite(trait_freq, parameters, replicate, data);
}

#endif

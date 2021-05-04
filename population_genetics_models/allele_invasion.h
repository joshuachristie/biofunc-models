/**
   @file allele_invasion.h
*/

#ifndef ALLELE_INVASION_H
#define ALLELE_INVASION_H

#include <random>
#include <vector>
#include "record_data.h"
#include "persistence_status.h"
#include "DataContainer.h"
/**
   @brief Runs a single invasion attempt of an allele
   @param[in] fitnesses Vector containing fitnesses
   @param[in] parameters.shared.population_size Number of individuals in the population
   @param[in] parameters.fixed.tolerance Tolerance for comparing equality of doubles
   @param[in, out] rng Random number generator
   @param[in] calculate_allele_freqs_function Template for method to calcluate allele frequency (one of HSE::calculate_allele_freqs, HTE::calculate_allele_freqs, DSE::calculate_allele_freqs, or HTEOE::calculate_allele_freqs)
   @param[in, out] data DataContainer class object
   @param[in] replicate Simulation number (for indexing into \p data)
   @param[in] reinvasions Equals -1 when the invasion is the initial one (i.e. trait invading resident)
   @return Nothing (but alters \p allele_A_freq and can alter final_A_freqs depending on parameters.shared.num_gens_to_output_pp)
*/
template <class P, class F>
void allele_invasion(const std::vector<double> &fitnesses, const P &parameters, std::mt19937 &rng,
		     std::vector<double> &trait_freq, F calculate_allele_freqs_function,
		     DataContainer &data, const int replicate, const int reinvasions){
  int gen = -1;
  bool extinct = false; // for formatting printing trait raw data

  bool allele_A_extinct, allele_A_fixed, reached_max_gen, output_pp_by_gen;
  do {
    calculate_allele_freqs_function(trait_freq, fitnesses, parameters, rng, gen);

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

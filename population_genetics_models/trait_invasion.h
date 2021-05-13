/**
   @file trait_invasion.h
*/

#ifndef TRAIT_INVASION_H
#define TRAIT_INVASION_H

#include <random>
#include <vector>
#include "conditional_existence_status.h"
#include "include/example.pb.h"
#include "record_data.h"

/**
   @brief Runs a single invasion attempt of a trait
   @param[in] fitnesses Vector containing fitnesses
   @param[in] parameters.shared.population_size Number of individuals in the population
   @param[in] parameters.fixed.tolerance Tolerance for comparing equality of doubles
   @param[in, out] rng Random number generator
   @param[in] calculate_trait_freqs Template for method to calcluate traitfrequency (one of HSE::calculate_trait_freqs, HTE::calculate_trait_freqs, DSE::calculate_trait_freqs, or HTEOE::calculate_trait_freqs)
   @param[in] reinvasions Equals -1 when the invasion is the initial one (i.e. trait invading resident)
   @return Nothing (but alters \p trait_freq)
*/
template <class P, class F>
void trait_invasion(const std::vector<double> &fitnesses, const P &parameters, std::mt19937 &rng,
		    std::vector<double> &trait_freq, F calculate_trait_freqs, int &gen){
  bool allele_A_extinct, allele_A_fixed, reached_max_gen;
  do {
    calculate_trait_freqs(trait_freq, fitnesses, parameters, rng, gen);

    allele_A_extinct = conditional_existence_status::allele_A_extinct(trait_freq, parameters);
    allele_A_fixed = conditional_existence_status::allele_A_fixed(trait_freq, parameters);
    reached_max_gen = conditional_existence_status::reached_max_gen(gen, parameters);
  }
  while ( !allele_A_extinct && !allele_A_fixed && !reached_max_gen );
}
// overloaded method for LSTM scenario
template <class P, class F>
void trait_invasion(const std::vector<double> &fitnesses, const P &parameters, std::mt19937 &rng,
		    std::vector<double> &trait_freq, F calculate_trait_freqs, int &gen,
		    data::FloatList* raw_trait_freq){
  bool allele_A_extinct, allele_A_fixed, reached_max_gen;
  record::raw_trait_freq(raw_trait_freq, trait_freq); // record initial freqs
  do {
    calculate_trait_freqs(trait_freq, fitnesses, parameters, rng, gen);
    record::raw_trait_freq(raw_trait_freq, trait_freq);
    
    allele_A_extinct = conditional_existence_status::allele_A_extinct(trait_freq, parameters);
    allele_A_fixed = conditional_existence_status::allele_A_fixed(trait_freq, parameters);
    reached_max_gen = conditional_existence_status::reached_max_gen(gen, parameters);
  }
  while ( !allele_A_extinct && !allele_A_fixed && !reached_max_gen );
}

#endif

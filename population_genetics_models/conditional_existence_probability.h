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
#include "trait_freq.h"
#include "include/example.pb.h"
#include "record_data.h"

/**
   @brief Template function to run replicates and calculate conditional existence probability for the pop gen models
   @param[in] params Template for HSE_Model_Parameters, DSE_Model_Parameters, HTE_Model_Parameters, or HTEOE_Model_Parameters
   @param[in] fitnesses Vector of allele or genotype fitnesses
   @param[in, out] rng Random number generator
   @param[in] calculate_trait_freqs Template for method to calcluate trait frequency (one of HSE::calculate_trait_freqs, HTE::calculate_trait_freqs, DSE::calculate_trait_freqs, or HTEOE::calculate_trait_freqs)
   @return Nothing (but modifies \p data)
*/
template <class P, class F>
void calculate_conditional_existence_probability(const P &params, std::mt19937 &rng, const std::vector<double>
						 &fitnesses, F calculate_trait_freqs, data::Int64List*
						 gen_extinct, data::Int64List* reinvasion_number){

  for (int i = 0; i < params.fixed.number_replicates_QEF; i++){
    std::vector<double> trait_freq = trait_freq::initialise_trait_freq(params);
    int reinvasions = -1;
    int gen = -1;
    // run simulation to see whether trait invades and either becomes fixed or withstands the max gens
    trait_invasion(fitnesses, params, rng, trait_freq, calculate_trait_freqs, gen);
    // record conditional existence status of trait
    record::generation_trait_extinction(gen_extinct, trait_freq, params, gen);
    // run reinvasion attempts by resident while trait remains (if number_reinvasions is non-zero)
    while (!conditional_existence_status::trait_extinct(trait_freq, params) &&
	   reinvasions < params.shared.number_reinvasions - 1){
      gen = -1;
      reinvasions++;
      // replace single individual carrying trait of interest with single individual carrying resident trait
      trait_freq[ params.shared.trait_info[0] ] -= params.shared.initial_trait_freq;
      // run simulation to see whether trait resists invasion
      trait_invasion(fitnesses, params, rng, trait_freq, calculate_trait_freqs, gen);
    }
    record::number_reinvasions_before_extinction(reinvasion_number, trait_freq, params, reinvasions);
  }
}
// overloaded method for LSTM scenario
template <class P, class F>
void calculate_conditional_existence_probability(const P &params, std::mt19937 &rng, const std::vector<double>
						 &fitnesses, F calculate_trait_freqs, data::Int64List*
						 gen_extinct, data::FeatureList &featurelist){

  for (int i = 0; i < params.fixed.number_replicates_LSTM; i++){
    std::vector<double> trait_freq = trait_freq::initialise_trait_freq(params);
    int gen = -1;
    data::Feature* raw_trait_frequencies = featurelist.add_feature();
    data::FloatList* raw_trait_freq = raw_trait_frequencies->mutable_float_list();
    // run replicate, record raw_trait_freq
    trait_invasion(fitnesses, params, rng, trait_freq, calculate_trait_freqs, gen, raw_trait_freq);
    // record conditional existence status of trait
    record::generation_trait_extinction(gen_extinct, trait_freq, params, gen);
  }

  for (int i = 0; i < params.fixed.number_replicates_QEF - params.fixed.number_replicates_LSTM; i++){
    std::vector<double> trait_freq = trait_freq::initialise_trait_freq(params);
    int gen = -1;
    // run replicate, don't record raw_trait_freq
    trait_invasion(fitnesses, params, rng, trait_freq, calculate_trait_freqs, gen);
    // record conditional existence status of trait
    record::generation_trait_extinction(gen_extinct, trait_freq, params, gen);
  }
  
}

#endif

/**
   @file persistence_status.h
*/
#ifndef PERSIST_STATUS_H
#define PERSIST_STATUS_H

#include <vector>
#include <numeric>

/**
   @brief Namespace for functions that determine whether a trait persists
**/
namespace persist_status {

  bool close_to_value(const double trait_freq, const double value, const double tolerance);
  double get_allele_A_freq(const std::vector<double> &trait_freq);

  /**
     @brief Checks whether allele A freq is equal to 0
     @param[in] trait_freq Vector storing allele A freq (if haploid) or AA and Aa genotype freq (if diploid)
     @param[in] params.fixed.tolerance Tolerance for comparison of doubles
     @return True if allele A freq is equal to 0; false otherwise
  */
  template<class T>
  bool allele_A_extinct(const std::vector<double> &trait_freq, const T &params){
    double allele_A_frequency = get_allele_A_freq(trait_freq);
    bool allele_A_is_extinct = close_to_value(allele_A_frequency, 0.0, params.fixed.tolerance);
    bool allele_A_negative_value = allele_A_frequency < 0.0; // shouldn't happen but to be extra sure
    return allele_A_is_extinct || allele_A_negative_value;
  }
  /**
     @brief Checks whether trait freq is equal to 0
     @param[in] trait_freq Vector storing allele A freq (if haploid) or AA and Aa genotype freq (if diploid)
     @param[in] params.fixed.tolerance Tolerance for comparison of doubles
     @param[in] params.shared.trait_info Contains index of our trait of interest (in \p trait_freq)
     @return True if trait freq is equal to 0; false otherwise
  */
  template<class T>
  bool trait_extinct(const std::vector<double> &trait_freq, const T &params){
    double trait_frequency = trait_freq[params.shared.trait_info[0]];
    bool trait_is_extinct = close_to_value(trait_frequency, 0.0, params.fixed.tolerance);
    bool trait_negative_value = trait_frequency < 0.0; // shouldn't happen but to be extra sure
    return trait_is_extinct || trait_negative_value;
  }
  /**
     @brief Checks whether allele A freq is equal to 1
     @param[in] trait_freq Vector storing allele A freq (if haploid) or AA and Aa genotype freq (if diploid)
     @param[in] params.fixed.tolerance Tolerance for comparison of doubles
     @return True if allele A freq is equal to 1; false otherwise
  */
  template<class T>
  bool allele_A_fixed(const std::vector<double> &trait_freq, const T &params){
    double allele_A_frequency = get_allele_A_freq(trait_freq);
    bool allele_A_is_fixed = close_to_value(allele_A_frequency, 1.0, params.fixed.tolerance);
    bool allele_A_exceeds_one = allele_A_frequency > 1.0; // shouldn't happen but to be extra sure
    return allele_A_is_fixed || allele_A_exceeds_one;
  }
  /**
     @brief Checks whether simulation has reached the max number of generations
     @param[in] gen Current generation of simulation
     @param[in] params.fixed.max_generations_per_sim Maximum number of generations for which a sim can run (approximation to time -> infinity)
     @return True if simulation has reached the max number of generations; false otherwise
  */
  template<class T>
  bool reached_max_gen(const double gen, const T &params){
    return gen >= params.fixed.max_generations_per_sim;
  }
  /**
     @brief Checks whether to output persistence probability by generation
     @param[in] gen Current generation of simulation
     @param[in] params.shared.number_gens_to_output_pp Number of generations to output persistency probability by generation
     @return True if simulation should output persistence probability for this generation; false otherwise
  */
  template<class T>
  bool output_pp_by_gen(const double gen, const T &params){
    return gen < params.shared.number_gens_to_output_pp;
  }

}

#endif 

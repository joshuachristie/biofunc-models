/**
   @file helper_functions.h
*/
#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <vector>

namespace help {

  bool close_to_value(const double allele_freq, const double value, const double tolerance);
  /**
     @brief Checks whether 0 < \p value < 1
     @param[in] gen Current generation of simulation
     @param[in] value Allele frequency
     @param[in] params.fixed.max_generations_per_sim Maximum number of generations for which a sim can run (approximation to time -> infinity)
     @param[in] params.fixed.tolerance Tolerance for comparison of doubles
     @return True if allele A is neither fixed nor extinct (nor has exceeded params.fixed.max_generations_per_sim); false otherwise
  */
  template<class T>
  bool is_neither_fixed_nor_extinct(const double gen, const double value, const T &params){
    return gen < params.fixed.max_generations_per_sim &&
      !(
	close_to_value(value, 0.0, params.fixed.tolerance) || close_to_value(value, 1.0, params.fixed.tolerance)
	);
  }
  /**
     @brief Checks whether \p value > 0
     @param[in] value Allele frequency
     @param[in] params.fixed.tolerance Tolerance for comparison of doubles
     @return True if \p value is > 0; false otherwise
  */
  template<class T>
  bool is_not_extinct(const double value, const T &params){
    return !(close_to_value(value, 0.0, params.fixed.tolerance) || value < 0.0);
  }
  /**
     @brief Records True if A allele is present; False if A allele is extinct
     @param[in] allele_A_freq Frequency of A allele
     @param[in] parameters Parameter struct
     @param[in, out] final_A_freqs Vector to store A allele frequencies in
     @return Nothing (but modifies \p final_A_freqs)
  */
  template<class T>
  void record_A_allele_presence(const double allele_A_freq, const T &parameters,
				std::vector<bool> &final_A_freqs){
    is_not_extinct(allele_A_freq, parameters) ? final_A_freqs.push_back(1) : final_A_freqs.push_back(0);
  }
  
}

#endif 

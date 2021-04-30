/**
   @file persistence_status.h
*/
#ifndef PERSIST_STATUS_H
#define PERSIST_STATUS_H

#include <vector>

/**
   @brief Namespace for functions that determine whether a trait persists
**/
namespace persist_status {

  bool close_to_value(const double trait_freq, const double value, const double tolerance);
  
  /**
     @brief Checks whether 0 < \p value < 1
     @param[in] gen Current generation of simulation
     @param[in] value Trait frequency
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
  
}

#endif 

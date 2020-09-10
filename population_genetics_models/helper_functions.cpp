/**
   @file helper_functions.cpp
   @brief Collection of miscellaneous functions
*/
#include <cmath>
#include "helper_functions.h"
#include "fixed_parameters.h"

/**
   @brief Namespace for helper functions
*/
namespace help {
  /**
     @brief Compare for approximate equality of doubles
     @param[in] allele_freq Allele frequency
     @param[in] value Value to compare \p allele_freq against
     @param[in] tolerance A double that determines how closely \p allele_freq and \p value must be
     @return True if the absolute difference between \p allele_freq and \p value is less than \p tolerance; false otherwise
  */
  bool close_to_value(const double allele_freq, const double value, const double tolerance){
    if (std::abs(allele_freq - value) < tolerance){
      return true;
    } else {
      return false;
    }
  }
  
  void record_A_allele_freq(const double allele_A_freq, const int replicate, DataContainer &data){
    
    data.append_allele_A_freq(replicate, allele_A_freq);
    // (re)reserve a large buffer if the replicate runs for a while (to avoid excessive memory allocations)
    std::size_t current_size = data._simulation_data[replicate]._allele_A_freq_by_gen.size();
    std::size_t current_capacity = data._simulation_data[replicate]._allele_A_freq_by_gen.capacity();
    if (current_size == current_capacity){
      data._simulation_data[replicate]._allele_A_freq_by_gen.reserve(current_capacity *
								     fixed_parameters::factor_to_expand_vector_memory);
    }
	
  }

}

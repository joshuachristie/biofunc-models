#include <cmath>
#include <cassert>
#include <vector>
#include "conditional_existence_status.h"

namespace conditional_existence_status {
  /**
     @brief Compare for approximate equality of doubles
     @param[in] trait_freq Trait frequency
     @param[in] value Value to compare \p trait_freq against
     @param[in] tolerance A double that determines how closely \p trait_freq and \p value must be
     @return True if the absolute difference between \p trait_freq and \p value is less than \p tolerance; false otherwise
  */
  bool close_to_value(const double trait_freq, const double value, const double tolerance){
    if (std::abs(trait_freq - value) < tolerance){
      return true;
    } else {
      return false;
    }
  }
  /**
     @brief Gets the frequency of the A allele (for the haploid and diploid models)
     @param[in] trait_freq Trait frequency vector (A allele freq for haploid; AA and Aa freqs for diploid)
     @return Frequency of the A allele
  */
  double get_allele_A_freq(const std::vector<double> &trait_freq){
    double allele_A_frequency;
    if (trait_freq.size() == 1){ // haploid
      allele_A_frequency = trait_freq[0];
    } else if (trait_freq.size() == 2){ // diploid
      allele_A_frequency = trait_freq[0] + 0.5 * trait_freq[1];
    } else {
      allele_A_frequency = -1; // supress compiler warning
      assert(trait_freq.size() != 0 && "trait_freq.size() = 0, should be 1 (haploid) or 2 (diploid)");
      assert(trait_freq.size() > 2 && "trait_freq.size() > 2, should be 1 (haploid) or 2 (diploid)");
      exit(EXIT_FAILURE); // should never trigger but serious problem if program got through the two asserts
    }
    return allele_A_frequency;
  }
  
}

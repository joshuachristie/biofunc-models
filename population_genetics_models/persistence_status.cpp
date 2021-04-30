#include <cmath>
#include "persistence_status.h"

namespace persist_status {
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

}

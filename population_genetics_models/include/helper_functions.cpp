#include <cmath>
#include "helper_functions.h"

bool closeToValue(double allele_freq, double value, double tolerance){
  // returns true if allele_freq is approx equal to value
  if (std::abs(allele_freq - value) < tolerance){
    return 1;
  } else {
    return 0;
  }
}

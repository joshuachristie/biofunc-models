// provides closeToValue(), which tests whether for equality of a double

#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <cmath>

bool closeToValue(double allele_freq, double value, double tolerance){
  // returns true if allele_freq is approx equal to value
  if (std::abs(allele_freq - value) < tolerance){
    return 1;
  } else {
    return 0;
  }
}

#endif 

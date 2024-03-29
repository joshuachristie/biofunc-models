/**
   @file rng.h
   @brief Provides \p initialise_rng
*/

#ifndef RNG_H
#define RNG_H

#include <random>

namespace rng {
  
  std::mt19937 initialise_rng();
 
}
#endif 

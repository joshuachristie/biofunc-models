// provides initialiseRNG(), which sets a random seed using both
// random_device and the high precision clock

#ifndef RNG_H
#define RNG_H

#include <random>
#include <chrono>

std::mt19937 initialiseRNG(){
  std::mt19937 temp_rng(std::random_device{}());
  std::uniform_int_distribution<> adjust_seed(0, 50000);
  int factor_to_adjust_seed = adjust_seed(temp_rng);
  auto seed =
    (std::chrono::high_resolution_clock::now().time_since_epoch().count()) * factor_to_adjust_seed;
  std::mt19937 rng(seed);
  return rng;
}

#endif 

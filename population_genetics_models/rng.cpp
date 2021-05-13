#include <random>
#include <chrono>
#include "rng.h"

namespace rng {
  
  /**
     @brief Initialises a Mersenne Twister rng object using both \p std::random_device and \p std::chrono::high_resolution_clock
     @return An rng object
  */
  std::mt19937 initialise_rng(){
    std::mt19937 temp_rng(std::random_device{}());
    std::uniform_int_distribution<> adjust_seed(0, 50000);
    int factor_to_adjust_seed = adjust_seed(temp_rng);
    auto seed =
      (std::chrono::high_resolution_clock::now().time_since_epoch().count()) * factor_to_adjust_seed;
    std::mt19937 rng(seed);
    return rng;
  }

}

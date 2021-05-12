#ifndef FIXED_PARAMETERS_H
#define FIXED_PARAMETERS_H

namespace fixed_parameters {

  inline constexpr double tolerance = 0.000000000001;
  inline constexpr int number_replicates_QEF = 1000000;
  inline constexpr int number_replicates_LSTM = 1000;
  inline constexpr int max_generations_per_sim = 1000000;
  inline constexpr std::size_t factor_to_expand_vector_memory = 10;
  inline constexpr int reserve_memory_trait_freq = 20;
  
}

#endif

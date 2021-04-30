#ifndef PATHS_H
#define PATHS_H

#include <string_view>

namespace paths {

  inline constexpr std::string_view persistence_infinite_data_dir = "../data/persistence_probability/infinite_approximation/";
  inline constexpr std::string_view persistence_finite_data_dir = "../data/persistence_probability/finite_generations/";
  inline constexpr std::string_view trait_data_dir = "../data/raw_trait_data/";
  inline constexpr std::string_view error_file_directory = "../run/error_logs/";

}

#endif

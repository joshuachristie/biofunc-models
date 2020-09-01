#ifndef PATHS_H
#define PATHS_H

#include <string_view>

namespace paths {

  inline constexpr std::string_view persistence_data_directory = "../data/persistence_probability/";
  inline constexpr std::string_view error_file_directory = "../run/error_logs/";

}

#endif

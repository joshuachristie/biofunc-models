/**
   @file print_results.cpp
   @brief Methods to print model results to output files
*/

#include <sstream>
#include <fstream>
#include <string>
#include "print_results.h"
#include "io.h"
#include "path_parameters.h"

namespace print {
  /**
     @brief Prints persistence probability to file
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @param[in] persistence_probability Probability that the allele persists in population
  */
  void print_persistence_probability(int argc, char* argv[], const double persistence_probability){
    const std::string file_path =
      io::create_dir_and_get_file_path(argc, argv, paths::persistence_data_directory, ".csv", argv[1]);
    write_value_to_file(persistence_probability, file_path);
  }
  /**
     @brief Checks whether file is empty
     @param[in] infile File to check
     @return True if file is empty; false otherwise
  */
  bool is_empty(std::ifstream& infile){
    return infile.peek() == std::ifstream::traits_type::eof();
  }
  
}

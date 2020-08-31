/**
   @file print_results.cpp
   @brief Methods to print model results to output files
*/

#include <sstream>
#include <fstream>
#include <string>
#include "print_results.h"
#include "io.h"
#include "file_parameters.h"

namespace print {
  /**
     @brief Prints persistence probability to file
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @param[in] persistence_probability Probability that the allele persists in population
  */
  void print_persistence_probability(int argc, char* argv[], const double persistence_probability){
    const std::string dir_path = io::create_dir(file::data_directory, argv[1]);
    const std::string file_path = io::get_file_path(argc, argv, dir_path, ".csv");
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

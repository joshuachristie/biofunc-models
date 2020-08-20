/**
   @file print_results.cpp
   @brief Methods to print model results to output files
*/
#include <filesystem>
#include <sstream>
#include <fstream>
#include <string>
#include "print_results.h"

namespace print {
  /**
     @brief Prints persistence probability to file
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @param[in] persistence_probability Probability that the allele persists in population
  */
  void print_persistence_probability(int argc, char* argv[], const double persistence_probability){
    std::ostringstream filename = create_dir_and_get_filename(argc, argv);
    write_value_to_file(persistence_probability, filename);
  }
  /**
     @brief Checks whether file is empty
     @param[in] infile File to check
     @return True if file is empty; false otherwise
  */
  bool is_empty(std::ifstream& infile){
    return infile.peek() == std::ifstream::traits_type::eof();
  }
  /**
     @brief Creates directory for results and returns filename
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @return filename \p std::ostringstream with path to csv file
*/
  std::ostringstream create_dir_and_get_filename(int argc, char* argv[]){
    std::string parent_directory = "../data/";
    std::ostringstream dir_path;
    dir_path << parent_directory << argv[1] << "/";
    std::filesystem::create_directories(dir_path.str());
    // generate filename for csv
    std::ostringstream filename (dir_path.str(), std::ios_base::ate);
    for (int i = 2; i < argc - 1; i++){
      filename << argv[i] << "_";
    }
    filename << argv[argc -1] << ".csv";
    return filename;
  }

}

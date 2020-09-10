/**
   @file print_results.h
*/
#ifndef PRINT_RESULTS_H
#define PRINT_RESULTS_H

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdio>
#include <numeric>
#include "path_parameters.h"
#include "io.h"
#include "DataContainers.h"
#include "fixed_parameters.h"

/** Namespace for print methods*/
namespace print {

  bool is_empty(std::ifstream& infile);
  void write_to_file(const double value_to_write, const std::string &filename);
  void write_to_file(const std::vector<double> &vector_to_write, const std::string &filename);

  /**
     @brief Prints persistence probability (approximation as time -> infinity) to file
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @param[in] persistence_prob Probability that the allele persists in population (\p double)
  */
  template<class T>
  void print_persistence_probability(int argc, char* argv[], const T persistence_prob,
				     const std::string_view &path){
    const std::string file_path =
      io::create_dir_and_get_file_path(argc, argv, path, ".csv", argv[1]);
    write_to_file(persistence_prob, file_path);
  }

  void print_results(int argc, char* argv[], DataContainer &data);

}

#endif

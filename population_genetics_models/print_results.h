/**
   @file print_results.h
   @brief Methods to print model results to output files
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
#include <cassert>
#include "path_parameters.h"
#include "io.h"
#include "DataContainer.h"

/** Namespace for print methods*/
namespace print {

  bool is_empty(std::ifstream& infile);
  void write_to_file(const double value_to_write, const std::string &filename);
  void write_to_file(const std::vector<double> &vector_to_write, const std::string &filename);

  /**
     @brief Prints persistence probability (approximation as time -> infinity) to file
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @param[in] persistence_prob Probability that the trait persists in population (\p double)
  */
  template<class T>
  void print_object(int argc, char* argv[], const T object, const std::string_view &path){
    const std::string file_path =
      io::create_dir_and_get_file_path(argc, argv, path, ".csv", argv[1]);
    write_to_file(object, file_path);
  }

  template<class P>
  void print_results(int argc, char* argv[], DataContainer &data, P &parameters){
    const double persistence_probability = data.get_persistence_infinite_approx();
    print_object(argc, argv, persistence_probability, paths::persistence_infinite_data_dir);
    if (parameters.shared.number_gens_to_output_pp != 0){
      assert(!data._simulation_data[0]._persistence_by_gen.empty());
      const std::vector<double> persistence_prob_by_gen = data.get_persistence_by_gen();
      print_object(argc, argv, persistence_prob_by_gen, paths::persistence_finite_data_dir);
    } else {
      assert(data._simulation_data[0]._persistence_by_gen.empty());
    }
    if (parameters.shared.print_trait_raw_data){
      assert(!data._simulation_data[0]._trait_freq_by_gen.empty());
      for (int i = 0; i < parameters.fixed.number_replicates; i++){
	const std::vector<double> trait_freqs = data.get_trait_freqs(i);
	print_object(argc, argv, trait_freqs, paths::trait_data_dir);
      }
    } else {
      assert(data._simulation_data[0]._trait_freq_by_gen.empty());
    }
  }

}

#endif

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

  template<class P>
  const std::vector<double> process_persistence_probability_finite(const P &parameters,
								   const std::vector<bool> &final_A_freqs){
    std::vector<double> persistence_probability(parameters.shared.number_gens_to_output_pp + 1, 0.0);
    for (int i = 0; i < parameters.fixed.number_replicates; i++){
      for (int j = 0; j < (parameters.shared.number_gens_to_output_pp + 1); j++){
	persistence_probability[j] +=
	  (static_cast<double>(final_A_freqs[j + i * (parameters.shared.number_gens_to_output_pp + 1)])
	   / static_cast<double>(parameters.fixed.number_replicates));
      }
    }
    return persistence_probability;
  }
  
  template<class P>
  const double process_persistence_probability_infinite(const P &parameters,
							const std::vector<bool> &final_A_freqs){
    return std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) /
      static_cast<double>(parameters.fixed.number_replicates);
  }
  
  template<class P>
  void print_results(int argc, char* argv[], const P &parameters, const std::vector<bool> &final_A_freqs){
    if (final_A_freqs.size() == static_cast<std::size_t>(parameters.fixed.number_replicates)){

      const double persistence_probability = process_persistence_probability_infinite(parameters, final_A_freqs);

      print_persistence_probability(argc, argv, persistence_probability, paths::persistence_infinite_data_dir);
    } else if (final_A_freqs.size() > static_cast<std::size_t>(parameters.fixed.number_replicates)) {

      const std::vector<double> persistence_probability = process_persistence_probability_finite(parameters, final_A_freqs);

      print_persistence_probability(argc, argv, persistence_probability, paths::persistence_finite_data_dir);
      print_persistence_probability(argc, argv, persistence_probability[persistence_probability.size() - 1],
				    paths::persistence_infinite_data_dir);
    } else { // should never execute; if it does then there's a major error
      try {
	throw std::runtime_error("final_A_freqs.size() < paramters.fixed.number_replicates...exiting program");
      } catch (const std::exception &e){
	const std::string error_file_path =
	  io::create_dir_and_get_file_path(argc, argv, paths::error_file_directory, "_error.txt");
	std::ofstream error_file(error_file_path, std::ostream::app);
	error_file << "exception: " << e.what() << "\n";
	exit(EXIT_FAILURE);
      }
      
    }
  }

}

#endif

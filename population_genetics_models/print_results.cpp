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
     @brief Method (overloaded) that writes a double to file
     @param[in] value_to_write Value to write to file
     @param[in] filename Name of file to write to
     @return Nothing (but writes \p value_to_write to \p filename)
  */
  void write_to_file(const double value_to_write, const std::string &filename){
    // check whether file is empty
    std::ifstream infile (filename);
    if (is_empty(infile)){ // if so, no comma
      infile.close();
      std::ofstream outfile (filename, std::ofstream::app);
      outfile << value_to_write;
    } else { // if not empty, add comma before value
      infile.close();
      std::ofstream outfile (filename, std::ofstream::app);
      outfile << "," << value_to_write;
    }
  }
  /**
     @brief Method (overloaded) that writes a std::vector<double> to file
     @param[in] vector_to_write Vector to write to file
     @param[in] filename Name of file to write to
     @return Nothing (but writes \p vector_to_write to \p filename)
  */
  void write_to_file(const std::vector<double> &vector_to_write, const std::string &filename){
    std::ofstream outfile (filename, std::ofstream::app);
    outfile << vector_to_write[0];
    for (std::size_t i = 1; i < vector_to_write.size(); i++){
      outfile << "," << vector_to_write[i];
    }
    outfile << "\n";
  }
  
  /**
     @brief Checks whether file is empty
     @param[in] infile File to check
     @return True if file is empty; false otherwise
  */
  bool is_empty(std::ifstream& infile){
    return infile.peek() == std::ifstream::traits_type::eof();
  }

  void print_results(int argc, char* argv[], const DataContainer<DataPersistenceInfinite, number_replicates> &data){
    
    const double persistence_probability = data.get_persistence_infinite_approx();
    print_persistence_probability(argc, argv, persistence_probability, paths::persistence_infinite_data_dir);
    
  }

  void print_results(int argc, char* argv[], const DataContainer<DataPersistenceByGen, number_replicates> &data){
    
    const double persistence_prob_infinite = data.get_persistence_infinite_approx();
    print_persistence_probability(argc, argv, persistence_prob_infinite, paths::persistence_infinite_data_dir);
    
    const std::vector<double> persistence_prob_by_gen = data.get_persistence_by_gen();
    print_persistence_probability(argc, argv, persistence_prob_by_gen, paths::persistence_finite_data_dir);
    
  }
  
  void print_results(int argc, char* argv[], const DataContainer<DataAlleleFreqAndPGB, number_replicates> &data){
    
    const double persistence_prob_infinite = data.get_persistence_infinite_approx();
    print_persistence_probability(argc, argv, persistence_prob_infinite, paths::persistence_infinite_data_dir);
    
    const std::vector<double> persistence_prob_by_gen = data.get_persistence_by_gen();
    print_persistence_probability(argc, argv, persistence_prob_by_gen, paths::persistence_finite_data_dir);
    
    for (int i = 0; i < number_replicates; i++){
      const std::vector<double> allele_A_freqs = data.get_allele_A_freqs(i);
      print_persistence_probability(argc, argv, allele_A_freqs, paths::allele_A_data_dir);
    }
    
  }

  void print_results(int argc, char* argv[], const DataContainer<DataAlleleFreq, number_replicates> &data){
      
    const double persistence_prob_infinite = data.get_persistence_infinite_approx();
    print_persistence_probability(argc, argv, persistence_prob_infinite, paths::persistence_infinite_data_dir);
    
    for (int i = 0; i < number_replicates; i++){
      const std::vector<double> allele_A_freqs = data.get_allele_A_freqs(i);
      print_persistence_probability(argc, argv, allele_A_freqs, paths::allele_A_data_dir);
    }
    
  }
  
}

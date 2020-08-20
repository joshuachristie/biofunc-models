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
  */
  void print_persistence_probability(int argc, char* argv[], const double persistence_probability){
    // make directory
    // std::string parent_directory = "./data/";
    // std::ostringstream dir_path;
    // dir_path << parent_directory << argv[1] << "/";
    // std::filesystem::create_directories(dir_path.str());
    // // generate filename for csv
    // std::ostringstream filename (dir_path.str(), std::ios_base::ate);
    // for (int i = 2; i < argc - 1; i++){
    //   filename << argv[i] << "_";
    // }
    // filename << argv[argc -1] << ".csv";
    std::ostringstream filename = create_dir_and_get_filename(argc, argv);
    write_value_to_file(persistence_probability, filename);


    
    // check whether file is empty
    // std::ifstream infile (filename.str());
    // if (is_empty(infile)){ // if so, no comma
    //   infile.close();
    //   std::ofstream outfile (filename.str(), std::ofstream::app);
    //   outfile << allele_A_freq;
    // } else { // if not empty, add comma before value
    //   infile.close();
    //   std::ofstream outfile (filename.str(), std::ofstream::app);
    //   outfile << "," << allele_A_freq;
    // }
  }

  
  /**
     @brief Checks whether file is empty
     @param[in] infile File to check
     @return True if file is empty; false otherwise
  */
  bool is_empty(std::ifstream& infile){
    return infile.peek() == std::ifstream::traits_type::eof();
  }

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

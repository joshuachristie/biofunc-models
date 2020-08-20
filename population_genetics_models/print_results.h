#ifndef PRINT_RESULTS_H
#define PRINT_RESULTS_H

#include <fstream>
#include <sstream>

namespace print {
  
  void print_persistence_probability(int argc, char* argv[], const double allele_A_freq);
  bool is_empty(std::ifstream& infile);
  std::ostringstream create_dir_and_get_filename(int argc, char* argv[]);

  template <class T>
  void write_value_to_file(const T value_to_write, std::ostringstream &filename){
    // check whether file is empty
    std::ifstream infile (filename.str());
    if (is_empty(infile)){ // if so, no comma
      infile.close();
      std::ofstream outfile (filename.str(), std::ofstream::app);
      outfile << value_to_write;
    } else { // if not empty, add comma before value
      infile.close();
      std::ofstream outfile (filename.str(), std::ofstream::app);
      outfile << "," << value_to_write;
    }
  }

}

#endif

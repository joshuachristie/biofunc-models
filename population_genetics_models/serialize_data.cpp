#include <fstream>
#include "serialize_data.h"
#include "include/example.pb.h"
#include "path_parameters.h"
#include "io.h"

namespace serialize {
  
  void data(data::Example& example, int argc, char* argv[]){
    std::string filename = io::setup_dir_and_file(argc, argv, paths::QEF_directory);
    std::fstream output(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    example.SerializeToOstream(&output);
  }

  void data(data::SequenceExample& seq_example, int argc, char* argv[]){
    std::string filename = io::setup_dir_and_file(argc, argv, paths::LSTM_directory);
    std::fstream output(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    seq_example.SerializeToOstream(&output);
  }
  
}

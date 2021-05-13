#include <functional>
#include <fstream>
#include "model_specification.h"
#include "HSE.h"
#include "HTE.h"
#include "DSE.h"
#include "HTEOE.h"
#include "io.h"
#include "path_parameters.h"

/**
   @brief Constructs a \p std::map so that models can be called via a \p std::string argument
   @return \p std::map with a \p std::string as "key" and a function that runs the model as "value"
*/
model_map get_model_map(){
  model_map map {
    {"HSE", HSE::run_model},
    {"HTE", HTE::run_model},
    {"DSE", DSE::run_model},
    {"HTEOE", HTEOE::run_model}
  };
  return map;
}

/**
   @brief Specifies and runs one of the four models (HSE, DSE, HTE, HTEOE)
   @param[in] argc Number of command line arguments
   @param[in] argv Array of command line arguments
   @return Nothing (but runs a specified model and combination of parameter values)
**/
void specify_and_run_model(int argc, char* argv[]){
  model_map map = get_model_map(); // hashmap/dict of available models
  try {
    map[ argv[1] ](argc, argv); // specify and run model
  }
  catch (const std::bad_function_call &e){
    const std::string error_file_path =
      io::setup_dir_and_file(argc, argv, paths::error_file_directory, "_error.txt");
    std::ofstream error_file(error_file_path, std::ofstream::app);
    error_file << "exception: " << e.what() << " - is the first cmdline argument the model identifier (e.g. HSE)?" << "\n";
  }
  
}

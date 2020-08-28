#include <stdexcept>
#include <iostream>
#include <sstream>
#include "exceptions.h"
#include "print_results.h"

const int check_parameter_value_compatibility(const int number_reinvasions, int argc, char* argv[],
					      const int output_pp_index){
  try {
    if (number_reinvasions > 0 && atoi(argv[output_pp_index]) > 0){
      // throw parameter_value_exception;
      std::ostringstream error_msg;
      error_msg << "number_reinvasions = " << number_reinvasions << " and " << \
	"number_generations_to_output_pp = " << argv[output_pp_index] << \
	"; if number_reinvasions > 0, number_generations_to_output_pp should = 0; " << \
	"setting number_generations_to_output_pp = 0";
      *argv[output_pp_index] = '0'; // set commandline parameter to 0
      throw std::runtime_error(error_msg.str());
    }
  }
  catch (const std::exception& e){
    std::cout << "exception: " << e.what() << "; renamed simulation name is: " << get_simulation_name(argc, argv).str() << std::endl;
  }
  return atoi(argv[output_pp_index]);
}

std::ostringstream get_simulation_name(int argc, char* argv[]){
  std::ostringstream simname;
  for (int i = 1; i < argc - 1; i++){
    simname << argv[i] << "_";
  }
  simname << argv[argc - 1];
  return simname;
}

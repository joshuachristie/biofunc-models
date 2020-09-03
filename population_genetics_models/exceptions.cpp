#include <stdexcept>
#include <istream>
#include <sstream>
#include <string>
#include "exceptions.h"
#include "print_results.h"
#include "io.h"
#include "path_parameters.h"

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
    const std::string error_file_path =
      io::create_dir_and_get_file_path(argc, argv, paths::error_file_directory, "_error.txt");
    std::ofstream error_file(error_file_path, std::ostream::app);
    error_file << "exception: " << e.what() << "; renamed simulation name is: " << \
      io::parameter_values_to_string(argc, argv) << "\n";
  }
  return atoi(argv[output_pp_index]);
}

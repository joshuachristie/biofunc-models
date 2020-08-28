#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <sstream>

const int check_parameter_value_compatibility(const int number_reinvasions, int argc, char* argv[], const int output_pp_index);

std::ostringstream get_simulation_name(int argc, char* argv[]);

#endif

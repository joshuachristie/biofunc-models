/**
   @file io.h
   @brief Functions dealing with input/output operations
**/
#ifndef IO_H
#define IO_H

#include <string_view>

/** Namespace for input/output functions **/
namespace io {
  
  std::string parameter_values_to_string(int argc, char* argv[]);
  
  const std::string create_dir(const std::string_view &parent_dir, const std::string &dir = "");
  
  const std::string get_file_path(int argc, char* argv[], const std::string &dir_path,
				  const std::string &extension);

  const std::string setup_dir_and_file(int argc, char* argv[], const std::string_view &parent_dir,
				       const std::string &extension = "", const std::string &dir = "");

}

#endif

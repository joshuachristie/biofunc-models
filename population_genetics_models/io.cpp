#include <string>
#include <sstream>
#include <filesystem>
#include "io.h"

namespace io {

  /**
     @brief Converts parameter values to a string object
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @return stringname String containing the simulation's parameter values
*/
  std::string parameter_values_to_string(int argc, char* argv[]){
    std::ostringstream stringname;
    for (int i = 1; i < argc - 1; i++){
      stringname << argv[i] << "_";
    }
    stringname << argv[argc - 1];
    return stringname.str();
  }
  /**
     @brief Creates directory
     @param[in] parent_dir String pointing to parent directory
     @param[in] dir String that names the directory (optional: default is empty string)
     @return dir_path String of path to directory
*/
  const std::string create_dir(const std::string_view &parent_dir, const std::string &dir){
    std::ostringstream dir_path;
    dir_path << parent_dir << dir << "/";
    std::filesystem::create_directories(dir_path.str());
    return dir_path.str();
  }
    /**
     @brief Produces string with path to file
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @param[in] dir_path String containing path to directory
     @param[in] extension File extension (e.g. ".csv" or ".txt")
     @return filename \p String with path to file
*/
  const std::string get_file_path(int argc, char* argv[], const std::string &dir_path,
				  const std::string &extension){
    std::ostringstream filename (dir_path, std::ios_base::ate);
    filename << io::parameter_values_to_string(argc, argv) << extension;
    return filename.str();
  }

  /**
*/
  const std::string create_dir_and_get_file_path(int argc, char* argv[], const std::string_view &parent_dir,
						 const std::string &extension, const std::string &dir){
    const std::string dir_path = create_dir(parent_dir, dir);
    const std::string file_path = get_file_path(argc, argv, dir_path, extension);
    return file_path;
  }
}

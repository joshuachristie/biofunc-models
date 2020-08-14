#ifndef MODEL_SPEC_H
#define MODEL_SPEC_H

#include <map>
#include <string>
#include <functional>

/** typedef for \p std::map that is used to choose which model to run */
typedef std::map<std::string, std::function<void(int, char*[])>> model_map;

model_map get_model_map();
void specify_and_run_model(int argc, char* argv[]);

#endif

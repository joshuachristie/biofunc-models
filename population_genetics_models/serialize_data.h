#ifndef SERIALISE_DATA_H
#define SERIALISE_DATA_H

#include "include/example.pb.h"

namespace serialize {
  
  void data(tensorflow::Example& example, int argc, char* argv[]);
  void data(tensorflow::SequenceExample& seq_example, int argc, char* argv[]);

}

#endif

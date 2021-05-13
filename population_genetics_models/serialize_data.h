#ifndef SERIALISE_DATA_H
#define SERIALISE_DATA_H

#include "include/example.pb.h"

namespace serialize {
  
  void data(data::Example& example, int argc, char* argv[]);
  void data(data::SequenceExample& seq_example, int argc, char* argv[]);

}

#endif

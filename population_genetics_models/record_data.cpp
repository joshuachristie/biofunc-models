#include <vector>
#include "record_data.h"
#include "include/example.pb.h"

namespace record_data {

  void raw_trait_freq(tensorflow::FloatList* raw_trait_freq, const std::vector<double> &trait_freq){
    // single value per gen for HSE (P(A); two values per gen for DSE (P(AA) and P(Aa))
    for (std::size_t i = 0; i < trait_freq.size(); i++){
      raw_trait_freq->add_value(trait_freq[i]);
    }
  }

}

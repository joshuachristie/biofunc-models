#ifndef RECORD_DATA_H
#define RECORD_DATA_H

#include <string>
#include <vector>
#include "include/example.pb.h"
#include "conditional_existence_status.h"

namespace record {

  void raw_trait_freq(data::FloatList* raw_trait_freq, const std::vector<double> &trait_freq);

  template <class P>
  void generation_trait_extinction(data::Int64List* gen_extinct, const std::vector<double> &trait_freq,
				   const P &params, const int gen){
    if (conditional_existence_status::trait_extinct(trait_freq, params)){
      gen_extinct->add_value(gen); // extinct trait, record generation of extinction
    } else if (!conditional_existence_status::trait_extinct(trait_freq, params)){
      gen_extinct->add_value(params.fixed.max_generations_per_sim); // trait persists, denote by max gen
    } else {
      exit(EXIT_FAILURE);
    }
  }

  template <class P>
  void number_reinvasions_before_extinction(data::Int64List* reinvasion_number, const std::vector<double>
					    &trait_freq, const P &params, const int reinvasions){
    // -1 if not looking at reinvasions
    // if the trait is extinct, then the value of `reinvasions` is correct
    // if the trait persists, however, we need to add 1 (otherwise I can't distinguish between traits
    // that go extinct in the final reinvasion and those that survive all attempts)
    if (conditional_existence_status::trait_extinct(trait_freq, params)){
      reinvasion_number->add_value(reinvasions);
    } else if (!conditional_existence_status::trait_extinct(trait_freq, params)){ 
      reinvasion_number->add_value(reinvasions + 1);
    } else {
      exit(EXIT_FAILURE);
    }
  }

}

#endif

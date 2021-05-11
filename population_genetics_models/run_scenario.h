#ifndef RUN_SCENARIO_H
#define RUN_SCENARIO_H

#include <iostream>
#include <fstream>
#include "include/example.pb.h"
#include "conditional_existence_probability.h"

namespace run_scenario {

  template <class P, class F>
  void QEF(const P &params, std::mt19937 &rng, const std::vector<double> &fitnesses, F calculate_trait_freqs){

    data::Example example = data::Example();
    data::Features* features = example.mutable_features();
    google::protobuf::Map<std::string, data::Feature>* feature_map = features->mutable_feature();

    std::string key_gen = "generation_of_extinction";
    data::Feature generation_of_extinction = data::Feature();
    data::Int64List* gen_extinct = generation_of_extinction.mutable_int64_list();

    std::string key_reinvasion = "number_reinvasions";
    data::Feature number_reinvasions = data::Feature();
    data::Int64List* reinvasion_number = number_reinvasions.mutable_int64_list();

    calculate_conditional_existence_probability(params, rng, fitnesses, calculate_trait_freqs,
						gen_extinct, reinvasion_number);
    
    (*feature_map)[key_gen] = generation_of_extinction;
    (*feature_map)[key_reinvasion] = number_reinvasions;

    // below will move into its own function (in the io.cpp/h files if I keep them)
    std::fstream output("test", std::ios::out | std::ios::trunc | std::ios::binary);
    example.SerializeToOstream(&output);
  }

  // template <class P, class F>
  // void LSTM(const P &params, std::mt19937 &rng, const std::vector<double> &fitnesses, F calculate_trait_freqs){
  //   data::SequenceExample seq_example = data::SequenceExample(); // qef will only use Example
  //   data::Features* features = seq_example.mutable_context();
  //   google::protobuf::Map<std::string, data::Feature>* feature_map = features->mutable_feature();
  //   data::Feature myfeature = data::Feature();
  //   data::Int64List* intlist = myfeature.mutable_int64_list();
    
  //   calculate_conditional_existence_probability(params, rng, fitnesses, calculate_trait_freqs, intlist);
    
  //   (*feature_map)["trait_presence"] = myfeature;

    
  //   std::fstream output("test", std::ios::out | std::ios::trunc | std::ios::binary);
  //   seq_example.SerializeToOstream(&output);

  // }

  
}

#endif

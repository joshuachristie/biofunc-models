#ifndef RUN_SCENARIO_H
#define RUN_SCENARIO_H

#include <string>
#include <vector>
#include "include/example.pb.h"
#include "conditional_existence_probability.h"
#include "record_context.h"
#include "serialize_data.h"

namespace run_scenario {

  template <class P, class F>
  void QEF(const P &params, std::mt19937 &rng, const std::vector<double> &fitnesses,
	   F calculate_trait_freqs, char* argv[], int argc){

    tensorflow::Example example = tensorflow::Example();
    tensorflow::Features* features = example.mutable_features();
    google::protobuf::Map<std::string, tensorflow::Feature>* feature_map = features->mutable_feature();

    std::string key_gen = "generation_of_extinction";
    tensorflow::Feature generation_of_extinction = tensorflow::Feature();
    tensorflow::Int64List* gen_extinct = generation_of_extinction.mutable_int64_list();
    std::string key_reinvasion = "number_reinvasions";
    tensorflow::Feature number_reinvasions = tensorflow::Feature();
    tensorflow::Int64List* reinvasion_number = number_reinvasions.mutable_int64_list();

    conditional_existence_probability::calculate(params, rng, fitnesses, calculate_trait_freqs,
						 gen_extinct, reinvasion_number);
    
    (*feature_map)[key_gen] = generation_of_extinction;
    (*feature_map)[key_reinvasion] = number_reinvasions;
    record_context::add_parameters_to_protobuf(feature_map, params, argv); // metadata, parameter values, etc.
    serialize::data(example, argc, argv);
  }

  template <class P, class F>
  void LSTM(const P &params, std::mt19937 &rng, const std::vector<double> &fitnesses,
	    F calculate_trait_freqs, char* argv[], int argc){
    tensorflow::SequenceExample seq_example = tensorflow::SequenceExample();
    // generation of extinction
    tensorflow::Features* features = seq_example.mutable_context();
    google::protobuf::Map<std::string, tensorflow::Feature>* feature_map = features->mutable_feature();
    
    std::string key_gen = "generation_of_extinction";
    tensorflow::Feature generation_of_extinction = tensorflow::Feature();
    tensorflow::Int64List* gen_extinct = generation_of_extinction.mutable_int64_list();
    
    // raw trait data
    tensorflow::FeatureLists* featurelists = seq_example.mutable_feature_lists();
  
    google::protobuf::Map<std::string, tensorflow::FeatureList>* featurelist_map =
      featurelists->mutable_feature_list();
    std::string key_freq = "raw_trait_frequencies";
    tensorflow::FeatureList featurelist = tensorflow::FeatureList();

    conditional_existence_probability::calculate(params, rng, fitnesses, calculate_trait_freqs,
						 gen_extinct, featurelist);
						
    (*feature_map)[key_gen] = generation_of_extinction;
    (*featurelist_map)[key_freq] = featurelist;
    record_context::add_parameters_to_protobuf(feature_map, params, argv); // metadata, parameter values, etc.
    serialize::data(seq_example, argc, argv);
  }

}

#endif

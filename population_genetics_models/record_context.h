#ifndef RECORD_CONTEXT_H
#define RECORD_CONTEXT_H

#include "include/example.pb.h"
#include "Parameters.h"

namespace record_context {

  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, tensorflow::Feature>* map,
					   parameters::HSE_Model_Parameters params);
  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, tensorflow::Feature>* map,
					   parameters::DSE_Model_Parameters params);
  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, tensorflow::Feature>* map,
					   parameters::HTE_Model_Parameters params);
  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, tensorflow::Feature>* map,
					   parameters::HTEOE_Model_Parameters params);

  template<class P>
  void add_shared_parameters_to_protobuf(google::protobuf::Map<std::string, tensorflow::Feature>* map,
					 P params, char* argv[]){
    tensorflow::Feature model = tensorflow::Feature();
    tensorflow::BytesList* model_name = model.mutable_bytes_list();
    model_name->add_value(std::string(argv[1]));
    (*map)["model"] = model;

    tensorflow::Feature scenario = tensorflow::Feature();
    tensorflow::BytesList* scenario_name = scenario.mutable_bytes_list();
    scenario_name->add_value(std::string(argv[2])); // QEF (all) or LSTM (HSE/DSE only)
    (*map)["scenario"] = scenario;

    tensorflow::Feature pop = tensorflow::Feature();
    tensorflow::Int64List* pop_size = pop.mutable_int64_list();
    pop_size->add_value(params.shared.population_size);
    (*map)["population_size"] = pop;

    tensorflow::Feature reinvasions = tensorflow::Feature();
    tensorflow::Int64List* num_reinvasions = reinvasions.mutable_int64_list();
    num_reinvasions->add_value(params.shared.number_reinvasions);
    (*map)["number_reinvasions"] = reinvasions;

    tensorflow::Feature trait = tensorflow::Feature();
    tensorflow::Int64List* trait_info = trait.mutable_int64_list();
    trait_info->add_value(params.shared.trait_info[0]);
    trait_info->add_value(params.shared.trait_info[1]);
    (*map)["trait_info"] = trait;
  }

  template<class P>
  void add_parameters_to_protobuf(google::protobuf::Map<std::string, tensorflow::Feature>* map,
				  P params, char* argv[]){
    add_specific_parameters_to_protobuf(map, params);
    add_shared_parameters_to_protobuf(map, params, argv);
  }
  
}

#endif

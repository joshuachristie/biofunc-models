#ifndef RECORD_CONTEXT_H
#define RECORD_CONTEXT_H

#include "include/example.pb.h"
#include "Parameters.h"

namespace record_context {

  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, data::Feature>* map,
					   HSE_Model_Parameters params, char* argv[]);
  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, data::Feature>* map,
					   DSE_Model_Parameters params, char* argv[]);
  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, data::Feature>* map,
					   HTE_Model_Parameters params, char* argv[]);
  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, data::Feature>* map,
					   HTEOE_Model_Parameters params, char* argv[]);

  template<class P>
  void add_shared_parameters_to_protobuf(google::protobuf::Map<std::string, data::Feature>* map,
					 P params, char* argv[]){
    data::Feature model = data::Feature();
    data::BytesList* model_name = model.mutable_bytes_list();
    model_name->add_value(std::string(argv[1]));
    (*map)["model"] = model;

    data::Feature pop = data::Feature();
    data::Int64List* pop_size = pop.mutable_int64_list();
    pop_size->add_value(params.shared.population_size);
    (*map)["population_size"] = pop;

    data::Feature reinvasions = data::Feature();
    data::Int64List* num_reinvasions = reinvasions.mutable_int64_list();
    num_reinvasions->add_value(params.shared.number_reinvasions);
    (*map)["number_reinvasions"] = reinvasions;

    data::Feature trait = data::Feature();
    data::Int64List* trait_info = trait.mutable_int64_list();
    trait_info->add_value(params.shared.trait_info[0]);
    trait_info->add_value(params.shared.trait_info[1]);
    (*map)["trait_info"] = trait;
  }

  template<class P>
  void add_parameters_to_protobuf(google::protobuf::Map<std::string, data::Feature>* map,
				  P params, char* argv[]){
    add_specific_parameters_to_protobuf(map, params, argv);
    add_shared_parameters_to_protobuf(map, params, argv);
  }
  
}

#endif

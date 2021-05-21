#include <string>
#include "record_context.h"
#include "include/example.pb.h"
#include "Parameters.h"

namespace record_context {
  
  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, tensorflow::Feature>* map,
					   parameters::HSE_Model_Parameters params){
    tensorflow::Feature selection = tensorflow::Feature();
    tensorflow::FloatList* selection_coefficient = selection.mutable_float_list();
    selection_coefficient->add_value(params.model.selection_coefficient);
    (*map)["selection_coefficient"] = selection;
  }
  
  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, tensorflow::Feature>* map,
					   parameters::DSE_Model_Parameters params){
    tensorflow::Feature selection_homo = tensorflow::Feature();
    tensorflow::FloatList* selection_homozygote = selection_homo.mutable_float_list();
    selection_homozygote->add_value(params.model.selection_coefficient_homozygote);
    (*map)["selection_coefficient_homozygote"] = selection_homo;

    tensorflow::Feature selection_hetero = tensorflow::Feature();
    tensorflow::FloatList* selection_heterozygote = selection_hetero.mutable_float_list();
    selection_heterozygote->add_value(params.model.selection_coefficient_heterozygote);
    (*map)["selection_coefficient_heterozygote"] = selection_hetero;
  }
  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, tensorflow::Feature>* map,
					   parameters::HTE_Model_Parameters params){
    tensorflow::Feature selection_Ae1 = tensorflow::Feature();
    tensorflow::FloatList* selection_A_env_1 = selection_Ae1.mutable_float_list();
    selection_A_env_1->add_value(params.model.selection_coefficient_A_env_1);
    (*map)["selection_coefficient_A_env_1"] = selection_Ae1;

    tensorflow::Feature selection_A2 = tensorflow::Feature();
    tensorflow::FloatList* selection_A_env_2 = selection_A2.mutable_float_list();
    selection_A_env_2->add_value(params.model.selection_coefficient_A_env_2);
    (*map)["selection_coefficient_A_env_2"] = selection_A2;

    tensorflow::Feature selection_ae1 = tensorflow::Feature();
    tensorflow::FloatList* selection_a_env_1 = selection_ae1.mutable_float_list();
    selection_a_env_1->add_value(params.model.selection_coefficient_a_env_1);
    (*map)["selection_coefficient_a_env_1"] = selection_ae1;

    tensorflow::Feature selection_a2 = tensorflow::Feature();
    tensorflow::FloatList* selection_a_env_2 = selection_a2.mutable_float_list();
    selection_a_env_2->add_value(params.model.selection_coefficient_a_env_2);
    (*map)["selection_coefficient_a_env_2"] = selection_a2;

    tensorflow::Feature gen_e1 = tensorflow::Feature();
    tensorflow::Int64List* gen_env_1 = gen_e1.mutable_int64_list();
    gen_env_1->add_value(params.model.gen_env_1);
    (*map)["gen_env_1"] = gen_e1;
  }

  void add_specific_parameters_to_protobuf(google::protobuf::Map<std::string, tensorflow::Feature>* map,
					   parameters::HTEOE_Model_Parameters params){
    tensorflow::Feature selection_A1 = tensorflow::Feature();
    tensorflow::FloatList* selection_coefficient_A1 = selection_A1.mutable_float_list();
    selection_coefficient_A1->add_value(params.model.selection_coefficient_A1);
    (*map)["selection_coefficient_A1"] = selection_A1;

    tensorflow::Feature selection_A2 = tensorflow::Feature();
    tensorflow::FloatList* selection_coefficient_A2 = selection_A2.mutable_float_list();
    selection_coefficient_A2->add_value(params.model.selection_coefficient_A2);
    (*map)["selection_coefficient_A2"] = selection_A2;

    tensorflow::Feature selection_a1 = tensorflow::Feature();
    tensorflow::FloatList* selection_coefficient_a1 = selection_a1.mutable_float_list();
    selection_coefficient_a1->add_value(params.model.selection_coefficient_a1);
    (*map)["selection_coefficient_a1"] = selection_a1;

    tensorflow::Feature selection_a2 = tensorflow::Feature();
    tensorflow::FloatList* selection_coefficient_a2 = selection_a2.mutable_float_list();
    selection_coefficient_a2->add_value(params.model.selection_coefficient_a2);
    (*map)["selection_coefficient_a2"] = selection_a2;
  }

}

/**
@file
@brief Definitions for parameter structs
*/

#ifndef PARAMETERS_H
#define PARAMETERS_H

/**
@brief Struct containing shared parameters (common to all models)
*/
struct Shared_Parameters {
  const int population_size;
  const int number_replicates;
};

/**
@brief Struct for parameters of the Haploid Single Environment model
*/
struct HSE_Model_Parameters {
  // add struct containing shared parameters
  Shared_Parameters shared;
  // add struct for HSE-specific parameters
  struct HSE_Specific_Parameters {
    const int number_generations;
    const double selection_coefficient;
  } model;
};

/**
@brief Struct for parameters of the Diploid Single Environment model
*/
struct DSE_Model_Parameters {
  // add struct containing shared parameters
  Shared_Parameters shared;
  // add struct for DSE-specific parameters
  struct DSE_Specific_Parameters {
    const int number_generations;
    const double selection_coefficient_homozygote;
    const double selection_coefficient_heterozygote;
  } model;
};
/**
@brief Struct for parameters of the Haploid Two Environments model
*/
struct HTE_Model_Parameters {
  // add struct containing shared parameters
  Shared_Parameters shared;
  // add struct for HTE-specific parameters
  struct HTE_Specific_Parameters {
    const double selection_coefficient_A_env_1;
    const double selection_coefficient_A_env_2;
    const double selection_coefficient_a_env_1;
    const double selection_coefficient_a_env_2;
    const int gen_env_1;
    const int gen_env_2;
  } model;
};

/**
@brief Struct for parameters of the Haploid Two Effects One Environment model
*/
struct HTEOE_Model_Parameters {
  // add struct containing shared parameters
  Shared_Parameters shared;
  // add struct for HTEOE-specific parameters
  struct HTEOE_Specific_Parameters {
    const double selection_coefficient_A1;
    const double selection_coefficient_A2;
    const double selection_coefficient_a1;
    const double selection_coefficient_a2;
    const int number_generations;
  } model;
};

#endif

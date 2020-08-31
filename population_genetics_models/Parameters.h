/**
   @file
   @brief Definitions for parameter structs
*/

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "fixed_parameters.h"

/**
   @brief Struct containing shared parameters (common to all models)
*/
struct Shared_Parameters {
  const int population_size;
  const double initial_A_freq;
  const int number_reinvasions;
  const int number_generations_to_output_pp;
};
/**
   @brief Struct containing constant parameter values (fixed for all models)
*/
struct Fixed_Parameters {
  inline static const double tolerance = fixed_parameters::tolerance;
  inline static const int number_replicates = fixed_parameters::number_replicates;
  inline static const int max_generations_per_sim = fixed_parameters::max_generations_per_sim;
};
/**
   @brief Struct for parameters of the Haploid Single Environment model
*/
struct HSE_Model_Parameters {
  // add struct containing shared parameters
  Shared_Parameters shared;
  // add struct for HSE-specific parameters
  struct HSE_Specific_Parameters {
    const double selection_coefficient;
  } model;
  // add struct containing fixed parameter values
  Fixed_Parameters fixed;
};

/**
   @brief Struct for parameters of the Diploid Single Environment model
*/
struct DSE_Model_Parameters {
  // add struct containing shared parameters
  Shared_Parameters shared;
  // add struct for DSE-specific parameters
  struct DSE_Specific_Parameters {
    const double selection_coefficient_homozygote;
    const double selection_coefficient_heterozygote;
  } model;
  // add struct containing fixed parameter values
  Fixed_Parameters fixed;
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
  } model;
  // add struct containing fixed parameter values
  Fixed_Parameters fixed;
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
  } model;
  // add struct containing fixed parameter values
  Fixed_Parameters fixed;
};

#endif

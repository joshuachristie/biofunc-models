/**
   @file
   @brief Stucts that store parameter values for the model
*/

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <vector>
#include "fixed_parameters.h"

/**
   @brief Struct containing shared parameters (i.e. parameters utilised by all of the population genetics models)
*/
struct Shared_Parameters {
  const int population_size; /**< Number of individuals in the population */
  const double initial_trait_freq; /**< Initial frequency of the trait  */
  const int number_reinvasions; /**< Number of reinvasion attempts that the trait faces */
  /** Sets the number of generations for which to calculate and print time-dependent persistence_probability */
  const int number_gens_to_output_pp;
  const bool print_trait_raw_data; /**< If true, will print raw trait data for each gen and replicate */
  /** First element is index of trait in \p trait_freq; second element is number of traits to track */
  const std::vector<int> trait_info;
};
/**
   @brief Struct containing constant parameter values (fixed for all models)
*/
struct Fixed_Parameters {
  inline static const double tolerance = fixed_parameters::tolerance; /**< Tolerance for comparing doubles */
   /** Number of independent stochastic replicates that are run */
  inline static const int number_replicates = fixed_parameters::number_replicates;
  /** Max number of generations to run a replicate (prevents a stable polymorphism causing an infinite loop) */
  inline static const int max_generations_per_sim = fixed_parameters::max_generations_per_sim;
  /** Memory allocation for Data::_trait_freq_by_gen */
  inline static const int reserve_memory_trait_freq = fixed_parameters::reserve_memory_trait_freq;
};
/**
   @brief Struct for parameters of the Haploid Single Environment model
*/
struct HSE_Model_Parameters {
  Shared_Parameters shared;
  struct HSE_Specific_Parameters {
    const double selection_coefficient; /**< Selection coefficient of allele A */
  } model;
  Fixed_Parameters fixed;
};

/**
   @brief Struct for parameters of the Diploid Single Environment model
*/
struct DSE_Model_Parameters {
  Shared_Parameters shared;
  struct DSE_Specific_Parameters {
    const double selection_coefficient_homozygote; /**< Selection coefficient of the AA genotype */
    const double selection_coefficient_heterozygote; /**< Selection coefficient of the Aa genotype */
  } model;
  Fixed_Parameters fixed;
};
/**
   @brief Struct for parameters of the Haploid Two Environments model
*/
struct HTE_Model_Parameters {
  Shared_Parameters shared;
  struct HTE_Specific_Parameters {
    const double selection_coefficient_A_env_1; /**< Selection coefficient of the A allele in environment 1 */
    const double selection_coefficient_A_env_2; /**< Selection coefficient of the A allele in environment 2 */
    const double selection_coefficient_a_env_1; /**< Selection coefficient of the a allele in environment 1 */
    const double selection_coefficient_a_env_2; /**< Selection coefficient of the a allele in environment 2 */
    const int gen_env_1; /**< Number of generations spent in environment 1 (the initial environment) */
  } model;
  Fixed_Parameters fixed;
};

/**
   @brief Struct for parameters of the Haploid Two Effects One Environment model
*/
struct HTEOE_Model_Parameters {
  Shared_Parameters shared;
  struct HTEOE_Specific_Parameters {
    const double selection_coefficient_A1; /**< Selection coefficient of the A allele's 1st effect */
    const double selection_coefficient_A2; /**< Selection coefficient of the A allele's 2nd effect */
    const double selection_coefficient_a1; /**< Selection coefficient of the a allele's 1st effect */
    const double selection_coefficient_a2; /**< Selection coefficient of the a allele's 2nd effect */
  } model;
  Fixed_Parameters fixed;
};

#endif

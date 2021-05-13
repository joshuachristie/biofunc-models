/**
   @file HTE.h
   @brief Contains functions to run the Haploid Two Environments model
*/

#ifndef HTE_H
#define HTE_H

#include <vector>
#include <random>
#include "Parameters.h"

/**
   @brief Namespace for Haploid Two Environments
   @details This model explores biological function in the context of a trait that experiences two different
   environments and has a separate fitness in each environment.
   The set-up is a Wright-Fisher haploid model: one locus with two alleles; one environment with two types.
   At any given time, there is one environment. The environment can change during the simulation however.
   Let's say that the resident population has trait (allele) a which is invaded by our allele of interest, A.
   Each allele (A/a) has two components to its fitness: fitness in env 1 and env 2.
   Trait/allele a is the resident trait and in environment 1 its fitness is 
   1 + parameters::HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_a_env_1 while in environment 2 
   its fitness is 1 + parameters::HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_a_env_2.
   Trait/allele A is our trait of interest and in environment 1 its fitness is 
   1 + parameters::HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_A_env_1 while in environment 2
   its fitness is 1 + parameters::HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_A_env_2.

   It requires specification of a particular selection regime (i.e. environment at each generation).
   For simplicity, I only consider a selection regime in which environment 1 exists for 
   parameters::HTE_Model_Parameters::HTE_Specific_Parameters::gen_env_1 generations, which is followed by environment 2 
   existing for the remaining generations (until fixation or loss). One interpretation of this model is the 
   evolved for/maintained by distinction (where evolved for is env 1 and maintained by is env 2).
*/
namespace HTE {

  const parameters::HTE_Model_Parameters parse_parameter_values(int argc, char* argv[]);
  
  const std::vector<double> get_fitness_function(const parameters::HTE_Model_Parameters &parameters);

  void calculate_allele_freqs(std::vector<double> &trait_freq, const std::vector<double> &fitnesses,
			      const parameters::HTE_Model_Parameters &parameters, std::mt19937 &rng, int &gen);
  
  void run_model(int argc, char* argv[]);

}

#endif

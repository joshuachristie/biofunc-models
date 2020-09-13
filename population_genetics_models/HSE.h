/**
   @file HSE.h
   @brief Contains functions to run the Haploid Single Environment model
*/

#ifndef HSE_H
#define HSE_H

#include <vector>
#include <random>
#include "Parameters.h"

/**
   @brief Namespace for Haploid Single Environment
   @details This model explores biological function in the simplest case: a single trait with a single effect.
   The set-up is a Wright-Fisher haploid model: one locus, two alleles, (and implicitly a single environment).
   For this model, there is no need to apportion the function metric between components using partial
   information decomposition because there is a simple mapping between trait -> selection_coefficient ->
   persistence_probability (i.e. function is wholly attributable to the trait's selection coefficient).
*/
namespace HSE {
  /**
     @brief Reads in parameter values from command line into an HSE_Model_Parameters object
     @param[in] argc Number of commandline arguments
     @param[in] argv Command line arguments
     @return params HSE_Model_Parameters struct
  */
  const HSE_Model_Parameters parse_parameter_values(int argc, char* argv[]);
  /**
     @brief Calculates fitness function
     @param[in] parameters HSE_Model_Parameters::HSE_Specific_Parameters::selection_coefficient
     @return fitnesses A vector of length 2 containing the fitnesses of the A and a alleles [wA, wa]
  */
  const std::vector<double> get_fitness_function(const HSE_Model_Parameters &parameters);
  /**
     @brief Calculates frequency of the A allele after selection
     @param[in, out] allele_A_freq The frequency of the A allele
     @param[in] fitnesses A vector containing the fitnesses of the A and a alleles [wA, wa]
     @param[in] HSE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @param[in, out] gen Current generation
     @return Nothing (but modifies \p allele_A_freq and \p gen)
  */
  void calculate_allele_freqs(double &allele_A_freq, const std::vector<double> &fitnesses,
			      const HSE_Model_Parameters &parameters, std::mt19937 &rng, int &gen);
  /**
     @brief Runs Haploid Single Environment model
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @return Nothing (but prints results)
  */
  void run_model(int argc, char* argv[]);

}

#endif

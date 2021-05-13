/**
   @file DSE.h
   @brief Contains functions to run the Diploid Single Environment model
**/
#ifndef DSE_H
#define DSE_H

#include <vector>
#include <random>
#include "Parameters.h"
/**
   @brief Namespace for Diploid Single Environment
   @details This model is a Wright-Fisher diploid model (one locus, two alleles, single environment).
   The resident population comprises individuals with trait (genotype) aa, which is "invaded" by our
   allele of interest, A (i.e. at initialisation, a single aa individual mutates to an Aa individual).
   This set-up allows us explores biological function in the context of a trait whose spread is influenced, not
   only by its own effects, but also by the effects of another trait. Let's say that our trait of interest is
   the heterozygote Aa (similar logic applies if we were interested in the homozygote AA). Whether Aa persists
   in the population depends not only on the fitness of Aa but also on the fitness of AA (since they are linked
   through allele A).
*/

namespace DSE {
  
  const parameters::DSE_Model_Parameters parse_parameter_values(int argc, char* argv[]);
  
  const std::vector<double> get_fitness_function(const parameters::DSE_Model_Parameters &parameters);

  void calculate_allele_freqs(std::vector<double> &trait_freq, const std::vector<double> &fitnesses,
			      const parameters::DSE_Model_Parameters &parameters, std::mt19937 &rng, int &gen);

  void run_model(int argc, char* argv[]);

}

#endif

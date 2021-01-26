/**
   @file DSE.h
   @brief Contains functions to run the Diploid Single Environment model
   @details This model explores biological function in the context of a trait whose spread is influenced, not
   only by its own effects, but also by the effects of another trait.
   The set-up is a Wright-Fisher diploid model: one locus, two alleles, (and implictly a single environment).
   The resident population has trait (genotype) aa, which is "invaded" by our allele of interest, A (i.e. an aa
   individual mutates to Aa).
   Let's say that our trait of interest is the heterozygote Aa.
   Whether A allele persists in the population depends not only on the fitness of Aa but also on the fitness of AA
   Since I use allele frequencies---not genotype frequencies---to calculate the function metric, the function of
   allele A needs to be apportioned over the traits (genotypes) Aa and AA (as well as their
   redundant/synergistic effects).
   For this I use PID, where the two "sources" are the fitnesses of Aa and AA.
   By combining the function metric and the PID, we can calculate the function of traits Aa and AA.

**/
#ifndef DSE_H
#define DSE_H

#include <vector>
#include <random>
#include "Parameters.h"
/**
   @brief Namespace for Diploid Single Environment
*/

namespace DSE {
  
  const DSE_Model_Parameters parse_parameter_values(int argc, char* argv[]);
  
  const std::vector<double> get_fitness_function(const DSE_Model_Parameters &parameters);

  void calculate_allele_freqs(double &allele_A_freq, const std::vector<double> &fitnesses,
			      const DSE_Model_Parameters &parameters, std::mt19937 &rng, int &gen);

  void run_model(int argc, char* argv[]);

}

#endif

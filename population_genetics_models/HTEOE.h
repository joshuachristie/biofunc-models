/**
   @file HTEOE.h
   @brief Contains functions to run the Haploid Two Effects One Environment model
*/

#ifndef HTEOE_H
#define HTEOE_H

#include <vector>
#include <random>
#include "Parameters.h"

/**
   @brief Namespace for Haploid Two Effects One Environment
   @details This model explores biological function in the context of a trait that causes two phenotypic
   (fitness) effects. The step-up is a Wright-Fisher haploid model: two loci, each with two alleles (and
   implicitly a single environment). (But as it is a haploid model, and there is no recombination between the
   loci, the two loci are linked, and there are only two possible configurations for a given simulation.)
   The resident population has trait (allele) a which is invaded by our allele of interest, A.
   Trait A has two effects, given by selection_coefficient_A1 and selection_coefficient_A2.
   The resident trait, a, also has two effects, selection_coefficient_a1 and selection_coefficient_a2.
   The fitness of the trait of interest is  wA = 1.0 + selection_coefficient_A1 + selection_coefficient_A2.
   The fitness of the resident trait is  wa = 1.0 + selection_coefficient_a1 + selection_coefficient_a2.
*/
namespace HTEOE {
  
  const HTEOE_Model_Parameters parse_parameter_values(int argc, char* argv[]);
  
  const std::vector<double> get_fitness_function(const HTEOE_Model_Parameters &parameters);

  void calculate_allele_freqs(double &allele_A_freq, const std::vector<double> &fitnesses,
			      const HTEOE_Model_Parameters &parameters, std::mt19937 &rng, int &gen);

  void run_model(int argc, char* argv[]);

}

#endif

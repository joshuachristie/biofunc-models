/**
   @file HTE.cpp
   @brief Contains functions to run the Haploid Two Environments model
   @detail This model explores biological function in the context of a trait that experiences two different
   environments and has a separate fitness in each environment.
   The set-up is a Wright-Fisher haploid model: one locus with two alleles; one environment with two types
   At any given time, there is one environment. The environment can change during the simulation however.
   Let's say that the resident population has trait (allele) a which is invaded by our allele of interest, A
   Each allele (A/a) has two components to its fitness: fitness in env 1 and env 2.
   Trait/allele a is the resident trait and in environment 1 its fitness is 1 + selection_coefficient_a_env_1
   while in environment 2 its fitness is 1 + selection_coefficient_a_env_2.
   Trait/allele A is our trait of interest and in environment 1 its fitness is 1 + selection_coefficient_A_env_1
   while in environment 2 its fitness is 1 + selection_coefficient_A_env_2.
   For PID purposes, I will fix the fitness of trait a (and thus also fix selection_coefficient_a_env_1 and
   selection_coefficient_a_env_2).
   The two "sources" are thus selection_coefficient_A_env_1 and selection_coefficient_A_env_2.
   PID is used to apportion trait A's function between selection_coefficient_A_env_1 and
   selection_coefficient_A_env_2 (and to their redundant/synergistic effects).

   It requires specification of a particular selection regime (i.e. environment at each generation).
   For simplicity, I only consider a selection regime in which environment 1 exists for gen_env_1 generations,
   which is followed by environment 2 existing for the remaining generations (until fixation or loss).
   One interpretation of this model is the evolved for/maintained by distinction (where evolved for is env 1
   and maintained by is env 2)
   
   Note that here PID is not strictly necessary in order to calculate the overall function of trait/allele A
   (since both sources are ultimately derived from trait/allele A), but it is necessary in order to apportion
   the function of trait/allele A between its effects.
   If a trait evolves in environment 1 but is maintained due to its effects in environment 2, we might want to
   quantify the function of the trait's effect in environment 1 separately from its effect in environment 2
*/

#include <cassert>
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include "Parameters.h"
#include "HTE.h"
#include "rng.h"
#include "persistence_probability.h"
#include "print_results.h"
#include "allele_invasion.h"
#include "DataContainer.h"

/**
   @brief Namespace for Haploid Two Environments
*/
namespace HTE {
  /**
     @brief Reads in parameter values from command line into a struct
     @param[in] argc Number of commandline arguments
     @param[in] argv Command line arguments
     @return params HTE_Model_Parameters struct
  */
  const HTE_Model_Parameters parse_parameter_values(int argc, char* argv[]){
    assert(std::string(argv[1]) == "HTE");
    assert(argc == 11 && "The HTE model must have 10 command line arguments (the first must be 'HTE')");
    const int population_size = atoi(argv[2]);
    const double selection_coefficient_A_env_1 = atof(argv[3]);
    const double selection_coefficient_A_env_2 = atof(argv[4]);
    const double selection_coefficient_a_env_1 = atof(argv[5]);
    const double selection_coefficient_a_env_2 = atof(argv[6]);
    const int gen_env_1 = atoi(argv[7]);
    const double initial_A_freq = 1.0 / static_cast<double>(population_size); // 1/N
    const int number_reinvasions = atoi(argv[8]);
    const int number_gens_to_output_pp = atoi(argv[9]);
    const bool print_allele_A_raw_data = static_cast<bool>(atoi(argv[10]));
    const HTE_Model_Parameters params {{population_size, initial_A_freq, number_reinvasions,
	number_gens_to_output_pp, print_allele_A_raw_data}, {selection_coefficient_A_env_1,
					 selection_coefficient_A_env_2, selection_coefficient_a_env_1,
					 selection_coefficient_a_env_2, gen_env_1}};
    return params;
  }
  /**
     @brief Calculates fitness function
     @param[in] parameters HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_A_env_1 - selection coefficient of the A allele in environment 1
     @param[in] parameters HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_A_env_2 - selection coefficient of the A allele in environment 2
     @param[in] parameters HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_a_env_1 - selection coefficient of the a allele in environment 1
     @param[in] parameters HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_a_env_2 - selection coefficient of the a allele in environment 2
     @return fitnesses A vector of length 4 containing the fitnesses of the A and a alleles in environments 1 and 2 [wA_1, wA_2, wa_1, wa_2]
  */
  const std::vector<double> get_fitness_function(const HTE_Model_Parameters &parameters){
    // haploid_fitnesses[w_A_1, w_A_2, w_a_1, w_a_2] gives relative fitness of A/a alleles in envs 1/2
    const std::vector<double> fitnesses {1.0 + parameters.model.selection_coefficient_A_env_1,
      1.0 + parameters.model.selection_coefficient_A_env_2,
      1.0 + parameters.model.selection_coefficient_a_env_1,
      1.0 + parameters.model.selection_coefficient_a_env_2};
    return fitnesses;
  }
  /**
     @brief Calculates frequency of the A allele after selection
     @param[in, out] allele_A_freq The frequency of the A allele
     @param[in] fitnesses A vector containing the fitnesses of the A and a alleles [wA, wa]
     @param[in] HTE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @param[in, out] gen Current generation of the simulation
     @return Nothing (but modifies \p allele_A_freq)
  */
  void calculate_allele_freqs(double &allele_A_freq, const std::vector<double> &fitnesses,
			      const HTE_Model_Parameters &parameters, std::mt19937 &rng, int &gen){
    std::vector<double> expected_allele_freq_raw(2);
    expected_allele_freq_raw[0] = allele_A_freq * fitnesses[gen >= parameters.model.gen_env_1];
    expected_allele_freq_raw[1] = (1.0 - allele_A_freq) * fitnesses[2 + (gen >= parameters.model.gen_env_1)];
    // calculate normalised expectation of allele_A_freq
    allele_A_freq = expected_allele_freq_raw[0] / (std::accumulate(expected_allele_freq_raw.begin(),
								   expected_allele_freq_raw.end(), 0.0));
    // sample to get realised allele_A_freq
    std::binomial_distribution<int> surviving_As(parameters.shared.population_size, allele_A_freq);
    allele_A_freq =
      static_cast<double>(surviving_As(rng)) / static_cast<double>(parameters.shared.population_size);
    ++gen;
  }
  /**
     @brief Runs Haploid Two Effects model
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @return Nothing (but prints results)
  */
  void run_model(int argc, char* argv[]){
    std::mt19937 rng = initialise_rng();
    const HTE_Model_Parameters params = parse_parameter_values(argc, argv);
    const std::vector<double> fitnesses = get_fitness_function(params);
    DataContainer data(params.fixed.number_replicates, params.shared.number_gens_to_output_pp,
		       params.fixed.reserve_memory_allele_freq);
    calculate_persistence_probability(params, rng, fitnesses, calculate_allele_freqs, data);
    print::print_results(argc, argv, data, params);
  }
  
}

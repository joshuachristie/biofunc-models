/**
   @file HTE.cpp
   @brief Contains functions to run the Haploid Two Environments model
*/

#include <cassert>
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include "Parameters.h"
#include "helper_functions.h"
#include "HTE.h"
#include "rng.h"
#include "persistence_probability.h"
#include "print_results.h"
#include "allele_invasion.h"

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
  HTE_Model_Parameters parse_parameter_values(int argc, char* argv[]){
    assert(std::string(argv[1]) == "HTE");
    assert(argc == 9 && "The HTE model must have 8 command line arguments (the first must be 'HTE')");
    const int population_size = atoi(argv[2]);
    const double selection_coefficient_A_env_1 = atof(argv[3]);
    const double selection_coefficient_A_env_2 = atof(argv[4]);
    const double selection_coefficient_a_env_1 = atof(argv[5]);
    const double selection_coefficient_a_env_2 = atof(argv[6]);
    const int gen_env_1 = atoi(argv[7]);
    const double initial_A_freq = 1.0 / static_cast<double>(population_size); // 1/N
    const int number_reinvasions = atoi(argv[8]);
    HTE_Model_Parameters params {{population_size}, {selection_coefficient_A_env_1,
	selection_coefficient_A_env_2, selection_coefficient_a_env_1,
	selection_coefficient_a_env_2, gen_env_1, initial_A_freq, number_reinvasions}};
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
  std::vector<double> get_fitness_function(const HTE_Model_Parameters &parameters){
    // haploid_fitnesses[w_A_1, w_A_2, w_a_1, w_a_2] gives relative fitness of A/a alleles in envs 1/2
    std::vector<double> fitnesses {1.0 + parameters.model.selection_coefficient_A_env_1,
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
    HTE_Model_Parameters params = parse_parameter_values(argc, argv);
    std::vector<double> fitnesses = get_fitness_function(params);
    std::vector<bool> final_A_freqs;
    final_A_freqs.reserve(params.fixed.number_replicates);
    double persistence_probability = calculate_persistence_probability(params, rng, fitnesses,
								       final_A_freqs, calculate_allele_freqs);
    print::print_persistence_probability(argc, argv, persistence_probability);
  }
  
}

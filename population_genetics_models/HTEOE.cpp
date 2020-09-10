/**
   @file HTEOE.cpp
   @brief Contains functions to run the Haploid Two Effects One Environment model
*/

#include <cassert>
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include "Parameters.h"
#include "helper_functions.h"
#include "HTEOE.h"
#include "rng.h"
#include "persistence_probability.h"
#include "print_results.h"
#include "allele_invasion.h"
#include "DataContainer.h"

/**
   @brief Namespace for Haploid Two Effects One Environment
*/
namespace HTEOE {
  /**
     @brief Reads in parameter values from command line into a struct
     @param[in] argc Number of commandline arguments
     @param[in] argv Command line arguments
     @return params HTEOE_Model_Parameters struct
  */
  const HTEOE_Model_Parameters parse_parameter_values(int argc, char* argv[]){
    assert(std::string(argv[1]) == "HTEOE");
    assert(argc == 10 && "The HTEOE model must have 9 command line arguments (the first must be 'HTEOE')");
    const int population_size = atoi(argv[2]);
    const double selection_coefficient_A1 = atof(argv[3]);
    const double selection_coefficient_A2 = atof(argv[4]);
    const double selection_coefficient_a1 = atof(argv[5]);
    const double selection_coefficient_a2 = atof(argv[6]);
    const double initial_A_freq = 1.0 / static_cast<double>(population_size); // 1/N
    const int number_reinvasions = atoi(argv[7]);
    const int number_gens_to_output_pp = atoi(argv[8]);
    const bool print_allele_A_raw_data = static_cast<bool>(atoi(argv[9]));
    const HTEOE_Model_Parameters params {{population_size, initial_A_freq, number_reinvasions,
	number_gens_to_output_pp, print_allele_A_raw_data}, {selection_coefficient_A1, selection_coefficient_A2,
					   selection_coefficient_a1, selection_coefficient_a2}};
    return params;
  }
  /**
     @brief Calculates fitness function
     @param[in] parameters HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_A1 - selection coefficient of the A allele's first effect
     @param[in] parameters HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_A2 - selection coefficient of the A allele's second effect
     @param[in] parameters HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_a1 - selection coefficient of the a allele's first effect
     @param[in] parameters HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_a2 - selection coefficient of the a allele's second effect
     @return fitnesses A vector of length 2 containing the fitnesses for the A and a alleles [wA, wa]
  */
  const std::vector<double> get_fitness_function(const HTEOE_Model_Parameters &parameters){
    const std::vector<double> fitnesses
      {1.0 + parameters.model.selection_coefficient_A1 + parameters.model.selection_coefficient_A2,
       1.0 + parameters.model.selection_coefficient_a1 + parameters.model.selection_coefficient_a2}; 
    return fitnesses;
  }
  /**
     @brief Calculates frequency of the A allele after selection
     @param[in, out] allele_A_freq The frequency of the A allele
     @param[in] fitnesses A vector containing the fitnesses of the A and a alleles [wA, wa]
     @param[in] HTEOE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @return Nothing (but modifies \p allele_A_freq)
  */
  void calculate_allele_freqs(double &allele_A_freq, const std::vector<double> &fitnesses,
			      const HTEOE_Model_Parameters &parameters, std::mt19937 &rng, int &gen){
    std::vector<double> expected_allele_freq_raw(2);
    expected_allele_freq_raw[0] = allele_A_freq * fitnesses[0];
    expected_allele_freq_raw[1] = (1.0 - allele_A_freq) * fitnesses[1];
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
     @brief Runs Haploid Two Effects One Environment model
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @return Nothing (but prints results)
  */
  void run_model(int argc, char* argv[]){
    std::mt19937 rng = initialise_rng();
    const HTEOE_Model_Parameters params = parse_parameter_values(argc, argv);
    const std::vector<double> fitnesses = get_fitness_function(params);
    DataContainer data(params.fixed.number_replicates, params.shared.number_gens_to_output_pp,
		       params.fixed.reserve_memory_allele_freq);
    calculate_persistence_probability(params, rng, fitnesses, calculate_allele_freqs, data);
    print::print_results(argc, argv, data, params);
  }

}

#include <cassert>
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include "Parameters.h"
#include "HTEOE.h"
#include "rng.h"
#include "conditional_existence_probability.h"
#include "trait_invasion.h"
#include "run_scenario.h"

namespace HTEOE {

  /**
     @brief Reads in parameter values from command line into a struct
     @param[in] argc Number of commandline arguments
     @param[in] argv Command line arguments
     @return params parameters::HTEOE_Model_Parameters struct
  */
  const parameters::HTEOE_Model_Parameters parse_parameter_values(int argc, char* argv[]){
    assert(std::string(argv[1]).compare("HTEOE") == 0);
    assert(argc == 9 && "The HTEOE model must have 8 command line arguments (the first must be 'HTEOE')");
    assert(std::string(argv[2]).compare("QEF") == 0 &&
	   "Incorrect model specification: for HTEOE, the second arg must be QEF");
    const int population_size = atoi(argv[3]);
    const double selection_coefficient_A1 = atof(argv[4]);
    const double selection_coefficient_A2 = atof(argv[5]);
    const double selection_coefficient_a1 = atof(argv[6]);
    const double selection_coefficient_a2 = atof(argv[7]);
    const double initial_trait_freq = 1.0 / static_cast<double>(population_size);
    const int number_reinvasions = atoi(argv[8]);
    const std::vector<int> trait_info {0, 1};
    const parameters::HTEOE_Model_Parameters params {{population_size, initial_trait_freq, number_reinvasions, trait_info},
					 {selection_coefficient_A1, selection_coefficient_A2,
					  selection_coefficient_a1, selection_coefficient_a2}};
    return params;
  }
  /**
     @brief Calculates fitness function
     @param[in] parameters parameters::HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_A1 - selection coefficient of the A allele's first effect
     @param[in] parameters parameters::HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_A2 - selection coefficient of the A allele's second effect
     @param[in] parameters parameters::HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_a1 - selection coefficient of the a allele's first effect
     @param[in] parameters parameters::HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_a2 - selection coefficient of the a allele's second effect
     @return fitnesses A vector of length 2 containing the fitnesses for the A and a alleles [wA, wa]
  */
  const std::vector<double> get_fitness_function(const parameters::HTEOE_Model_Parameters &parameters){
    const std::vector<double> fitnesses
      {1.0 + parameters.model.selection_coefficient_A1 + parameters.model.selection_coefficient_A2,
       1.0 + parameters.model.selection_coefficient_a1 + parameters.model.selection_coefficient_a2}; 
    return fitnesses;
  }
  /**
     @brief Calculates frequency of the trait after selection
     @param[in, out] trait_freq The frequency of the trait (A allele)
     @param[in] fitnesses A vector containing the fitnesses of the A and a alleles [wA, wa]
     @param[in] parameters::HTEOE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @param[in, out] gen The current generation
     @return Nothing (but modifies \p trait_freq and increments \p gen)
  */
  void calculate_trait_freqs(std::vector<double> &trait_freq, const std::vector<double> &fitnesses,
			     const parameters::HTEOE_Model_Parameters &parameters, std::mt19937 &rng, int &gen){
    std::vector<double> expected_allele_freq_raw(2);
    expected_allele_freq_raw[0] = trait_freq[0] * fitnesses[0];
    expected_allele_freq_raw[1] = (1.0 - trait_freq[0]) * fitnesses[1];
    // calculate normalised expectation of trait_freq
    double expectation = expected_allele_freq_raw[0] / (std::accumulate(expected_allele_freq_raw.begin(),
									expected_allele_freq_raw.end(), 0.0));
    // sample to get realised trait_freq
    std::binomial_distribution<int> surviving_As(parameters.shared.population_size, expectation);
    trait_freq[0] =
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
    std::mt19937 rng = rng::initialise_rng();
    const parameters::HTEOE_Model_Parameters params = parse_parameter_values(argc, argv);
    const std::vector<double> fitnesses = get_fitness_function(params);
    run_scenario::QEF(params, rng, fitnesses, calculate_trait_freqs, argv, argc);
  }

}

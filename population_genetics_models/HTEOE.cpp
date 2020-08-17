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
#include <iostream> // will probably need to be deleted once after refactoring and addition of print methods

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
  HTEOE_Model_Parameters parse_parameter_values(int argc, char* argv[]){
    assert(std::string(argv[1]) == "HTEOE");
    assert(argc == 8 && "The HTEOE model must have 7 command line arguments (the first must be 'HTEOE')");
    const int population_size = atoi(argv[2]);
    const double selection_coefficient_A1 = atof(argv[3]);
    const double selection_coefficient_A2 = atof(argv[4]);
    const double selection_coefficient_a1 = atof(argv[5]);
    const double selection_coefficient_a2 = atof(argv[6]);
    const int number_generations = atoi(argv[7]);
    HTEOE_Model_Parameters params {{population_size}, {selection_coefficient_A1, selection_coefficient_A2,
	selection_coefficient_a1, selection_coefficient_a2, number_generations}};
    return params;
  }
  /**
     @brief Calculates fitness function
     @param[in] parameters HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_A1 - selection coefficient of the A allele's first effect
     @param[in] parameters HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_A2 - selection coefficient of the A allele's second effect
     @param[in] parameters HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_a1 - selection coefficient of the a allele's first effect
     @param[in] parameters HTEOE_Model_Parameters::HTEOE_Specific_Parameters::selection_coefficient_a2 - selection coefficient of the a allele's second effect
     @return haploid_fitnesses A vector of length 2 containing the fitnesses for the A and a alleles [wA, wa]
  */
  std::vector<double> get_fitness_function(const HTEOE_Model_Parameters &parameters){
    std::vector<double> haploid_fitnesses
      {1.0 + parameters.model.selection_coefficient_A1 + parameters.model.selection_coefficient_A2,
      1.0 + parameters.model.selection_coefficient_a1 + parameters.model.selection_coefficient_a2}; 
    return haploid_fitnesses;
  }
  /**
     @brief Calculates expected frequency of the A allele after selection
     @param[in, out] allele_A_freq The frequency of the A allele
     @param[in] haploid_fitnesses A vector containing the fitnesses of the A and a alleles [wA, wa]
     @return Nothing (but modifies \p allele_A_freq)
*/
  void expected_allele_freqs(double &allele_A_freq, const std::vector<double> &haploid_fitnesses){
      std::vector<double> expected_allele_freq_raw(2);
  expected_allele_freq_raw[0] = allele_A_freq * haploid_fitnesses[0];
  expected_allele_freq_raw[1] = (1.0 - allele_A_freq) * haploid_fitnesses[1];
  // replace allele_A_freq with normalised expectation
  allele_A_freq = expected_allele_freq_raw[0] / (std::accumulate(expected_allele_freq_raw.begin(),
								 expected_allele_freq_raw.end(), 0.0));
  }
  /**
     @brief Calculates realised frequency of the A allele (by binomial sampling of expectated number of As)
     @param[in, out] allele_A_freq Frequency of allele A
     @param[in] HTEOE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @return Nothing (but modifies \p allele_A_freq)
*/
  void realised_allele_freqs(double &allele_A_freq, const HTEOE_Model_Parameters &parameters, std::mt19937 &rng){
    std::binomial_distribution<int> surviving_As(parameters.shared.population_size, allele_A_freq);
    allele_A_freq =
      static_cast<double>(surviving_As(rng)) / static_cast<double>(parameters.shared.population_size);
  }
  /**
     @brief Runs one replicate of the simulation
     @param[in] haploid_fitnesses Vector containing the fitnesses of the A and a alleles in the two environments [wA_1, wA_2, wa_1, wa_2]
     @param[in] HTEOE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in] HTEOE_Model_Parameters::Shared_Parameters::number_generations Maximum number of generations for which the simulation will run (for the case in which a trait invades the ancestral population and is not challenged afterwards; for the case in which the fixed trait must withstand invaders, each step runs until an absorbing state is reached).
     @param[in, out] rng Random number generator
     @param[in, out] final_A_freqs Vector of bools storing whether allele A persists (true/false) at census
     @param[in] HTEOE_Model_Parameters::Fixed_Parameters::tolerance Tolerance for comparing equality of doubles
     @return Nothing (but alters \p final_A_freqs)
*/
  void run_simulation(const std::vector<double> &haploid_fitnesses, const HTEOE_Model_Parameters &parameters,
		      std::mt19937 &rng, std::vector<bool> &final_A_freqs){
    double allele_A_freq = 1.0 / static_cast<double>(parameters.shared.population_size); // initial freq is 1/N
    int gen = 0;
    while (gen < parameters.model.number_generations &&
	   !(close_to_value(allele_A_freq, 0.0, parameters.fixed.tolerance) ||
					 close_to_value(allele_A_freq, 1.0, parameters.fixed.tolerance))){
      expected_allele_freqs(allele_A_freq, haploid_fitnesses);
      realised_allele_freqs(allele_A_freq, parameters, rng);
      ++gen;
    }
    close_to_value(allele_A_freq, 0.0, parameters.fixed.tolerance) ? final_A_freqs.push_back(0) :
      final_A_freqs.push_back(1);
  }
    /**
     @brief Runs Haploid Two Effects One Environment model
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @return Nothing (but prints results)
*/
  void run_model(int argc, char* argv[]){
    std::mt19937 rng = initialise_rng();
    HTEOE_Model_Parameters params = parse_parameter_values(argc, argv);
    std::vector<double> haploid_fitnesses = get_fitness_function(params);
    std::vector<bool> final_A_freqs;
    final_A_freqs.reserve(params.fixed.number_replicates);
    // will change how I implement running of multiple replicates, but leaving it for post-refactor extension
    for (int rep = 0; rep < params.fixed.number_replicates; rep++){
      run_simulation(haploid_fitnesses, params, rng, final_A_freqs);
    }
    // will alter output in a later extension (will print to file directly from c++ rather than via the python run script)
    std::cout << std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) / static_cast<double>(params.fixed.number_replicates) << std::endl;

  }

}

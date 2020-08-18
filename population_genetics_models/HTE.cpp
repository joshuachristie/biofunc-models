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
#include <iostream> // will probably need to be deleted once after refactoring and addition of print methods

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
    const int gen_env_2 = atoi(argv[8]);
    HTE_Model_Parameters params {{population_size}, {selection_coefficient_A_env_1,
	selection_coefficient_A_env_2, selection_coefficient_a_env_1,
	selection_coefficient_a_env_2, gen_env_1, gen_env_2}};
    return params;
  }
  /**
     @brief Calculates fitness function
     @param[in] parameters HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_A_env_1 - selection coefficient of the A allele in environment 1
     @param[in] parameters HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_A_env_2 - selection coefficient of the A allele in environment 2
     @param[in] parameters HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_a_env_1 - selection coefficient of the a allele in environment 1
     @param[in] parameters HTE_Model_Parameters::HTE_Specific_Parameters::selection_coefficient_a_env_2 - selection coefficient of the a allele in environment 2
     @return haploid_fitnesses A vector of length 4 containing the fitnesses of the A and a alleles in environments 1 and 2 [wA_1, wA_2, wa_1, wa_2]
  */
  std::vector<double> get_fitness_function(const HTE_Model_Parameters &parameters){
    // haploid_fitnesses[w_A_1, w_A_2, w_a_1, w_a_2] gives relative fitness of A/a alleles in envs 1/2
    std::vector<double> haploid_fitnesses {1.0 + parameters.model.selection_coefficient_A_env_1,
      1.0 + parameters.model.selection_coefficient_A_env_2,
      1.0 + parameters.model.selection_coefficient_a_env_1,
      1.0 + parameters.model.selection_coefficient_a_env_2};
    return haploid_fitnesses;
  }
  /**
     @brief Calculates expected frequency of the A allele after selection
     @param[in, out] allele_A_freq The frequency of the A allele
     @param[in] haploid_fitnesses A vector containing the fitnesses of the A and a alleles [wA, wa]
     @return Nothing (but modifies \p allele_A_freq)
*/
  void expected_allele_freqs(double &allele_A_freq, const std::vector<double> &haploid_fitnesses,
			     int env_state){
  std::vector<double> expected_allele_freq_raw(2);
  expected_allele_freq_raw[0] = allele_A_freq * haploid_fitnesses[env_state];
  expected_allele_freq_raw[1] = (1.0 - allele_A_freq) * haploid_fitnesses[2 + env_state];
  // replace allele_A_freq with normalised expectation
  allele_A_freq = expected_allele_freq_raw[0] / (std::accumulate(expected_allele_freq_raw.begin(),
								 expected_allele_freq_raw.end(), 0.0));
}
  /**
     @brief Calculates realised frequency of the A allele (by binomial sampling of expectated number of As)
     @param[in, out] allele_A_freq Frequency of allele A
     @param[in] HTE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @return Nothing (but modifies \p allele_A_freq)
*/
  void realised_allele_freqs(double &allele_A_freq, const HTE_Model_Parameters &parameters, std::mt19937 &rng){
    std::binomial_distribution<int> surviving_As(parameters.shared.population_size, allele_A_freq);
    allele_A_freq =
      static_cast<double>(surviving_As(rng)) / static_cast<double>(parameters.shared.population_size);
  }
  /**
     @brief Runs one replicate of the simulation
     @param[in] haploid_fitnesses Vector containing the fitnesses of the A and a alleles in the two environments [wA_1, wA_2, wa_1, wa_2]
     @param[in] HTE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in] HTE_Model_Parameters::HTE_Specific_Parameters::gen_env_1 Number of generations spent in environment 1
     @param[in] HTE_Model_Parameters::HTE_Specific_Parameters::gen_env_2 Number of generations spent in environment 2
     @param[in] HTE_Model_Parameters::Fixed_Parameters::tolerance Tolerance for comparing equality of doubles
     @param[in, out] rng Random number generator
     @param[in, out] final_A_freqs Vector of bools storing whether allele A persists (true/false) at census
     @return Nothing (but alters \p final_A_freqs)
*/
  void run_simulation(const std::vector<double> &haploid_fitnesses, const HTE_Model_Parameters &parameters,
		      std::mt19937 &rng, std::vector<bool> &final_A_freqs){
    double allele_A_freq = 1.0 / static_cast<double>(parameters.shared.population_size); // initial freq is 1/N
    int gen = 0;
    int env_state = 0; // env 1 = 0; env 2 = 1
    while (gen < (parameters.model.gen_env_1 + parameters.model.gen_env_2) &&
	   !(close_to_value(allele_A_freq, 0.0, parameters.fixed.tolerance) ||
					 close_to_value(allele_A_freq, 1.0, parameters.fixed.tolerance))){
      expected_allele_freqs(allele_A_freq, haploid_fitnesses, env_state);
      realised_allele_freqs(allele_A_freq, parameters, rng);
      ++gen;
      if (gen == parameters.model.gen_env_1 - 1) {env_state = 1;}
    }
    close_to_value(allele_A_freq, 0.0, parameters.fixed.tolerance) ? final_A_freqs.push_back(0) :
      final_A_freqs.push_back(1);
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
    std::vector<double> haploid_fitnesses = get_fitness_function(params);
    std::vector<bool> final_A_freqs;
    final_A_freqs.reserve(params.fixed.number_replicates);
    calculate_persistence_probability(params, run_simulation, rng, haploid_fitnesses, final_A_freqs);
    // will alter output in a later extension (will print to file directly from c++ rather than via the python run script)
    std::cout << std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) / static_cast<double>(params.fixed.number_replicates) << std::endl;

  }
  
}

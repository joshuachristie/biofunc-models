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
#include "include/helper_functions.h"
#include "HTE.h"
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
    assert(argc == 9);
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
     @param[in] allele_A_freq Frequency of allele A
     @param[in] haploid_fitnesses Vector containing the fitnesses of the A and a alleles in the two environments [wA_1, wA_2, wa_1, wa_2]
     @param[in] HTE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in] HTE_Model_Parameters::HTE_Specific_Parameters::gen_env_1 Number of generations spent in environment 1
     @param[in] HTE_Model_Parameters::HTE_Specific_Parameters::gen_env_2 Number of generations spent in environment 2
     @param[in, out] rng Random number generator
     @param[in] HTE_Model_Parameters::Fixed_Parameters::tolerance Tolerance for comparing equality of doubles
     @return
*/
  void run_simulation(double &allele_A_freq, const std::vector<double> &haploid_fitnesses,
		      const HTE_Model_Parameters &parameters, std::mt19937 &rng){
    // iterate through generations until one allele is fixed or the max number of generations is reached
    int gen = 0;
    int env_state = 0; // env 1 = 0; env 2 = 1
    while (gen < (parameters.model.gen_env_1 + parameters.model.gen_env_2) &&
	   !(closeToValue(allele_A_freq, 0.0, parameters.fixed.tolerance) ||
					 closeToValue(allele_A_freq, 1.0, parameters.fixed.tolerance))){
      expected_allele_freqs(allele_A_freq, haploid_fitnesses, env_state);
      realised_allele_freqs(allele_A_freq, parameters, rng);
      ++gen;
      if (gen == parameters.model.gen_env_1 - 1) {env_state = 1;}
    }
  }
  
}

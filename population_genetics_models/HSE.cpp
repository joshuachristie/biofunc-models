/**
   @file HSE.cpp
   @brief Contains functions to run the Haploid Single Environment model
*/

#include <cassert>
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include "Parameters.h"
#include "include/helper_functions.h"
#include "HSE.h"
#include "include/rng.h"
#include <iostream> // will probably need to be deleted once after refactoring and addition of print methods

/**
   @brief Namespace for Haploid Single Environment
*/
namespace HSE {
  /**
     @brief Reads in parameter values from command line into a struct
     @param[in] argc Number of commandline arguments
     @param[in] argv Command line arguments
     @return params HSE_Model_Parameters struct
  */
  HSE_Model_Parameters parse_parameter_values(int argc, char* argv[]){
    assert(std::string(argv[1]) == "HSE");
    assert(argc == 5 && "The HSE model must have 4 command line arguments (the first must be HSE)");
    const int population_size = atoi(argv[2]);
    const int number_generations = atoi(argv[3]);
    const double selection_coefficient = atof(argv[4]);
    HSE_Model_Parameters params {{population_size}, {number_generations, selection_coefficient}};
    return params;
  }
  /**
     @brief Calculates fitness function
     @param[in] parameters HSE_Model_Parameters::HSE_Specific_Parameters::selection_coefficient - selection coefficient of the A allele
     @return haploid_fitnesses A vector of length 2 containing the fitnesses of the A and a alleles [wA, wa]
  */
  std::vector<double> get_fitness_function(const HSE_Model_Parameters &parameters){
    // haploid_fitnesses[w_A, w_a] gives relative fitness of A allele
    std::vector<double> haploid_fitnesses {1.0 + parameters.model.selection_coefficient, 1.0};
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
     @param[in] HSE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @return Nothing (but modifies \p allele_A_freq)
*/
  void realised_allele_freqs(double &allele_A_freq, const HSE_Model_Parameters &parameters, std::mt19937 &rng){
    std::binomial_distribution<int> surviving_As(parameters.shared.population_size, allele_A_freq);
    allele_A_freq =
      static_cast<double>(surviving_As(rng)) / static_cast<double>(parameters.shared.population_size);
  }
  /**
     @brief Runs one replicate of the simulation
     @param[in] allele_A_freq Frequency of allele A
     @param[in] haploid_fitnesses A vector containing the fitnesses of the A and a alleles [wA, wa]
     @param[in] HSE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in] HSE_Model_Parameters::Shared_Parameters::number_generations Maximum number of generations for which the simulation will run (for the case in which a trait invades the ancestral population and is not challenged afterwards; for the case in which the fixed trait must withstand invaders, each step runs until an absorbing state is reached).
     @param[in, out] rng Random number generator
     @param[in] HSE_Model_Parameters::Fixed_Parameters::tolerance Tolerance for comparing equality of doubles
     @return
*/
  void run_simulation(double &allele_A_freq, const std::vector<double> &haploid_fitnesses,
		      const HSE_Model_Parameters &parameters, std::mt19937 &rng){
    // iterate through generations until one allele is fixed or the max number of generations is reached
    int gen = 0;
    while (gen < parameters.model.number_generations &&
	   !(closeToValue(allele_A_freq, 0.0, parameters.fixed.tolerance) ||
					 closeToValue(allele_A_freq, 1.0, parameters.fixed.tolerance))){
      expected_allele_freqs(allele_A_freq, haploid_fitnesses);
      realised_allele_freqs(allele_A_freq, parameters, rng);
      ++gen;
    }
  }
  /**
     @brief Runs Haploid Single Environment model
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @return Nothing (but prints results)
*/
  void run_model(int argc, char* argv[]){
    // should update initialiseRNG to snake case for consistency
    std::mt19937 rng = initialiseRNG();
    HSE_Model_Parameters params = parse_parameter_values(argc, argv);
    std::vector<double> haploid_fitnesses = get_fitness_function(params);
    std::vector<bool> final_A_freqs;
    // will change how I implement running of multiple replicates, but leaving it for post-refactor extension
    for (int rep = 0; rep < params.fixed.number_replicates; rep++){
      double allele_A_freq = 1.0 / static_cast<double>(params.shared.population_size); // initial freq is 1/N
      run_simulation(allele_A_freq, haploid_fitnesses, params, rng);
      // note that I also need to update closeToValue() (at the very least, change to snake_case; might also be worth hardcoding the tolerance given that it's now a fixed parameter)
      closeToValue(allele_A_freq, 0.0, params.fixed.tolerance) ? final_A_freqs.push_back(0) : final_A_freqs.push_back(1);
    }
    // will alter output in a later extension (will print to file directly from c++ rather than via the python run script)
    std::cout << std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) / static_cast<double>(params.fixed.number_replicates) << std::endl;

  }

}

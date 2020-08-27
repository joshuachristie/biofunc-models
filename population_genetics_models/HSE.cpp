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
#include "helper_functions.h"
#include "HSE.h"
#include "rng.h"
#include "persistence_probability.h"
#include "print_results.h"

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
    assert(argc == 4 && "The HSE model must have 3 command line arguments (the first must be HSE)");
    const int population_size = atoi(argv[2]);
    const double selection_coefficient = atof(argv[3]);
    const double initial_A_freq = 1.0 / static_cast<double>(population_size); // 1/N
    HSE_Model_Parameters params {{population_size}, {selection_coefficient, initial_A_freq}};
    return params;
  }
  /**
     @brief Calculates fitness function
     @param[in] parameters HSE_Model_Parameters::HSE_Specific_Parameters::selection_coefficient - selection coefficient of the A allele
     @return fitnesses A vector of length 2 containing the fitnesses of the A and a alleles [wA, wa]
  */
  std::vector<double> get_fitness_function(const HSE_Model_Parameters &parameters){
    // fitnesses[w_A, w_a] gives relative fitness of A allele
    std::vector<double> fitnesses {1.0 + parameters.model.selection_coefficient, 1.0};
    return fitnesses;
  }
  /**
     @brief Calculates frequency of the A allele after selection
     @param[in, out] allele_A_freq The frequency of the A allele
     @param[in] fitnesses A vector containing the fitnesses of the A and a alleles [wA, wa]
     @param[in] HSE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @return Nothing (but modifies \p allele_A_freq)
  */
  void calculate_allele_freqs(double &allele_A_freq, const std::vector<double> &haploid_fitnesses,
			      const HSE_Model_Parameters &parameters, std::mt19937 &rng){
    std::vector<double> expected_allele_freq_raw(2);
    expected_allele_freq_raw[0] = allele_A_freq * fitnesses[0];
    expected_allele_freq_raw[1] = (1.0 - allele_A_freq) * fitnesses[1];
    // get normalised expectation for allele_A_freq
    allele_A_freq = expected_allele_freq_raw[0] / (std::accumulate(expected_allele_freq_raw.begin(),
								   expected_allele_freq_raw.end(), 0.0));
    // sample to get realised outcome for allele_A_freq
    std::binomial_distribution<int> surviving_As(parameters.shared.population_size, allele_A_freq);
    allele_A_freq =
      static_cast<double>(surviving_As(rng)) / static_cast<double>(parameters.shared.population_size);
  }
  /**
     @brief Runs one replicate of the simulation
     @param[in] haploid_fitnesses A vector containing the fitnesses of the A and a alleles [wA, wa]
     @param[in] HSE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in] HSE_Model_Parameters::Fixed_Parameters::tolerance Tolerance for comparing equality of doubles
     @param[in, out] rng Random number generator
     @param[in, out] final_A_freqs Vector of bools storing whether allele A persists (true/false) at census
     @return Nothing (but alters \p final_A_freqs)
  */
  void run_simulation(const std::vector<double> &haploid_fitnesses, const HSE_Model_Parameters &parameters,
		      std::mt19937 &rng, std::vector<bool> &final_A_freqs){
    double allele_A_freq = 1.0 / static_cast<double>(parameters.shared.population_size); // initial freq is 1/N
    int gen = 0;
    while (help::is_neither_fixed_nor_extinct(gen, allele_A_freq, parameters)){
      calculate_allele_freqs(allele_A_freq, haploid_fitnesses, parameters, rng);
      ++gen;
    }
    help::close_to_value(allele_A_freq, 0.0, parameters.fixed.tolerance) ? final_A_freqs.push_back(0) :
      final_A_freqs.push_back(1);
  }
  /**
     @brief Runs Haploid Single Environment model
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @return Nothing (but prints results)
  */
  void run_model(int argc, char* argv[]){
    std::mt19937 rng = initialise_rng();
    HSE_Model_Parameters params = parse_parameter_values(argc, argv);
    std::vector<double> fitnesses = get_fitness_function(params);
    std::vector<bool> final_A_freqs;
    final_A_freqs.reserve(params.fixed.number_replicates);
    double persistence_probability = calculate_persistence_probability(params, run_simulation, rng,
								       haploid_fitnesses, final_A_freqs);
    print::print_persistence_probability(argc, argv, persistence_probability);
  }

} 

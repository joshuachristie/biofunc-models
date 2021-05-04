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
#include "trait_invasion.h"
#include "DataContainer.h"

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
    const double initial_trait_freq = 1.0 / static_cast<double>(population_size);
    const int number_reinvasions = atoi(argv[8]);
    const int number_gens_to_output_pp = atoi(argv[9]);
    const bool print_trait_raw_data = static_cast<bool>(atoi(argv[10]));
    const std::vector<int> trait_info {0, 1};
    const HTE_Model_Parameters params {{population_size, initial_trait_freq, number_reinvasions,
	number_gens_to_output_pp, print_trait_raw_data, trait_info}, {selection_coefficient_A_env_1,
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
     @brief Calculates frequency of the trait after selection
     @param[in, out] trait_freq The frequency of the trait (A allele)
     @param[in] fitnesses A vector containing the fitnesses of the A and a alleles [wA, wa]
     @param[in] HTE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @param[in, out] gen The current generation
     @return Nothing (but modifies \p trait_freq and increments \p gen)
  */
  void calculate_trait_freqs(std::vector<double> &trait_freq, const std::vector<double> &fitnesses,
			     const HTE_Model_Parameters &parameters, std::mt19937 &rng, int &gen){
    std::vector<double> expected_allele_freq_raw(2);
    expected_allele_freq_raw[0] = trait_freq[0] * fitnesses[gen >= parameters.model.gen_env_1];
    expected_allele_freq_raw[1] = (1.0 - trait_freq[0]) * fitnesses[2 + (gen >= parameters.model.gen_env_1)];
    // calculate normalised expectation of trait_freq
    double expectation = expected_allele_freq_raw[0] / (std::accumulate(expected_allele_freq_raw.begin(),
									expected_allele_freq_raw.end(), 0.0));
    // sample to get realised allele_A_freq
    std::binomial_distribution<int> surviving_As(parameters.shared.population_size, expectation);
    trait_freq[0] =
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
		       params.fixed.reserve_memory_trait_freq);
    calculate_persistence_probability(params, rng, fitnesses, calculate_trait_freqs, data);
    print::print_results(argc, argv, data, params);
  }
  
}

#include <cassert>
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include "Parameters.h"
#include "DSE.h"
#include "rng.h"
#include "conditional_existence_probability.h"
#include "trait_invasion.h"
#include "run_scenario.h"

namespace DSE {

  /**
     @brief Reads in parameter values from command line into a struct
     @param[in] argc Number of commandline arguments
     @param[in] argv Command line arguments
     @return params DSE_Model_Parameters struct
  */
  const DSE_Model_Parameters parse_parameter_values(int argc, char* argv[]){
    assert(std::string(argv[1]).compare("DSE"));
    assert(argc == 8 && "The DSE model must have 7 command line arguments (the first must be 'DSE')");
    assert((std::string(argv[2]).compare("LSTM") == 0 || std::string(argv[2]).compare("QEF") == 0) &&
	   "Incorrect model specification: specify whether the model type is LSTM or QEF in the second arg");
    const int population_size = atoi(argv[3]);
    const double selection_coefficient_homozygote = atof(argv[4]);
    const double selection_coefficient_heterozygote = atof(argv[5]);
    double initial_trait_freq = 1.0 / static_cast<double>(population_size);
    const int number_reinvasions = atoi(argv[6]);
    const std::vector<int> trait_info {atoi(argv[7]), 2};
    const DSE_Model_Parameters params {{population_size, initial_trait_freq, number_reinvasions, trait_info},
				       {selection_coefficient_homozygote, selection_coefficient_heterozygote}};
    return params;
  }
  /**
     @brief Calculates fitness function
     @param[in] parameters DSE_Model_Parameters::DSE_Specific_Parameters::selection_coefficient_homozygote - selection coefficient of the AA genotype
     @param[in] parameters DSE_Model_Parameters::DSE_Specific_Parameters::selection_coefficient_heterozygte - selection coefficient of the Aa genotype
     @return fitnesses A vector of length 3 containing the fitnesses of the AA, Aa, and aa genotypes [wAA, WAa, Waa]
  */
  const std::vector<double> get_fitness_function(const DSE_Model_Parameters &parameters){
    // diploid_fitnesses[w_AA, w_Aa, w_aa] gives relative fitness of genotypes
    const std::vector<double> fitnesses {1.0 + parameters.model.selection_coefficient_homozygote,
      1.0 + parameters.model.selection_coefficient_heterozygote, 1.0};
    return fitnesses;
  }
  /**
     @brief Calculates frequency of the trait after selection and random mating
     @param[in, out] trait_freq The frequency of the trait
     @param[in] fitnesses Vector containing AA, Aa, and aa genotype fitnesses [wAA, wAa, waa]
     @param[in] DSE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @param[in, out] gen The current generation
     @return Nothing (but modifies \p trait_freq and increments \p gen)
  */
  void calculate_trait_freqs(std::vector<double> &trait_freq, const std::vector<double> &fitnesses,
			     const DSE_Model_Parameters &parameters, std::mt19937 &rng, int &gen){
    double allele_A_freq = trait_freq[0] + 0.5 * trait_freq[1];
    std::vector<double> expected_genotype_freq(3);
    expected_genotype_freq[0] = std::pow(allele_A_freq, 2.0) * fitnesses[0]; // AA
    expected_genotype_freq[1] = 2 * allele_A_freq * (1.0 - allele_A_freq) * fitnesses[1]; // Aa
    expected_genotype_freq[2] = std::pow((1 - allele_A_freq), 2.0) * fitnesses[2]; // aa
    // multinomial sample to get realised outcome for trait_freq (discrete_distribution normalises probs)
    std::discrete_distribution<int> multinom {expected_genotype_freq.begin(), expected_genotype_freq.end()};
    // sample surviving (individuals with) traits
    std::vector<int> surviving_traits(3);
    for (int i = 0; i < parameters.shared.population_size; i++){
      ++surviving_traits[ multinom(rng) ];
    }
    // convert counts into proportions for individuals with traits AA and Aa
    trait_freq[0] = static_cast<double>(surviving_traits[0]) / parameters.shared.population_size;
    trait_freq[1] = static_cast<double>(surviving_traits[1]) / parameters.shared.population_size;
    ++gen;
  }

  /**
     @brief Runs Diploid Single Environment model
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @return Nothing (but prints results)
  */
  void run_model(int argc, char* argv[]){
    std::mt19937 rng = initialise_rng();
    const DSE_Model_Parameters params = parse_parameter_values(argc, argv);
    const std::vector<double> fitnesses = get_fitness_function(params);

    if (std::string(argv[2]).compare("QEF") == 0){
      run_scenario::QEF(params, rng, fitnesses, calculate_trait_freqs, argv);
    } else if (std::string(argv[2]).compare("LSTM") == 0){
      run_scenario::LSTM(params, rng, fitnesses, calculate_trait_freqs, argv);
    }

  }

}

/**
   @file DSE.cpp
   @brief Contains functions to run the Diploid Single Environment model
*/

#include <cassert>
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include "Parameters.h"
#include "helper_functions.h"
#include "DSE.h"
#include "rng.h"
#include <iostream> // will probably need to be deleted once after refactoring and addition of print methods

/**
   @brief Namespace for Diploid Single Environment
*/
namespace DSE {

    /**
     @brief Reads in parameter values from command line into a struct
     @param[in] argc Number of commandline arguments
     @param[in] argv Command line arguments
     @return params DSE_Model_Parameters struct
  */
  DSE_Model_Parameters parse_parameter_values(int argc, char* argv[]){
    assert(std::string(argv[1]) == "DSE");
    assert(argc == 6 && "The DSE model must have 5 command line arguments (the first must be 'DSE')");
    const int population_size = atoi(argv[2]);
    const int number_generations = atoi(argv[3]);
    const double selection_coefficient_homozygote = atof(argv[4]);
    const double selection_coefficient_heterozygote = atof(argv[5]);
    DSE_Model_Parameters params {{population_size}, {number_generations, selection_coefficient_homozygote,
						     selection_coefficient_heterozygote}};
    return params;
  }
  /**
     @brief Calculates fitness function
     @param[in] parameters DSE_Model_Parameters::DSE_Specific_Parameters::selection_coefficient_homozygote - selection coefficient of the AA genotype
     @param[in] parameters DSE_Model_Parameters::DSE_Specific_Parameters::selection_coefficient_heterozygte - selection coefficient of the Aa genotype
     @return genotype_fitnesses A vector of length 3 containing the fitnesses of the AA, Aa, and aa genotypes [wAA, WAa, Waa]
  */
  std::vector<double> get_fitness_function(const DSE_Model_Parameters &parameters){
    // diploid_fitnesses[w_AA, w_Aa, w_aa] gives relative fitness of genotypes
    std::vector<double> genotype_fitnesses {1.0 + parameters.model.selection_coefficient_homozygote,
					    1.0 + parameters.model.selection_coefficient_heterozygote, 1.0};
    return genotype_fitnesses;
  }
  /**
     @brief Calculates expected frequency of the A allele after selection
     @param[in, out] allele_A_freq The frequency of the A allele
     @param[in] genotype_fitnesses Vector containing AA, Aa, and aa genotype fitnesses [wA, wAa, wa]
     @return Nothing (but modifies \p allele_A_freq)
*/
  void expected_allele_freqs(double &allele_A_freq, const std::vector<double> &genotype_fitnesses){
    std::vector<double> expected_genotype_freq_raw(3);
    expected_genotype_freq_raw[0] = std::pow(allele_A_freq, 2.0) * genotype_fitnesses[0]; // AA
    expected_genotype_freq_raw[1] = 2 * allele_A_freq * (1.0 - allele_A_freq) * genotype_fitnesses[1]; // Aa
    expected_genotype_freq_raw[2] = std::pow((1 - allele_A_freq), 2.0) * genotype_fitnesses[2]; // aa
    // replace allele_A_freq with normalised expectation
    allele_A_freq = (expected_genotype_freq_raw[0] + 0.5 * expected_genotype_freq_raw[1]) /		(std::accumulate(expected_genotype_freq_raw.begin(), expected_genotype_freq_raw.end(), 0.0));
  }
  /**
     @brief Calculates realised frequency of the A allele (by binomial sampling of expectated number of As)
     @param[in, out] allele_A_freq Frequency of allele A
     @param[in] DSE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @return Nothing (but modifies \p allele_A_freq)
  */
  void realised_allele_freqs(double &allele_A_freq, const DSE_Model_Parameters &parameters, std::mt19937 &rng){
    // population_size multiplied by 2 because organisms are diploid
    std::binomial_distribution<int> surviving_As(parameters.shared.population_size * 2, allele_A_freq);
    allele_A_freq =
      static_cast<double>(surviving_As(rng)) / static_cast<double>(parameters.shared.population_size * 2);
  }
  /**
     @brief Runs one replicate of the simulation
     @param[in] allele_A_freq Frequency of allele A
     @param[in] @param[in] genotype_fitnesses Vector containing AA, Aa, and aa genotype fitnesses [wA, wAa, wa]
     @param[in] DSE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in] DSE_Model_Parameters::Shared_Parameters::number_generations Maximum number of generations for which the simulation will run (for the case in which a trait invades the ancestral population and is not challenged afterwards; for the case in which the fixed trait must withstand invaders, each step runs until an absorbing state is reached).
     @param[in, out] rng Random number generator
     @param[in] DSE_Model_Parameters::Fixed_Parameters::tolerance Tolerance for comparing equality of doubles
     @return tbd (nothing right now and will probably remain void but might call a print function from here)
*/
  void run_simulation(double &allele_A_freq, const std::vector<double> &genotype_fitnesses,
		      const DSE_Model_Parameters &parameters, std::mt19937 &rng){
    // iterate through generations until one allele is fixed or the max number of generations is reached
    int gen = 0;
    while (gen < parameters.model.number_generations &&
	   !(close_to_value(allele_A_freq, 0.0, parameters.fixed.tolerance) ||
					 close_to_value(allele_A_freq, 1.0, parameters.fixed.tolerance))){
      expected_allele_freqs(allele_A_freq, genotype_fitnesses);
      realised_allele_freqs(allele_A_freq, parameters, rng);
      ++gen;
    }
  }
  /**
     @brief Runs Diploid Single Environment model
     @param[in] argc Number of command line arguments
     @param[in] argv Array of command line arguments
     @return Nothing (but prints results)
*/
  void run_model(int argc, char* argv[]){
    std::mt19937 rng = initialise_rng();
    DSE_Model_Parameters params = parse_parameter_values(argc, argv);
    std::vector<double> genotype_fitnesses = get_fitness_function(params);
    std::vector<bool> final_A_freqs;

    // will change how I implement running of multiple replicates, but leaving it for post-refactor extension
    for (int rep = 0; rep < params.fixed.number_replicates; rep++){
      double allele_A_freq = 1.0 / static_cast<double>(params.shared.population_size * 2); // initial freq is 1/2N
      run_simulation(allele_A_freq, genotype_fitnesses, params, rng);
      // note that I also need to update close_to_value() (at the very least, change to snake_case; might also be worth hardcoding the tolerance given that it's now a fixed parameter)
      close_to_value(allele_A_freq, 0.0, params.fixed.tolerance) ? final_A_freqs.push_back(0) : final_A_freqs.push_back(1);
    }
    // will alter output in a later extension (will print to file directly from c++ rather than via the python run script)
    std::cout << std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) / static_cast<double>(params.fixed.number_replicates) << std::endl;

  }
}

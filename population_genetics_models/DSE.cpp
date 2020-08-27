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
#include "persistence_probability.h"
#include "print_results.h"

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
    assert(argc == 5 && "The DSE model must have 4 command line arguments (the first must be 'DSE')");
    const int population_size = atoi(argv[2]);
    const double selection_coefficient_homozygote = atof(argv[3]);
    const double selection_coefficient_heterozygote = atof(argv[4]);
    double allele_A_freq = 1.0 / static_cast<double>(population_size * 2); // 1/2N
    DSE_Model_Parameters params {{population_size}, {selection_coefficient_homozygote,
	selection_coefficient_heterozygote, allele_A_freq}};
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
     @brief Calculates frequency of the A allele after selection and random mating
     @param[in, out] allele_A_freq The frequency of the A allele
     @param[in] genotype_fitnesses Vector containing AA, Aa, and aa genotype fitnesses [wA, wAa, wa]
     @param[in] DSE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in, out] rng Random number generator
     @return Nothing (but modifies \p allele_A_freq)
  */
  void calculate_allele_freqs(double &allele_A_freq, const std::vector<double> &genotype_fitnesses,
			      const DSE_Model_Parameters &parameters, std::mt19937 &rng){
    std::vector<double> expected_genotype_freq_raw(3);
    expected_genotype_freq_raw[0] = std::pow(allele_A_freq, 2.0) * genotype_fitnesses[0]; // AA
    expected_genotype_freq_raw[1] = 2 * allele_A_freq * (1.0 - allele_A_freq) * genotype_fitnesses[1]; // Aa
    expected_genotype_freq_raw[2] = std::pow((1 - allele_A_freq), 2.0) * genotype_fitnesses[2]; // aa
    // get normalised expectation for allele_A_freq
    allele_A_freq = (expected_genotype_freq_raw[0] + 0.5 * expected_genotype_freq_raw[1]) /
      (std::accumulate(expected_genotype_freq_raw.begin(), expected_genotype_freq_raw.end(), 0.0));
    // sample to get realised outcome for allele_A_freq (diploid, so population_size multiplied by 2)
    std::binomial_distribution<int> surviving_As(parameters.shared.population_size * 2, allele_A_freq);
    allele_A_freq =
      static_cast<double>(surviving_As(rng)) / static_cast<double>(parameters.shared.population_size * 2);
  }
  /**
     @brief Runs one replicate of the simulation
     @param[in] genotype_fitnesses Vector containing AA, Aa, and aa genotype fitnesses [wA, wAa, wa]
     @param[in] DSE_Model_Parameters::Shared_Parameters::population_size Number of individuals in the population
     @param[in] DSE_Model_Parameters::Fixed_Parameters::tolerance Tolerance for comparing equality of doubles
     @param[in, out] rng Random number generator
     @param[in, out] final_A_freqs Vector of bools storing whether allele A persists (true/false) at census
     @return Nothing (but alters \p final_A_freqs)
  */
  void run_simulation(const std::vector<double> &genotype_fitnesses, const DSE_Model_Parameters &parameters,
		      std::mt19937 &rng, std::vector<bool> &final_A_freqs){
    double allele_A_freq = 1.0 / static_cast<double>(parameters.shared.population_size * 2); // initial freq is 1/2N
    int gen = 0;
    while (help::is_neither_fixed_nor_extinct(gen, allele_A_freq, parameters)){
      calculate_allele_freqs(allele_A_freq, genotype_fitnesses, parameters, rng);
      ++gen;
    }
    help::close_to_value(allele_A_freq, 0.0, parameters.fixed.tolerance) ? final_A_freqs.push_back(0) :
      final_A_freqs.push_back(1);
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
    final_A_freqs.reserve(params.fixed.number_replicates);
    double persistence_probability = calculate_persistence_probability(params, run_simulation, rng,
								       genotype_fitnesses, final_A_freqs);
    print::print_persistence_probability(argc, argv, persistence_probability);
  }


}

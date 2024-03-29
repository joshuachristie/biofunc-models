#include <cassert>
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include "Parameters.h"
#include "HSE.h"
#include "rng.h"
#include "conditional_existence_probability.h"
#include "trait_invasion.h"
#include "run_scenario.h"

namespace HSE {
  
  const parameters::HSE_Model_Parameters parse_parameter_values(int argc, char* argv[]){
    assert(std::string(argv[1]).compare("HSE") == 0);
    assert(argc == 6 && "The HSE model must have 5 command line arguments (the first must be HSE)");
    assert((std::string(argv[2]).compare("LSTM") == 0 || std::string(argv[2]).compare("QEF") == 0) &&
	   "Incorrect model specification: specify whether the model type is LSTM or QEF in the second arg");
    const int population_size = atoi(argv[3]);
    const double selection_coefficient = atof(argv[4]);
    const double initial_trait_frequency = 1.0 / static_cast<double>(population_size);
    const int number_reinvasions = atoi(argv[5]);
    const std::vector<int> trait_info {0, 1};
    const parameters::HSE_Model_Parameters params {{population_size, initial_trait_frequency,
	number_reinvasions, trait_info}, {selection_coefficient}};
    return params;
  }
  /**
     @details Calculates and returns a 1-by-2 vector of allele fitnesses [w_A, w_a].
  */
  const std::vector<double> get_fitness_function(const parameters::HSE_Model_Parameters &parameters){
    const std::vector<double> fitnesses {1.0 + parameters.model.selection_coefficient, 1.0};
    return fitnesses;
  }
  /**
     @details The function first calculates the expected (deterministic) allele frequency due to selection.
     It then uses this expectation as the probability for a (random) binomial sampling process to get a new
     \p trait_freq. It also increments the current \p gen.
  */
  void calculate_trait_freqs(std::vector<double> &trait_freq, const std::vector<double> &fitnesses,
			     const parameters::HSE_Model_Parameters &parameters, std::mt19937 &rng, int &gen){
    std::vector<double> expected_allele_freq_raw(2);
    expected_allele_freq_raw[0] = trait_freq[0] * fitnesses[0];
    expected_allele_freq_raw[1] = (1.0 - trait_freq[0]) * fitnesses[1];
    // get normalised expectation for trait_freq
    double expectation = expected_allele_freq_raw[0] / (std::accumulate(expected_allele_freq_raw.begin(),
									expected_allele_freq_raw.end(), 0.0));
    // sample to get realised outcome for trait_freq
    std::binomial_distribution<int> surviving_As(parameters.shared.population_size, expectation);
    trait_freq[0] =
      static_cast<double>(surviving_As(rng)) / static_cast<double>(parameters.shared.population_size);
    ++gen;
  }
  /**
     @details Calls initialise_rng(), HSE::parse_parameter_values(), HSE::get_fitness_function(), calls calculate_conditional_existence_probability(), and finally calls
     print::print_results().
  */
  void run_model(int argc, char* argv[]){

    std::mt19937 rng = rng::initialise_rng();
    const parameters::HSE_Model_Parameters params = parse_parameter_values(argc, argv);
    const std::vector<double> fitnesses = get_fitness_function(params);

    if (std::string(argv[2]).compare("QEF") == 0){
      run_scenario::QEF(params, rng, fitnesses, calculate_trait_freqs, argv, argc);
    } else if (std::string(argv[2]).compare("LSTM") == 0){
      run_scenario::LSTM(params, rng, fitnesses, calculate_trait_freqs, argv, argc);
    }
    
  }

}

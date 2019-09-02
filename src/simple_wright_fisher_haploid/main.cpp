#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include "functions.h"
int main(int argc, char* argv[]){
  // parameters
  const double selection_coefficient = atof(argv[1]);
  const int population_size = atoi(argv[2]);
  const int number_generations = atoi(argv[3]);
  const int number_replicates = atoi(argv[4]);
  // set rng
  std::mt19937 temp_rng(std::random_device{}());
  std::uniform_int_distribution<> adjust_seed(0, 5000);
  int factor_to_adjust_seed = adjust_seed(temp_rng);
  auto seed =
    (std::chrono::high_resolution_clock::now().time_since_epoch().count()) * factor_to_adjust_seed;
  std::mt19937 rng(seed);
  const std::vector<double> haploid_fitnesses = getFitnessFunction(selection_coefficient);
  std::vector<bool> final_A_freqs;
  // iterate over generations and replicates
  for (int rep = 0; rep < number_replicates; rep++){
    double allele_A_freq = 1.0 / static_cast<double>(population_size); // initial freq is 1/N
    iterateOverGenerations(allele_A_freq, haploid_fitnesses, population_size, number_generations, rng);
    // record 0 if extinct, 1 otherwise (indicating the allele exists at a non-zero proportion)
    closeToValue(allele_A_freq, 0.0) ? final_A_freqs.push_back(0) : final_A_freqs.push_back(1);
  }
  std::cout << std::accumulate(final_A_freqs.begin(), final_A_freqs.end(), 0.0) / static_cast<double>(number_replicates) << std::endl;
  return 0;
}

#include <vector>
#include <random>
#include "functions.h"
int main(int argc, char* argv[]){
  // parameters
  const double selection_coefficient = atof(argv[1]);
  const int population_size = atoi(argv[2]);
  const int number_generations = atoi(argv[3]);
  // set rng
  std::mt19937 temp_rng(std::random_device{}());
  std::uniform_int_distribution<> adjust_seed(0, 5000);
  int factor_to_adjust_seed = adjust_seed(temp_rng);
  auto seed =
    (std::chrono::high_resolution_clock::now().time_since_epoch().count()) * factor_to_adjust_seed;
  std::mt19937 rng(seed);

  
  const std::vector<double> haploid_fitnesses = getFitnessFunction(selection_coefficient);
  double allele_A_freq = 1.0 / static_cast<double>(number_generations); // initial freq is 1/N

  expectedAlleleFreqs(allele_A_freq, haploid_fitnesses);
  realisedAlleleFreqs(allele_A_freq, population_size);
  return 0;
}
  
  

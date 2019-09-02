#include <vector>
#include "functions.h"
int main(int argc, char* argv[]){
  // parameters
  const double selection_coefficient = atof(argv[1]);
  const int population_size = atoi(argv[2]);
  const int number_generations = atoi(argv[3]);
  
  const std::vector<double> haploid_fitnesses = getFitnessFunction(selection_coefficient);
  double allele_A_freq = 1.0 / static_cast<double>(number_generations); // initial freq is 1/N

  

  return 0;
}
  
  

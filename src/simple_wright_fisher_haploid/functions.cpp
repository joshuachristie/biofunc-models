#include <cmath>
#include "functions.h"

std::vector<double> getFitnessFunction(const double selection_coefficient){
  // haploid_fitnesses[w_A, w_a] gives relative fitness of A allele
  std::vector<double> haploid_fitnesses {1.0 + selection_coefficient, 1.0};
  return haploid_fitnesses;
}

void expectedAlleleFreqs(double &allele_A_freq, const std::vector<double> haploid_fitnesses){
  // expected frequency of allele A after selection
  std::vector<double> expected_allele_freq_raw(2);
  expected_allele_freq_raw[0] = allele_A_freq * haploid_fitnesses[0];
  expected_allele_freq_raw[1] = (1.0 - allele_A_freq) * haploid_fitnesses[1];
  // replace allele_A_freq with normalised expectation
  allele_A_freq = expected_allele_freq_raw[0] / (std::accumulate(expected_allele_freq_raw.begin(),
								 expected_allele_freq_raw.end(), 0.0));
}

void realisedAlleleFreqs(double &allele_A_freq, const int population_size, std::mt19937 &rng){
  // realised frequency of allele A after selection
  std::binomial_distribution<int> surviving_As(population_size, allele_A_freq); // use expectation of allele A
  allele_A_freq = static_cast<double>(surviving_As(rng)) / static_cast<double>(population_size); // proportion
}

void iterateOverGenerations(double &allele_A_freq, const std::vector<double> haploid_fitnesses,
			    const int population_size, const int number_generations, std::mt19937 &rng){
  for (int gen = 0; gen < number_generations; gen++){
    expectedAlleleFreqs(allele_A_freq, haploid_fitnesses);
    realisedAlleleFreqs(allele_A_freq, population_size, rng);
  }
}

bool closeToValue(double allele_freq, double value){
  if (std::abs(allele_freq - value) < 0.00000001){
    return 1;
  } else {
    return 0;
  }
}

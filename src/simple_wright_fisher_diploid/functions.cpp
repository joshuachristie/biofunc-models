#include <cmath>
#include <random>
#include <numeric>
#include "functions.h"

std::vector<double> getFitnessFunction(const double selection_coefficient_homozygote, const double selection_coefficient_heterozygote){
  // diploid_fitnesses[w_AA, w_Aa, w_aa] gives relative fitness of genotypes
  // std::vector<double> genotype_fitnesses {1.0 + selection_coefficient * 2,
  // 					  1.0 + selection_coefficient * 2 * dominance_coefficient, 1.0};
  std::vector<double> genotype_fitnesses {1.0 + selection_coefficient_homozygote,
					  1.0 + selection_coefficient_heterozygote, 1.0};

  return genotype_fitnesses;
}

void expectedAlleleFreqs(double &allele_A_freq, const std::vector<double> genotype_fitnesses){
  // expected frequency of the A allele after selection
  std::vector<double> expected_genotype_freq_raw(3);
  expected_genotype_freq_raw[0] = std::pow(allele_A_freq, 2.0) * genotype_fitnesses[0]; // AA
  expected_genotype_freq_raw[1] = 2 * allele_A_freq * (1.0 - allele_A_freq) * genotype_fitnesses[1]; // Aa
  expected_genotype_freq_raw[2] = std::pow((1 - allele_A_freq), 2.0) * genotype_fitnesses[2]; // aa
  // replace allele_A_freq with normalised expectation
  allele_A_freq = (expected_genotype_freq_raw[0] + 0.5 * expected_genotype_freq_raw[1]) /
    (std::accumulate(expected_genotype_freq_raw.begin(), expected_genotype_freq_raw.end(), 0.0));
}

void realisedAlleleFreqs(double &allele_A_freq, const int population_size, std::mt19937 &rng){
  // realised frequency of allele A after selection
  std::binomial_distribution<int> surviving_As(population_size * 2, allele_A_freq);
  allele_A_freq = static_cast<double>(surviving_As(rng)) / static_cast<double>(population_size * 2);
}

void iterateOverGenerations(double &allele_A_freq, const std::vector<double> genotype_fitnesses,
			    const int population_size, const int number_generations, std::mt19937 &rng){
  // iterate through generations until one allele is fixed or the max number of generations is reached
  int gen = 0;
  while (gen < number_generations && !(closeToValue(allele_A_freq, 0.0) || closeToValue(allele_A_freq, 1.0))){
    expectedAlleleFreqs(allele_A_freq, genotype_fitnesses);
    realisedAlleleFreqs(allele_A_freq, population_size, rng);
    ++gen;
  }
}

bool closeToValue(double allele_freq, double value){
  // returns true if allele_freq is approx equal to value
  if (std::abs(allele_freq - value) < 0.00000001){
    return 1;
  } else {
    return 0;
  }
}

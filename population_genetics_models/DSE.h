#ifndef DSE_H
#define DSE_H

namespace DSE {
  
  DSE_Model_Parameters parse_parameter_values(int argc, char* argv[]);
  
  std::vector<double> get_fitness_function(const DSE_Model_Parameters &parameters);
  
  void expected_allele_freqs(double &allele_A_freq, const std::vector<double> &genotype_fitnesses);
  
  void realised_allele_freqs(double &allele_A_freq, const DSE_Model_Parameters &parameters, std::mt19937 &rng);
  
  void run_simulation(double &allele_A_freq, const std::vector<double> &genotype_fitnesses,
		      const DSE_Model_Parameters &parameters, std::mt19937 &rng);
  
}

#endif

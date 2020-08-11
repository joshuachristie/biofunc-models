#ifndef HSE_H
#define HSE_H

namespace HSE{
HSE_Model_Parameters parse_parameter_values(int argc, char* argv[]);
 
std::vector<double> get_fitness_function(const HSE_Model_Parameters &parameters);
 
void expected_allele_freqs(double &allele_A_freq, const std::vector<double> &haploid_fitnesses);
 
void realised_allele_freqs(double &allele_A_freq, const HSE_Model_Parameters &parameters, std::mt19937 &rng);
 
void run_simulation(double &allele_A_freq, const std::vector<double> &haploid_fitnesses,
		    const HSE_Model_Parameters &parameters, std::mt19937 &rng);
}

#endif

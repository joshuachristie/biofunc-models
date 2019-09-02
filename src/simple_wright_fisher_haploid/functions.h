std::vector<double> getFitnessFunction(const double selection_coefficient);
void expectedAlleleFreqs(double &allele_A_freq, std::vector<double> &haploid_fitnesses);
void realisedAlleleFreqs(double &allele_A_freq, const int population_size, std::mt19937 &rng);


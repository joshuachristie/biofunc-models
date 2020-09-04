#ifndef DATA_H
#define DATA_H

#include <vector>

class DataPersistenceInfiniteApprox {
private:
  bool _persistence;
public:
  void set_infinite_approx_persistence(bool);
  void print_infinite_approx_persistence();
  /** constructors */
  DataPersistenceInfiniteApprox() {}
  DataPersistenceInfiniteApprox(bool persistence) : _persistence{persistence} {}
};

class DataPersistenceMultipleGens : public DataPersistenceInfiniteApprox {
private:
  std::vector<bool> _persistence_by_gen;
public:
  void append_persistence_by_gen(bool);
  void print_persistence_by_gen();
  /** constructors */
  DataPersistenceMultipleGens() = default;
  DataPersistenceMultipleGens(int reserve_length_mg) { _persistence_by_gen.reserve(reserve_length_mg); }
};

class DataAlleleFreq : public DataPersistenceMultipleGens {
private:
  std::vector<double> _allele_A_freq_by_gen;
public:
  void append_allele_A_freq_by_gen(double);
  void print_allele_A_freq_by_gen();
  /** constructors */
  DataAlleleFreq(int reserve_length_af) { _allele_A_freq_by_gen.reserve(reserve_length_af); }
  DataAlleleFreq(int reserve_length_mg, int reserve_length_af) :
    DataPersistenceMultipleGens {reserve_length_mg} { _allele_A_freq_by_gen.reserve(reserve_length_af); }
};

#endif

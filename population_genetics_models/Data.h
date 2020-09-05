#ifndef DATA_H
#define DATA_H

#include <vector>

class DataPersistenceInfinite {
public:
  bool _persistence;
  /** constructors */
  DataPersistenceInfinite() = default;
  DataPersistenceInfinite(bool persistence) : _persistence{persistence} {}
};

class DataPersistenceByGen : public DataPersistenceInfinite {
public:
  std::vector<bool> _persistence_by_gen;
  void append_persistence_by_gen(bool);
  /** constructors */
  DataPersistenceByGen() = default;
  DataPersistenceByGen(int reserve_length_mg) { _persistence_by_gen.reserve(reserve_length_mg); }
};

class DataAlleleFreq : public DataPersistenceByGen {
public:
  std::vector<double> _allele_A_freq_by_gen;
  void append_allele_A_freq_by_gen(double);
  /** constructors */
  DataAlleleFreq(int reserve_length_af) { _allele_A_freq_by_gen.reserve(reserve_length_af); }
  DataAlleleFreq(int reserve_length_mg, int reserve_length_af) :
    DataPersistenceByGen {reserve_length_mg} { _allele_A_freq_by_gen.reserve(reserve_length_af); }
};

#endif

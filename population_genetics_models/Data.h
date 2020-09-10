#ifndef DATA_H
#define DATA_H

#include <vector>

class Data {
public:
  bool _persistence;
  std::vector<bool> _persistence_by_gen;
  std::vector<double> _allele_A_freq_by_gen;

  void append_persistence_by_gen(bool);
  void append_allele_A_freq_by_gen(double);
  
  Data(int number_gens_to_record_pp, int reserve_length_af) {
    _persistence_by_gen.reserve(number_gens_to_record_pp);
    _allele_A_freq_by_gen.reserve(reserve_length_af);
  }
    
};
  

// class DataPersistenceInfinite {
// public:
//   bool _persistence;
//   /** constructors */
//   DataPersistenceInfinite() = default;
//   DataPersistenceInfinite(bool persistence) : _persistence{persistence} {}
// };

// class DataPersistenceByGen : public DataPersistenceInfinite {
// public:
//   std::vector<bool> _persistence_by_gen;
//   void append_persistence_by_gen(bool);
//   /** constructor */
//   DataPersistenceByGen(int reserve_length_mg) { _persistence_by_gen.reserve(reserve_length_mg); }
// };

// class DataAlleleFreq : public DataPersistenceInfinite {
// public:
//   std::vector<double> _allele_A_freq_by_gen;
//   void append_allele_A_freq_by_gen(double);
//   /** constructor */
//   DataAlleleFreq(int reserve_length_af) { _allele_A_freq_by_gen.reserve(reserve_length_af); }
// };

// class DataAlleleFreqAndPBG : public DataPersistenceByGen {
// public:
//   std::vector<double> _allele_A_freq_by_gen;
//   void append_allele_A_freq_by_gen(double);
//   /** constructor */
//   DataAlleleFreqAndPBG(int reserve_length_mg, int reserve_length_af) :
//     DataPersistenceByGen {reserve_length_mg} { _allele_A_freq_by_gen.reserve(reserve_length_af); }
// };

#endif

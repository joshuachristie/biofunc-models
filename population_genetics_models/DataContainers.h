#ifndef DATA_CONTAINERS_H
#define DATA_CONTAINERS_H

#include <vector>
#include "Data.h"

// keep in mind that it might be possible to make it cleaner by inheriting from templated DataContainer
// but it might only be the one line with _simulation_data so maybe more work than worth
// but whether I inherit from a base template class, i will want to be doing some inheritance from the
// specialised classes as they build on each other (e.g. I need to be able to print the infinite approx fixation
// from all three classes)


// also just a note that it is probably worth having some generic print_results() member function for all of the
// classes but for it to be polymorphic (for DataPersistenceInfiniteApprox it just prints the infinite approx, for DataPersistenceMultipleGens it prints the infinite approx (which will be inherited) and the multiple gen print (which will be a member function of DataPersistenceMultipleGens; for DataAlleleFreq it prints the infinite approx (inherited from DataPersistenceInfiniteApprox), multiple gen (inherited from DataPersistenceMultipleGens), and the raw data (member function of DataAlleleFreq)

// actually changed my mind on this - I don't think I should be printing from within the member functions at all
// I can keep some (simplified) free print functions but the member functions here should just be returning the
// data 

template <class, int number_replicates>
class DataContainer {};

template <int number_replicates>
class DataContainer <DataPersistenceInfinite, number_replicates> {
// private:
//   std::vector<DataPersistenceInfinite> _simulation_data;
public:
  std::vector<DataPersistenceInfinite> _simulation_data;
  DataContainer(){
    for (int i = 0; i < number_replicates; i++){
      _simulation_data.push_back(DataPersistenceInfinite());
    }
  }
  const double get_persistence_infinite_approx(){
    int persistence_count = 0;
    for (auto it = _simulation_data.begin(); it != _simulation_data.end(); it++){
      persistence_count += it->_persistence;
    }
    return static_cast<double>(persistence_count) / static_cast<double>(_simulation_data.size());
  };

};

template <int number_replicates>
class DataContainer <DataPersistenceByGen, number_replicates> {
// private:
//   std::vector<DataPersistenceByGen> _simulation_data;
public:
  std::vector<DataPersistenceByGen> _simulation_data;
  DataContainer(int reserve_length_mg){
    for (int i = 0; i < number_replicates; i++){
      _simulation_data.push_back(DataPersistenceByGen(reserve_length_mg));
    }
  }
  
  const std::vector<double> get_persistence_by_gen(){
    
    std::vector<double> persistence_probability(_simulation_data.size(), 0.0);
    for (int i = 0; i < number_replicates; i++){
      for (int j = 0; j < persistence_probability.size(); j++){
	persistence_probability[j] += (static_cast<double>(_simulation_data[i]._persistence_by_gen[j]) /
				       static_cast<double>(number_replicates));
      }
    }
    return persistence_probability;
  }
  
};

template <int number_replicates>
class DataContainer <DataAlleleFreq, number_replicates> {
// private:
//   std::vector<DataAlleleFreq> _simulation_data;
public:
  std::vector<DataAlleleFreq> _simulation_data;
  DataContainer(int reserve_length_mg, int reserve_length_af){
    for (int i = 0; i < number_replicates; i++){
      _simulation_data.push_back(DataAlleleFreq(reserve_length_mg, reserve_length_af));
    }
  }
  DataContainer(int reserve_length_af){
    for (int i = 0; i < number_replicates; i++){
      _simulation_data.push_back(DataAlleleFreq(reserve_length_af));
    }
  }

  double get_allele_A_freqs(int replicate, int gen){
    return _simulation_data[replicate]._allele_A_freq_by_gen[gen];
  }
};
  
#endif

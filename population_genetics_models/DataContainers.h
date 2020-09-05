#ifndef DATA_CONTAINERS_H
#define DATA_CONTAINERS_H

#include <vector>
#include "Data.h"

template <class T>
class DataContainer_Base_Infinite_Approx {
public:
  
  std::vector<T> _simulation_data;
  
  virtual const double get_persistence_infinite_approx(){
    int persistence_count = 0;
    for (auto it = _simulation_data.begin(); it != _simulation_data.end(); it++){
      persistence_count += it->_persistence;
    }
    return static_cast<double>(persistence_count) / static_cast<double>(_simulation_data.size());
  };

};

template <class T, int number_replicates>
class DataContainer_Base_Persistence_By_Gen : public virtual DataContainer_Base_Infinite_Approx<T>
{
public:
  
  virtual const std::vector<double> get_persistence_by_gen(){
    std::vector<double> persistence_probability(this->_simulation_data.size(), 0.0);
    for (int i = 0; i < number_replicates; i++){
      for (int j = 0; j < persistence_probability.size(); j++){
	persistence_probability[j] += (static_cast<double>(this->_simulation_data[i]._persistence_by_gen[j]) /
				       static_cast<double>(number_replicates));
      }
    }
    return persistence_probability;
  }
 
};

template <class T>
class DataContainer_Base_Allele_A_Freq : public virtual  DataContainer_Base_Infinite_Approx<T> {
public:
  const double get_allele_A_freqs(int replicate, int gen){
    return this->_simulation_data[replicate]._allele_A_freq_by_gen[gen];
  }
};



template <class T, int number_replicates>
class DataContainer : public DataContainer_Base_Infinite_Approx<T> {};

template <int number_replicates>
class DataContainer <DataPersistenceInfinite, number_replicates> : public
DataContainer_Base_Infinite_Approx<DataPersistenceInfinite> {
public:
  // constructor
  DataContainer(){
    for (int i = 0; i < number_replicates; i++){
      this->_simulation_data.push_back(DataPersistenceInfinite());
    }
  }

};

template <int number_replicates>
class DataContainer <DataPersistenceByGen, number_replicates> :
  public DataContainer_Base_Persistence_By_Gen<DataPersistenceByGen, number_replicates> {
public:
  // constructor
  DataContainer(int reserve_length_mg){
    for (int i = 0; i < number_replicates; i++){
      this->_simulation_data.push_back(DataPersistenceByGen(reserve_length_mg));
    }
  }
  
};

template <int number_replicates>
class DataContainer <DataAlleleFreq, number_replicates> :
  public DataContainer_Base_Persistence_By_Gen<DataAlleleFreq, number_replicates>, public DataContainer_Base_Allele_A_Freq<DataAlleleFreq> {
public:
  // constructor for infinite approx and to output raw allele data
  // DataContainer(int reserve_length_af){
  //   for (int i = 0; i < number_replicates; i++){
  //     this->_simulation_data.push_back(DataAlleleFreq(reserve_length_af));
  //   }
  // }
  // constructor for storing (in addition to infinite approx) both raw allele data and persistence_prob by gen
  DataContainer(int reserve_length_mg, int reserve_length_af){
    for (int i = 0; i < number_replicates; i++){
      this->_simulation_data.push_back(DataAlleleFreq(reserve_length_mg, reserve_length_af));
    }
  }

  double get_allele_A_freqs(int replicate, int gen){
    return this->_simulation_data[replicate]._allele_A_freq_by_gen[gen];
  }
  
};
  
#endif

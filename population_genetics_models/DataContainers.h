#ifndef DATA_CONTAINERS_H
#define DATA_CONTAINERS_H

#include <vector>
#include "Data.h"

template <class T>
class DataContainer_Base_Infinite_Approx {
protected:
  std::vector<T> _simulation_data;
public:
  virtual const double get_persistence_infinite_approx(){
    int persistence_count = 0;
    for (auto it = _simulation_data.begin(); it != _simulation_data.end(); it++){
      persistence_count += it->_persistence;
    }
    return static_cast<double>(persistence_count) / static_cast<double>(_simulation_data.size());
  };
  
  virtual void set_persistence_outcome_infinite(const int replicate, const bool outcome){
    this->_simulation_data[replicate]._persistence = outcome;
  }
  
};

template <class T, int number_replicates>
class DataContainer_Base_Persistence_By_Gen : public virtual DataContainer_Base_Infinite_Approx<T> {
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
  
  virtual const void append_persistence(const int replicate, const bool persistence){
    this->_simulation_data[replicate].append_persistence_by_gen(persistence);
  }
  
};

template <class T>
class DataContainer_Base_Allele_A_Freq : public virtual  DataContainer_Base_Infinite_Approx<T> {
public:
  const double get_allele_A_freq(const int replicate, const int gen){
    return this->_simulation_data[replicate]._allele_A_freq_by_gen[gen];
  }
  
  void append_allele_A_freq(const int replicate, const double allele_A_freq){
    this->_simulation_data[replicate].append_allele_A_freq_by_gen(allele_A_freq);
  }
  
};

template <class T, int number_replicates>
class DataContainer : public DataContainer_Base_Infinite_Approx<T> {};

template <int number_replicates>
class DataContainer <DataPersistenceInfinite, number_replicates> : public
DataContainer_Base_Infinite_Approx<DataPersistenceInfinite> {
public:

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

  DataContainer(int reserve_length_mg){
    for (int i = 0; i < number_replicates; i++){
      this->_simulation_data.push_back(DataPersistenceByGen(reserve_length_mg));
    }
  }
  
};

template <int number_replicates>
class DataContainer <DataAlleleFreqAndPBG, number_replicates> :
  public DataContainer_Base_Persistence_By_Gen<DataAlleleFreqAndPBG, number_replicates>, public DataContainer_Base_Allele_A_Freq<DataAlleleFreqAndPBG> {
public:
  
  DataContainer(int reserve_length_mg, int reserve_length_af){
    for (int i = 0; i < number_replicates; i++){
      this->_simulation_data.push_back(DataAlleleFreqAndPBG(reserve_length_mg, reserve_length_af));
    }
  }
  
};

template <int number_replicates>
class DataContainer <DataAlleleFreq, number_replicates> :
  public virtual DataContainer_Base_Infinite_Approx<DataAlleleFreq>, public DataContainer_Base_Allele_A_Freq<DataAlleleFreq> {
public:
  
  DataContainer(int reserve_length_af){
    for (int i = 0; i < number_replicates; i++){
      this->_simulation_data.push_back(DataAlleleFreq(reserve_length_af));
    }
  }
};

#endif

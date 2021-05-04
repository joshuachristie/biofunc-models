#ifndef RECORD_DATA_H
#define RECORD_DATA_H

#include "DataContainer.h"
#include "persistence_status.h"
#include "fixed_parameters.h"

namespace record {

  template<class P>
  void trait_freq(const std::vector<double> &trait_freq, const P &parameters, const int replicate,
		  DataContainer &data){
    data.append_trait_freq(replicate, trait_freq[parameters.shared.trait_info[0]]);
    // (re)reserve a large buffer if the replicate runs for a while (to avoid excessive memory allocations)
    std::size_t current_size = data._simulation_data[replicate]._trait_freq_by_gen.size();
    std::size_t current_capacity = data._simulation_data[replicate]._trait_freq_by_gen.capacity();
    if (current_size == current_capacity){
      data._simulation_data[replicate]._trait_freq_by_gen.reserve(current_capacity *
								  fixed_parameters::factor_to_expand_vector_memory);
    }
  }

  template<class P>
  void trait_presence_infinite(const std::vector<double> &trait_freq, const P &parameters, const int replicate,
			       DataContainer &data){
    if (persist_status::trait_extinct(trait_freq, parameters)){
      data.set_persistence_outcome_infinite(replicate, false);
    } else { // trait is either fixed or exists at a non-zero proportion
      data.set_persistence_outcome_infinite(replicate, true);
    }
  }

  template<class P>
  void trait_presence_by_gen(const std::vector<double> &trait_freq, const P &parameters, const int replicate,
			     DataContainer &data){
    if (persist_status::trait_extinct(trait_freq, parameters)){
      data.append_persistence(replicate, false);
    } else { // trait is either fixed or exists at a non-zero proportion
      data.append_persistence(replicate, true);
    }
  }
  
  void trait_freq(const double allele_A_freq, const int replicate, DataContainer &data);
  
}

#endif

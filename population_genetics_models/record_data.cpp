#include "record_data.h"
#include "DataContainer.h"
#include "fixed_parameters.h"

namespace record {
  
  void A_allele_freq(const double allele_A_freq, const int replicate, DataContainer &data){
    
    data.append_allele_A_freq(replicate, allele_A_freq);
    // (re)reserve a large buffer if the replicate runs for a while (to avoid excessive memory allocations)
    std::size_t current_size = data._simulation_data[replicate]._allele_A_freq_by_gen.size();
    std::size_t current_capacity = data._simulation_data[replicate]._allele_A_freq_by_gen.capacity();
    if (current_size == current_capacity){
      data._simulation_data[replicate]._allele_A_freq_by_gen.reserve(current_capacity *
								     fixed_parameters::factor_to_expand_vector_memory);
    }
  }
  
}

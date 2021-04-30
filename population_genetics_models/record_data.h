#ifndef RECORD_DATA_H
#define RECORD_DATA_H

#include "DataContainer.h"
#include "persistence_status.h"

namespace record {
  
  template<class P>
  void trait_presence_infinite(const double allele_A_freq, const P &parameters, const int replicate,
			       DataContainer &data){
    if (persist_status::is_not_extinct(allele_A_freq, parameters)){
      data.set_persistence_outcome_infinite(replicate, true);
    } else {
      data.set_persistence_outcome_infinite(replicate, false);
    }
  }

  template<class P>
  void trait_presence_by_gen(const double allele_A_freq, const P &parameters, const int replicate,
			     DataContainer &data){
    persist_status::is_not_extinct(allele_A_freq, parameters) ? data.append_persistence(replicate, true) :
      data.append_persistence(replicate, false);
  }
  
  void trait_freq(const double allele_A_freq, const int replicate, DataContainer &data);
  
}

#endif

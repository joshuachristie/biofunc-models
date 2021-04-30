#ifndef TRAIT_FREQ_H
#define TRAIT_FREQ_H

#include <vector>
#include "Parameters.h"

namespace trait_freq {
  
  template <class P>
  std::vector<double> initialise_trait_freq(const P &params){
    // params.model.trait_info[0] is index of trait of interest
    // params.model.trait_info[1] is the number of traits to track
    // haploid trait_freq is { P(A) }; diploid trait freq is { P(AA), P(Aa) }
    std::vector<double> trait_freq(params.shared.trait_info[1]);
    trait_freq[ params.shared.trait_info[0] ] = params.shared.initial_trait_freq;
    return trait_freq;
  }

}
#endif

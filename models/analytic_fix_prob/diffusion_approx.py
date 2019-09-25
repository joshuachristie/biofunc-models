import numpy as np

def fixation_probability(ploidy_level, effective_pop_size,
                         selection_coefficient, frequency):
    return (1 - np.exp(-2 * ploidy_level * effective_pop_size *
                          selection_coefficient * frequency, dtype=np.float128)) / (
                              1 - np.exp(-2 * ploidy_level * effective_pop_size *
                                            selection_coefficient,dtype=np.float128))

fix_prob = np.array([fixation_probability(1, N, s, 1/N) for s in np.linspace(-1,1,998)])

import numpy as np

def haploid_fix_prob(eff_pop_size, select_coef, initial_freq):
    num = 1 - np.exp(-2 * eff_pop_size * select_coef * initial_freq, dtype=np.float128)
    dem = 1 - np.exp(-2 * eff_pop_size * select_coef, dtype=np.float128)
    return num / dem


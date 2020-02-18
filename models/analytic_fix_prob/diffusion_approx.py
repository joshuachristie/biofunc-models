import numpy as np

def fix_prob(eff_pop_size, select_coef, ploidy):
    num = 1 - np.exp(-2 * (ploidy * eff_pop_size) * select_coef * (1 / (ploidy * eff_pop_size)), dtype=np.float128)
    dem = 1 - np.exp(-2 * (ploidy * eff_pop_size) * select_coef, dtype=np.float128)
    return num / dem

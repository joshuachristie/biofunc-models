""" 
Starting with a basic model with the following properties:
- haploid with two alleles (p and q)
- two environmental types (E_p and E_q)
- "infinite" number of environments with frequency P(E_p) and P(E_q)
- variable fitness of each allele in each environment (p matched to E_p; q matched to E_q)
- model can switch between deterministic and stochastic
""" 

import numpy as np
import matplotlib.pyplot as plt

def get_fitness_function(selection_coefficient):
    """
    fitness function relating haplotype to environment
    row 0 stores fitness in E_p; row 1 stores fitness in E_q
    col 0 stores fitness of allele p; col 1 stores fitness of allele q
    w_pp w_pq
    w_qp w_qq
    """
    haplotype_fitness = np.zeros((2,2))
    for i in range(2):
        haplotype_fitness[0, i] = 1 - selection_coefficient * i
        haplotype_fitness[1, i] = 1 - selection_coefficient * (1 - i)
    return haplotype_fitness


# PARAMETER VALUES
selection_coefficient = 0.1

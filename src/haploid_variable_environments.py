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

def get_environment_distribution(environment_distribution, environment_transition_matrix):
    """
    inputs: 
    environment_transition_matrix, which is of the form 
    e_00 e_01
    e_10 e_11
    where e_xy gives the probability of transitioning from env x to env y
    environment_distribution, which stores the frequency of E_p in index 0 and E_q in index 1

    output: the environment_distribution for the next generation
    """
    return np.matmul(environment_distribution, environment_transition_matrix)

def deterministic_allele_freq_after_selection(allele_freq, environment_distribution, haplotype_fitness):
    """
    inputs: 
    allele_freq (pre-selection)
    environment_distribution
    haplotype_fitness

    output: allele_freq (post-selection)
    """
    # multiply the allele freq by the fitnesses of the allele in each environment
    # weighted by the proportion of each environmental type
    # e.g. for p, the calc is p( P(E_p)w_pp + P(E_q)w_qp )
    unnormalised_allele_freq = np.multiply(
        allele_freq, np.matmul(environment_distribution, haplotype_fitness)
        )
    return np.divide(unnormalised_allele_freq, unnormalised_allele_freq.sum())

# PARAMETERS
selection_coefficient = 0.1
p_starting_freq = 0.5 # frequency of p allele at initialisation
env_p_starting_freq = 0.2 # frequency of E_p at initialisation
env_trans_p_to_q = 0.5 # transition rate of of E_p to E_q
env_trans_q_to_p = 0.5 # transition rate of of E_q to E_p
environment_transition_matrix = np.array([(1 - env_trans_p_to_q, env_trans_p_to_q),
                                          (env_trans_q_to_p, 1 - env_trans_q_to_p)])
haplotype_fitness = get_fitness_function(selection_coefficient)
# VARIABLES
allele_freq = np.array([p_starting_freq, 1 - p_starting_freq])
environment_distribution = np.array([env_p_starting_freq, 1 - env_p_starting_freq])

# store information from each run of a simulation
freq_p_allele_over_time = np.empty(0)
freq_environment_a_over_time = np.empty(0)


# script

    
    

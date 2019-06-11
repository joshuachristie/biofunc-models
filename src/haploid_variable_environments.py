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
from sys import argv
import math

def get_fitness_function(selection_coefficient_p1, selection_coefficient_p2,
                         selection_coefficient_q1, selection_coefficient_q2):
    """
    fitness function relating haplotype to environment
    row 0 stores fitness in E_1; row 1 stores fitness in E_2
    col 0 stores fitness of allele p; col 1 stores fitness of allele q
    w_p1 w_q1
    w_p2 w_q2
    where w_xy is read as "the fitness of x in env y"
    """
    haplotype_fitness = np.zeros((2,2))
    
    haplotype_fitness[0, 0] = 1 + selection_coefficient_p1
    haplotype_fitness[1, 0] = 1 + selection_coefficient_p2
    haplotype_fitness[0, 1] = 1 + selection_coefficient_q1
    haplotype_fitness[1, 1] = 1 + selection_coefficient_q2
        
    return haplotype_fitness

def get_environment_distribution(environment_distribution, environment_transition_matrix):
    """
    inputs: 
    environment_transition_matrix, which is of the form 
    r_11 r_12
    r_21 r_22
    where r_xy gives the transition rate from env x to env y
    environment_distribution, which stores the frequency of E_1 in index 0 (e_1) and E_2 in index 1 (e_2)

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
    # e.g. for Ap, the calc is p( e_1w_p1 + e_2w_p2 )
    unnormalised_allele_freq = np.multiply(
        allele_freq, np.matmul(environment_distribution, haplotype_fitness)
        )
    return np.divide(unnormalised_allele_freq, unnormalised_allele_freq.sum())

def has_allele_fixed(threshold, frequency_tm1, frequency_t):
    return abs(frequency_tm1 - frequency_t) < threshold

def num_gens_to_pass_threshold(threshold, gen, threshold_gen):
    if allele_freq[0] >= threshold and math.isnan(threshold_gen):
        threshold_gen = int(gen)
    return threshold_gen


# PARAMETERS
script, scpp, scqq, scpq, scqp, psf, epsf, et12, et21, fn = argv

selection_coefficient_p1 = float(scpp) # benefit to p of being in env 1
selection_coefficient_p2 = float(scpq) # benefit to p of being in env 2
selection_coefficient_q1 = float(scqp) # benefit to q of being in env 1
selection_coefficient_q2 = float(scqq) # benefit to q of being in env 2
p_starting_freq = float(psf) # frequency of p allele at initialisation (p_0)
env_p_starting_freq = float(epsf) # frequency of E_p at initialisation (e_p0)
env_trans_12 = float(et12) # transition rate of E_1 to E_2
env_trans_21 = float(et21) # transition rate of E_2 to E_1

environment_transition_matrix = np.array([(1 - env_trans_12, env_trans_12),
                                          (env_trans_21, 1 - env_trans_21)])
haplotype_fitness = get_fitness_function(selection_coefficient_p1, selection_coefficient_p2,
                         selection_coefficient_q1, selection_coefficient_q2)
fixation_threshold = 0.000000001
# VARIABLES
allele_freq = np.array([p_starting_freq, 1 - p_starting_freq])
environment_distribution = np.array([env_p_starting_freq, 1 - env_p_starting_freq]) # (e_1, e_2)

# store information from each run of a simulation
freq_p_allele_over_time = np.empty(0)
freq_env_p_over_time = np.empty(0)
allele_freq_tm1 = -1 # dummy value to ensure it enters while loop
gen = 0
threshold_gen = float('nan')
# script
while not has_allele_fixed(fixation_threshold, allele_freq_tm1, allele_freq[0]):
    # record information about sim
    freq_p_allele_over_time = np.append(freq_p_allele_over_time, allele_freq[0])
    freq_env_p_over_time = np.append(freq_env_p_over_time, environment_distribution[0])
    allele_freq_tm1 = allele_freq[0]
    
    # allele frequencies after selection
    allele_freq = deterministic_allele_freq_after_selection(
        allele_freq, environment_distribution, haplotype_fitness)
    # environment distribution after change
    environment_distribution = get_environment_distribution(
        environment_distribution, environment_transition_matrix)
    gen += 1
    threshold_gen = num_gens_to_pass_threshold(1 - p_starting_freq, gen, threshold_gen)

fig, ax1 = plt.subplots(1, 1, sharex=True, sharey=True)
fig.text(0.5, 0.02, 'generations', ha='center')
fig.text(0.02, 0.5, 'frequency', va='center', rotation='vertical')
ax1.plot(freq_p_allele_over_time, '--k', label='$p$')
ax1.plot(freq_env_p_over_time, 'r', label = '$e_p$')
legend1 = ax1.legend(loc='center right', shadow=True, fontsize='x-large')
legend1.get_frame().set_facecolor('C0')

filename = '/home/joshua/insync-gdrive/projects/active/seleff/models/docs/figs/{}.png'.format(fn)
_ = plt.savefig(filename)
# _ = plt.show()

if not math.isnan(threshold_gen):
    print('Time to fixation = {} generations'.format(threshold_gen))
else:
    print('p did not fix; p = {0:1.2f}'.format(allele_freq[0]))

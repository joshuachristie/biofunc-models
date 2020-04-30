"""
The wright-fisher haploid single environment model doesn't require partial information decomposition.
As such, the function is wholly apportioned to the haploid trait.
"""
from subprocess import Popen, PIPE

selection_coefficient = '0.0'
population_size = '1000'
number_generations = '100000'
number_replicates = '100000'
binary = '/home/joshua/projects/metric/models/wright_fisher/haploid_single_environment/hse'

process = Popen([binary, selection_coefficient, population_size, number_generations, number_replicates], stdout=PIPE)

probability_persistence = float(process.stdout.read().strip())

filename = '../../data/HSE/s_{}_N_{}_g_{}_r_{}.csv'.format(
    selection_coefficient, population_size, number_generations, number_replicates)
with open(filename,'w') as f:
    f.write(str(probability_persistence))
    f.write('\n')

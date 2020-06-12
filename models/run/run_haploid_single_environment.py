"""
The wright-fisher haploid single environment model doesn't require partial information decomposition.
As such, the function metric is wholly apportioned to the haploid trait.
I need to output two persistence probabilities:
    0th element: neutral drift (i.e. selection_coefficient absent)
    1st element: selection_coefficient present
"""
from subprocess import Popen, PIPE

population_size = '1000'
number_generations = '100000'
number_replicates = '100000'
binary = '../wright_fisher/haploid_single_environment/hse'

selection_coefficient = '0.01'
SC = [0.0, selection_coefficient]

persistence_probabilities = [0] * 2
counter = 0
for i in range(2):
    process = Popen([binary, str(SC[i]), population_size,
                     number_generations, number_replicates], stdout=PIPE)

    persistence_probabilities[counter] = float(process.stdout.read().strip())
    counter +=1
    
filename = '../../data/persistence_probabilities/HSE/s_{}_N_{}_g_{}_r_{}.csv'.format(
    selection_coefficient, population_size, number_generations, number_replicates)

with open(filename,'w') as f:
    [f.write(str(persistence_probabilities[i]) + ',') if i < len(persistence_probabilities) - 1
    else f.write(str(persistence_probabilities[i])) for i in range(len(persistence_probabilities))]
    f.write('\n')

"""
The wright-fisher diploid single environment model requires partial information decomposition.
I need to run four combinations of parameter values to get the full set of inputs.
In total, this will give eight outputs (the four combinations for probability_persistence and 
1 - probability_persistence).

The probabilities that the allele is not lost will be stored in persistence_probabilities.
The order is:
    0th element: neutral drift
    1st element: heterozygote effect alone
    2nd element: homozygote effect alone
    3rd element: focal situation (effect of both heterozygote and homozygote)
"""
from subprocess import Popen, PIPE

population_size = '1000'
number_generations = '10000'
number_replicates = '10000'
binary = '/home/joshua/projects/metric/models/wright_fisher/diploid_single_environment/dse'

selection_coefficient_homozygote = 0.2
selection_coefficient_heterozygote = 0.1
SC_homo = [0.0, selection_coefficient_homozygote]
SC_hetero = [0.0, selection_coefficient_heterozygote]

persistence_probabilities = [0] * 4
counter = 0
for i in range(2):
    for j in range(2):
        
        process = Popen([binary, str(SC_homo[i]), str(SC_hetero[j]),
                         population_size, number_generations, number_replicates], stdout=PIPE)

        persistence_probabilities[counter] = float(process.stdout.read().strip())
        counter += 1
        
print(persistence_probabilities)

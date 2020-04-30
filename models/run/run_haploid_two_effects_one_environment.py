"""
The wright-fisher haploid two effects one environment model requires partial information decomposition.
The two alleles are A and a whose fitnesses are given by wA = sA1 + sA2 and wa = sa1 + sa2
where 1 and 2 indicate the environment.
I will fix the value of wa (and thus both sa1 and sa2), which means I need to vary SA1 and SA2.
I thus will run four combinations of parameter values to get the full set of inputs (fully factorial
SA1/SA2 with their values and 0.0).
In total, this will give eight outputs (the four combinations for probability_persistence and 
1 - probability_persistence).

The probabilities that the allele is not lost will be stored in persistence_probabilities.
The order is:
    0th element: neutral drift
    1st element: selection_coefficient_A2 alone
    2nd element: selection_coefficient_A1 alone
    3rd element: focal situation (effect of both selection_coefficient_A1 and selection_coefficient_A2)
"""
from subprocess import Popen, PIPE

population_size = '1000'
number_generations = '100000'
number_replicates = '100000'
binary = '/home/joshua/projects/metric/models/wright_fisher/haploid_two_effects_one_environment/hteoe'

selection_coefficient_A1 = '0.015'
selection_coefficient_A2 = '0.03'
selection_coefficient_a1 = '0.01'
selection_coefficient_a2 = '0.02'
SC_A1 = [0.0, selection_coefficient_A1]
SC_A2 = [0.0, selection_coefficient_A2]

persistence_probabilities = [0] * 4
counter = 0
for i in range(2):
    for j in range(2):
        process = Popen([binary, str(SC_A1[i]), str(SC_A2[j]), selection_coefficient_a1,
                         selection_coefficient_a2, population_size,
                         number_generations, number_replicates], stdout=PIPE)

        persistence_probabilities[counter] = float(process.stdout.read().strip())
        counter += 1

filename = '../../data/HTEOE/A1_{}_A2_{}_a1_{}_a2_{}_N_{}_g_{}_r_{}.csv'.format(
    selection_coefficient_A1, selection_coefficient_A2, selection_coefficient_a1,
    selection_coefficient_a2, population_size, number_generations, number_replicates)

with open(filename,'w') as f:
    [f.write(str(persistence_probabilities[i]) + ',') if i < len(persistence_probabilities) - 1
    else f.write(str(persistence_probabilities[i])) for i in range(len(persistence_probabilities))]
    f.write('\n')

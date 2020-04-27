
from subprocess import Popen, PIPE

selection_coefficient_A1 = '0.0'
selection_coefficient_A2 = '0.0'
selection_coefficient_a1 = '0.0'
selection_coefficient_a2 = '0.0'
population_size = '1000'
number_generations = '100000'
number_replicates = '100000'
binary = '/home/joshua/projects/metric/models/wright_fisher/haploid_two_effects_one_environment/hteoe'

process = Popen([binary, selection_coefficient_A1, selection_coefficient_A2, selection_coefficient_a1,
                 selection_coefficient_a2, population_size, number_generations, number_replicates], stdout=PIPE)

probability_persistence = float(process.stdout.read().strip())
print(probability_persistence)

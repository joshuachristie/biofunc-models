from subprocess import Popen, PIPE

selection_coefficient = '0.0'
population_size = '1000'
number_generations = '100000'
number_replicates = '100000'
binary = '/home/joshua/projects/metric/models/wright_fisher/haploid_single_environment/hse'

process = Popen([binary, selection_coefficient, population_size, number_generations, number_replicates], stdout=PIPE)

probability_persistence = float(process.stdout.read().strip())
print(probability_persistence)

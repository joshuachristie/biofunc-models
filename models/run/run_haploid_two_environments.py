from subprocess import Popen, PIPE

selection_coefficient_A_env_1 = '0.0'
selection_coefficient_A_env_2 = '0.0'
selection_coefficient_a_env_1 = '0.0'
selection_coefficient_a_env_2 = '0.0'
population_size = '1000'
gen_env_1 = '5000'
gen_env_2 = '5000'
number_replicates = '100000'
binary = '/home/joshua/projects/metric/models/wright_fisher/haploid_two_environments/hte'

process = Popen([binary, selection_coefficient_A_env_1, selection_coefficient_A_env_2,
                 selection_coefficient_a_env_1, selection_coefficient_a_env_2, population_size,
                 gen_env_1, gen_env_2, number_replicates], stdout=PIPE)
                 
probability_persistence = float(process.stdout.read().strip())
print(probability_persistence)

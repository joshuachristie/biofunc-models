"""
The wright-fisher haploid two environment model requires partial information decomposition.
There are two alleles (traits), A and a, and two environments.
First environment 1 (for gen_env_1 generations) and then environment 2 (for gen_env_2 generations).
Each allele has a fitness for each environment.
For allele A, sA_env_1 and sA_env_2 and for allele a, sa_env_1 and sa_env_2.
I will fix the value of wa (and thus both sa_env_1 and sa_env_2), which means I need to vary SA_env_1 and SA_env_2.
I thus will run four combinations of parameter values to get the full set of inputs (fully factorial
SA_env_1/SA_env_2 with their values and 0.0).
In total, this will give eight outputs (the four combinations for probability_persistence and 
1 - probability_persistence).

The probabilities that the allele is not lost will be stored in persistence_probabilities.
The order is:
    0th element: neutral drift
    1st element: selection_coefficient_A2 alone
    2nd element: selection_coefficient_A1 alone
    3rd element: focal situation (effect of both selection_coefficient_A_env_1 and selection_coefficient_A_env_2)
"""
from subprocess import Popen, PIPE
import argparse

# parse parameters
parser = argparse.ArgumentParser()
parser.add_argument("population_size", type=int,
                    help="number of individuals in the population")
parser.add_argument("gen_env_1", type=int,
                    help="number of generations that environment 1 persists for")
parser.add_argument("gen_env_2", type=int,
                    help="number of generations that environment 2 persists for")
parser.add_argument("number_replicates", type=int,
                    help="number of sims to run for each parameter value set")
parser.add_argument("selection_coefficient_A_env_1", type=float,
                    help="selection coefficient of allele A in environment 1")
parser.add_argument("selection_coefficient_A_env_2", type=float,
                    help="selection coefficient of allele A in environment 2")
parser.add_argument("selection_coefficient_a_env_1", type=float,
                    help="selection coefficient of allele a in environment 1")
parser.add_argument("selection_coefficient_a_env_2", type=float,
                    help="selection coefficient of allele a in environment 2")
args = parser.parse_args()

SC_A_e_1 = [0.0, args.selection_coefficient_A_env_1]
SC_A_e_2 = [0.0, args.selection_coefficient_A_env_2]
binary = '../../../population_genetics_models/haploid_two_environments/hte'

persistence_probabilities = [0] * 4
counter = 0
for i in range(2):
    for j in range(2):
        process = Popen([binary, str(SC_A_e_1[i]), str(SC_A_e_2[j]), str(args.selection_coefficient_a_env_1),
                 str(args.selection_coefficient_a_env_2), str(args.population_size), str(args.gen_env_1),
                         str(args.gen_env_2), str(args.number_replicates)], stdout=PIPE)
        persistence_probabilities[counter] = float(process.stdout.read().strip())
        counter += 1

filename = '../../../data/persistence_probabilities/HTE/A1_{}_A2_{}_a1_{}_a2_{}_N_{}_g1_{}_g2_{}_r_{}.csv'.format(
    args.selection_coefficient_A_env_1, args.selection_coefficient_A_env_2, args.selection_coefficient_a_env_1,
    args.selection_coefficient_a_env_2, args.population_size, args.gen_env_1, args.gen_env_2,
    args.number_replicates)

with open(filename,'w') as f:
    [f.write(str(persistence_probabilities[i]) + ',') if i < len(persistence_probabilities) - 1
    else f.write(str(persistence_probabilities[i])) for i in range(len(persistence_probabilities))]
    f.write('\n')

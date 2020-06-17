"""
The wright-fisher haploid two effects one environment model requires partial information decomposition.
The two alleles are A and a whose fitnesses are given by wA = sA1 + sA2 and wa = sa1 + sa2
where 1 and 2 indicate the two effects of the trait in the environment.
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
import argparse

# parse parameters
parser = argparse.ArgumentParser()
parser.add_argument("population_size", type=int,
                    help="number of individuals in the population")
parser.add_argument("number_generations", type=int,
                    help="number of generations to run the simulation for")
parser.add_argument("number_replicates", type=int,
                    help="number of sims to run for each parameter value set")
parser.add_argument("selection_coefficient_A1", type=float,
                    help="selection coefficient of effect 1 of allele A")
parser.add_argument("selection_coefficient_A2", type=float,
                    help="selection coefficient of effect 2 of allele A")
parser.add_argument("selection_coefficient_a1", type=float,
                    help="selection coefficient of effect 1 of allele a")
parser.add_argument("selection_coefficient_a2", type=float,
                    help="selection coefficient of effect 2 of allele a")
args = parser.parse_args()

SC_A1 = [0.0, args.selection_coefficient_A1]
SC_A2 = [0.0, args.selection_coefficient_A2]
binary = '../../../population_genetics_models/haploid_two_effects_one_environment/hteoe'

persistence_probabilities = [0] * 4
counter = 0
for i in range(2):
    for j in range(2):
        process = Popen([binary, str(SC_A1[i]), str(SC_A2[j]), str(args.selection_coefficient_a1),
                         str(args.selection_coefficient_a2), str(args.population_size),
                         str(args.number_generations), str(args.number_replicates)], stdout=PIPE)
        persistence_probabilities[counter] = float(process.stdout.read().strip())
        counter += 1

filename = '../../../data/persistence_probabilities/HTEOE/A1_{}_A2_{}_a1_{}_a2_{}_N_{}_g_{}_r_{}.csv'.format(
    args.selection_coefficient_A1, args.selection_coefficient_A2, args.selection_coefficient_a1,
    args.selection_coefficient_a2, args.population_size, args.number_generations, args.number_replicates)

with open(filename,'w') as f:
    [f.write(str(persistence_probabilities[i]) + ',') if i < len(persistence_probabilities) - 1
    else f.write(str(persistence_probabilities[i])) for i in range(len(persistence_probabilities))]
    f.write('\n')

"""
The wright-fisher diploid single environment model requires partial information decomposition.
I need to run four combinations of parameter values to get the full set of inputs.
In total, this will give eight outputs (the four combinations for probability_persistence and 
1 - probability_persistence).

The probabilities that the allele is not lost will be stored in persistence_probabilities.
The order is:
    0th element: neutral drift
    1st element: heterozygote effect alone
    2nd element: homozygote (AA) effect alone
    3rd element: focal situation (effect of both heterozygote and homozygote)
"""
from subprocess import Popen, PIPE
import argparse

# parse parameters
parser = argparse.ArgumentParser()
parser.add_argument("population_size", type=int,
                    help="number of individuals in the population")
parser.add_argument("number_generations", type=int,
                    help="number of generations to run simulation for")
parser.add_argument("number_replicates", type=int,
                    help="number of sims to run for each parameter value set")
parser.add_argument("selection_coefficient_homozygote", type=float,
                    help="selection coefficient of AA homozygote")
parser.add_argument("selection_coefficient_heterozygote", type=float,
                    help="selection coefficient of allele A in environment 2")
args = parser.parse_args()

binary = '../../../population_genetics_models/diploid_single_environment/dse'
SC_homo = [0.0, args.selection_coefficient_homozygote]
SC_hetero = [0.0, args.selection_coefficient_heterozygote]

persistence_probabilities = [0] * 4
counter = 0
for i in range(2):
    for j in range(2):
        
        process = Popen([binary, str(SC_homo[i]), str(SC_hetero[j]),
                         str(args.population_size), str(args.number_generations), str(args.number_replicates)],
                        stdout=PIPE)
        persistence_probabilities[counter] = float(process.stdout.read().strip())
        counter += 1
filename = '../../../data/persistence_probabilities/DSE/ho_{}_he_{}_N_{}_g_{}_r_{}.csv'.format(
    args.selection_coefficient_homozygote, args.selection_coefficient_heterozygote, args.population_size,
    args.number_generations, args.number_replicates)
with open(filename,'w') as f:
    [f.write(str(persistence_probabilities[i]) + ',') if i < len(persistence_probabilities) - 1
    else f.write(str(persistence_probabilities[i])) for i in range(len(persistence_probabilities))]
    f.write('\n')


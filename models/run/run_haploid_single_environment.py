"""
The wright-fisher haploid single environment model doesn't require partial information decomposition.
As such, the function metric is wholly apportioned to the haploid trait.
I need to output two persistence probabilities:
    0th element: neutral drift (i.e. selection_coefficient absent)
    1st element: selection_coefficient present
"""
from subprocess import Popen, PIPE
import argparse

# parse parameters
parser = argparse.ArgumentParser()
parser.add_argument("population_size", type=int,
                    help="number of individuals in the population")
parser.add_argument("number_generations", type=int,
                    help="number of generations to run each simulation for")
parser.add_argument("number_replicates", type=int,
                    help="number of sims to run for each parameter value set")
parser.add_argument("selection_coefficient", type=float,
                    help="selection coefficient of allele A")
args = parser.parse_args()

SC = [0.0, args.selection_coefficient]
binary = '../wright_fisher/haploid_single_environment/hse'

persistence_probabilities = [0] * 2
counter = 0
for i in range(2):
    process = Popen([binary, str(SC[i]), str(args.population_size),
                     str(args.number_generations), str(args.number_replicates)], stdout=PIPE)
    persistence_probabilities[counter] = float(process.stdout.read().strip())
    counter +=1
    
filename = '../../data/persistence_probabilities/HSE/s_{}_N_{}_g_{}_r_{}.csv'.format(
    args.selection_coefficient, args.population_size, args.number_generations, args.number_replicates)

with open(filename,'w') as f:
    [f.write(str(persistence_probabilities[i]) + ',') if i < len(persistence_probabilities) - 1
    else f.write(str(persistence_probabilities[i])) for i in range(len(persistence_probabilities))]
    f.write('\n')

# model for kimura's diffusion approximation
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../models/analytic_fix_prob"))
from diffusion_approx import haploid_fix_prob
# modules for running the C++ model
import subprocess
# some general packages
import numpy as np
import matplotlib.pyplot as plt

# parameter values for Kimura's diffusion approximation
N = 100 # N is effective population size; s is selection coefficient
fix_prob_vs_s = np.array([haploid_fix_prob(N, s, 1/N) for s in np.linspace(-1,1,1000)]) 
# set up parallel subprocesses to run C++ wright-fisher model
# parameter values for wright-fisher (in addition to those specified in the diffusion approx)
max_gen = 10000
num_rep = 1000000

cmd_line = ("./wright_fisher", "{}", "{}", "{}", "{}".format(N, s, max_gen, num_rep))


# plot diffusion approx vs simulation fixation prob
plt.plot(np.linspace(-1,1,1000), fix_prob_vs_s)
plt.xlabel('selection coefficient')
plt.ylabel('fixation probability')

plt.savefig('../figs/test_{}.jpg'.format(N))

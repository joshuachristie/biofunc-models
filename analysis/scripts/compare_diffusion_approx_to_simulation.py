# import model for kimura's diffusion approximation
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../models/analytic_fix_prob"))
from diffusion_approx import haploid_fix_prob
#
# import some general packages
import numpy as np
import matplotlib.pyplot as plt

# first, plot Kimura's diffusion approximation (for a given N and a range of s)
# N is effective population size; s is selection coefficient
N = 100

fix_prob_vs_s = np.array([haploid_fix_prob(N, s, 1/N) for s in np.linspace(-1,1,1000)])

plt.plot(np.linspace(-1,1,1000), fix_prob_vs_s)
plt.xlabel('selection coefficient')
plt.ylabel('fixation probability')

plt.savefig('../figs/test_{}.jpg'.format(N))


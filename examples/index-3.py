from astro_traj.galaxy import Hernquist_NFW
from astro_traj import constr_dict
from astro_traj.sample import Sample
import numpy as np
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
font = {'size': 22}

samples = 'posterior_samples.dat'
Galaxy = constr_dict.galaxy('NGC', samples, 2.8, 5, 0.679)
gal = Hernquist_NFW(Galaxy['Mspiral'], Galaxy['Mbulge'], Galaxy['Mhalo'], Galaxy['R_eff'], 0.679, rcut=100)
samp = Sample(gal)
Nsys = 1000
bins = int(np.round(np.sqrt(Nsys)))

Apre_dist_log = samp.sample_Apre(Amin=0.1, Amax=10.0, method='log', size=Nsys)
Apre_dist_uniform = samp.sample_Apre(Amin=0.1, Amax=10.0, method='uniform', size=Nsys)

plot, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(18.5, 10.5))
ax1.hist(Apre_dist_log, bins=bins)
ax2.hist(Apre_dist_uniform, bins=bins)
ax1.set_xlabel('Apre: Log', fontdict=font)
ax2.set_xlabel('Apre: Uniform', fontdict=font)
plot.show()

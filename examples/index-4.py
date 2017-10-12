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

epre_dist_circularized = samp.sample_epre(method='circularized', size=Nsys)

plot, ax1 = plt.subplots(1, sharex=True)
ax1.hist(epre_dist_circularized, bins=bins)
ax1.set_xlabel('Eccentricity Pre-Supernova: Circularized', fontdict=font)
plot.show()

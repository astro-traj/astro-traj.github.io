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

d_dist_posterior = samp.sample_distance(samples, method='posterior', size=Nsys)
d_dist_median = samp.sample_distance(samples, method='median', size=Nsys)
d_dist_mean = samp.sample_distance(samples, method='mean', size=Nsys)
d_dist_gaussian = samp.sample_distance(samples, method='gaussian', size=Nsys)

plot, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, figsize=(18.5, 10.5))
ax1.hist(d_dist_posterior, bins=bins)
ax2.hist(d_dist_median, bins=bins)
ax3.hist(d_dist_mean, bins=bins)
ax4.hist(d_dist_gaussian, bins=bins)
ax1.set_xlabel('Distance: Posterior', fontdict=font)
ax2.set_xlabel('Distance: Median', fontdict=font)
ax3.set_xlabel('Distance: Mean', fontdict=font)
ax4.set_xlabel('Distance: Gaussian', fontdict=font)
plot.show()

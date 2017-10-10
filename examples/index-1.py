from astro_traj.galaxy import Hernquist_NFW
from astro_traj import constr_dict
from astro_traj.sample import Sample
import numpy as np
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt

samples = 'posterior_samples.dat'
Galaxy = constr_dict.galaxy('NGC', samples, 100, 5, 0.73)
gal = Hernquist_NFW(Galaxy['Mspiral'], Galaxy['Mbulge'], Galaxy['Mhalo'], Galaxy['R_eff'], 0.73, rcut=100)
samp = Sample(gal)
Nsys = 1000
bins = int(np.round(np.sqrt(Nsys)))

Mcomp_dist_posterior, Mns_dist_posterior = samp.sample_masses(samples, method='posterior', size=Nsys)
Mcomp_dist_median, Mns_dist_median = samp.sample_masses(samples, method='median', size=Nsys)
Mcomp_dist_mean, Mns_dist_mean = samp.sample_masses(samples, method='mean', size=Nsys)
Mcomp_dist_gaussian, Mns_dist_gaussian = samp.sample_masses(samples, method='gaussian', size=Nsys)


plot, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, sharex=True, figsize=(18.5, 10.5))
ax1.hist(Mcomp_dist_posterior, bins=bins)
ax2.hist(Mcomp_dist_median, bins=bins)
ax3.hist(Mcomp_dist_mean, bins=bins)
ax4.hist(Mcomp_dist_gaussian, bins=bins)
ax5.hist(Mns_dist_posterior, bins=bins)
ax6.hist(Mns_dist_median, bins=bins)
ax7.hist(Mns_dist_mean, bins=bins)
ax8.hist(Mns_dist_gaussian, bins=bins)
plot.show()

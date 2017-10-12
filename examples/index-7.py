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
dumrand = np.random.uniform(0,1,size=Nsys)

ECS,CCS = samp.initialize_Vkick()
Vkick_dist_uniform = samp.sample_Vkick(method='uniform', size=Nsys, ECSPDF=ECS, CCSPDF=CCS, Mhe=None, irand=dumrand)
Vkick_dist_maxwellian = samp.sample_Vkick(method='maxwellian', size=Nsys, ECSPDF=ECS, CCSPDF=CCS, Mhe=None, irand=dumrand)
Vkick_dist_beniamini2 = samp.sample_Vkick(method='beniamini2', size=Nsys, ECSPDF=ECS, CCSPDF=CCS, Mhe=None, irand=dumrand)

plot, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(18.5, 10.5))
ax1.hist(Vkick_dist_uniform, bins=bins)
ax2.hist(Vkick_dist_maxwellian, bins=bins)
ax3.hist(Vkick_dist_beniamini2, bins=bins)
ax1.set_xlabel('Supernova Kick Velocity: Uniform', fontdict=font)
ax2.set_xlabel('Supernova Kick Velocity: Maxwellian', fontdict=font)
ax3.set_xlabel('Supernova Kick Velocity: Beniamini2', fontdict=font)
plot.show()

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
Mcomp_dist, Mns_dist = samp.sample_masses(samples, method='posterior', size=Nsys)

ECSPDFMhe = samp.initialize_Mhe(0.1)
CCSPDFMhe = samp.initialize_Mhe(1.0)
dumrand = np.random.uniform(0,1,size=Nsys)
Mhe_dist_uniform = samp.sample_Mhe(Mmin=Mns_dist, method='uniform', size=Nsys, PDF=None, ECSPDF=ECSPDFMhe, CCSPDF=CCSPDFMhe, irand=dumrand) # (Msun)
Mhe_dist_power = samp.sample_Mhe(Mmin=Mns_dist, method='power', size=Nsys, PDF=None, ECSPDF=ECSPDFMhe, CCSPDF=CCSPDFMhe, irand=dumrand) # (Msun)
Mhe_dist_beniamini2 = samp.sample_Mhe(Mmin=Mns_dist, method='beniamini2', size=Nsys, PDF=None, ECSPDF=ECSPDFMhe, CCSPDF=CCSPDFMhe, irand=dumrand) # (Msun)

plot, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(18.5, 10.5))
ax1.hist(Mhe_dist_uniform, bins=bins)
ax2.hist(Mhe_dist_power, bins=bins)
ax3.hist(Mhe_dist_beniamini2, bins=bins)
ax1.set_xlabel('Mass Helium Star: Uniform', fontdict=font)
ax2.set_xlabel('Mass Helium Star: Power', fontdict=font)
ax3.set_xlabel('Mass Helium Star: Beniamini2', fontdict=font)
plot.show()

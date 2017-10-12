from astro_traj.galaxy import Hernquist_NFW
from astro_traj import constr_dict
from astro_traj.sample import Sample
from astro_traj.system import System
import numpy as np
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
font = {'size': 16}

samples = 'posterior_samples.dat'
Galaxy = constr_dict.galaxy('NGC', samples, 2.8, 5, 0.679)
gal = Hernquist_NFW(Galaxy['Mspiral'], Galaxy['Mbulge'], Galaxy['Mhalo'], Galaxy['R_eff'], 0.679, rcut=100)
samp = Sample(gal)
Nsys = 1000
bins = int(np.round(np.sqrt(Nsys)))
dumrand = np.random.uniform(0,1,size=Nsys)


Mcomp_dist, Mns_dist = samp.sample_masses(samples, method='posterior', size=Nsys)
d_dist = samp.sample_distance(samples, method='median', size=Nsys)
Apre_dist_log = samp.sample_Apre(Amin=0.1, Amax=10.0, method='log', size=Nsys)
Apre_dist_uniform = samp.sample_Apre(Amin=0.1, Amax=10.0, method='uniform', size=Nsys)
epre_dist = samp.sample_epre(method='circularized', size=Nsys)
PDFR = samp.initialize_R()
R_dist = samp.sample_R(PDFR, Nsys)
ECSPDFMhe = samp.initialize_Mhe(0.1)
CCSPDFMhe = samp.initialize_Mhe(1.0)
Mhe_dist_uniform = samp.sample_Mhe(Mmin=Mns_dist, method='uniform', size=Nsys, PDF=None, ECSPDF=ECSPDFMhe, CCSPDF=CCSPDFMhe, irand=dumrand) # (Msun)
Mhe_dist_power = samp.sample_Mhe(Mmin=Mns_dist, method='power', size=Nsys, PDF=None, ECSPDF=ECSPDFMhe, CCSPDF=CCSPDFMhe, irand=dumrand) # (Msun)
Mhe_dist_beniamini2 = samp.sample_Mhe(Mmin=Mns_dist, method='beniamini2', size=Nsys, PDF=None, ECSPDF=ECSPDFMhe, CCSPDF=CCSPDFMhe, irand=dumrand) # (Msun)
ECS, CCS = samp.initialize_Vkick()
Vkick_dist_uniform = samp.sample_Vkick(method='uniform', size=Nsys, ECSPDF=ECS, CCSPDF=CCS, Mhe=None, irand=dumrand)
Vkick_dist_maxwellian = samp.sample_Vkick(method='maxwellian', size=Nsys, ECSPDF=ECS, CCSPDF=CCS, Mhe=None, irand=dumrand)
Vkick_dist_beniamini2 = samp.sample_Vkick(method='beniamini2', size=Nsys, ECSPDF=ECS, CCSPDF=CCS, Mhe=None, irand=dumrand)

successful_binaries_log_ben2_ben2 = []
successful_binaries_uni_uni_uni = []
successful_binaries_uni_power_maxwell = []

for R, d, Mcomp, Mns, Apre, epre, Mhe, Vkick in zip(R_dist, d_dist, Mcomp_dist, Mns_dist, Apre_dist_log, epre_dist, Mhe_dist_beniamini2, Vkick_dist_beniamini2):
    T = System(gal, R, Mns, Mcomp, Mhe, Apre, epre, d, Vkick, sys_flag=None)
    T.SN()
    if T.flag != 3:
        successful_binaries_log_ben2_ben2.append(True)
    else:
        successful_binaries_log_ben2_ben2.append(False)

for R, d, Mcomp, Mns, Apre, epre, Mhe, Vkick in zip(R_dist, d_dist, Mcomp_dist, Mns_dist, Apre_dist_uniform, epre_dist, Mhe_dist_uniform, Vkick_dist_uniform):
    T = System(gal, R, Mns, Mcomp, Mhe, Apre, epre, d, Vkick, sys_flag=None)
    T.SN()
    if T.flag != 3:
        successful_binaries_uni_uni_uni.append(True)
    else:
        successful_binaries_uni_uni_uni.append(False)


for R, d, Mcomp, Mns, Apre, epre, Mhe, Vkick in zip(R_dist, d_dist, Mcomp_dist, Mns_dist, Apre_dist_uniform, epre_dist, Mhe_dist_power, Vkick_dist_maxwellian):
    T = System(gal, R, Mns, Mcomp, Mhe, Apre, epre, d, Vkick, sys_flag=None)
    T.SN()
    if T.flag != 3:
        successful_binaries_uni_power_maxwell.append(True)
    else:
        successful_binaries_uni_power_maxwell.append(False)

plot, ((ax1, ax4, ax7), (ax2, ax5, ax8), (ax3, ax6, ax9)) = plt.subplots(3, 3, sharex='col', figsize=(18.5, 10.5))
ax1.hist(Mhe_dist_beniamini2, bins=bins)
ax2.hist(Mhe_dist_uniform, bins=bins)
ax3.hist(Mhe_dist_power, bins=bins)
ax1.hist(Mhe_dist_beniamini2[successful_binaries_log_ben2_ben2], bins=bins)
ax2.hist(Mhe_dist_uniform[successful_binaries_uni_uni_uni], bins=bins)
ax3.hist(Mhe_dist_power[successful_binaries_uni_power_maxwell], bins=bins)
ax1.set_xlabel('Mhe: Apre: Log Mhe: ben2 Vkick: ben2', fontdict=font)
ax2.set_xlabel('Mhe: Apre: uni Mhe: uni Vkick: uni', fontdict=font)
ax3.set_xlabel('Mhe: Apre: uni Mhe: power Vkick: maxwell', fontdict=font)

ax4.hist(Vkick_dist_beniamini2, bins=bins)
ax5.hist(Vkick_dist_uniform, bins=bins)
ax6.hist(Vkick_dist_maxwellian, bins=bins)
ax4.hist(Vkick_dist_beniamini2[successful_binaries_log_ben2_ben2], bins=bins)
ax5.hist(Vkick_dist_uniform[successful_binaries_uni_uni_uni], bins=bins)
ax6.hist(Vkick_dist_maxwellian[successful_binaries_uni_power_maxwell], bins=bins)
ax4.set_xlabel('Vkick: Apre: Log Mhe: ben2 Vkick: ben2', fontdict=font)
ax5.set_xlabel('Vkick: Apre: uni Mhe: uni Vkick: uni', fontdict=font)
ax6.set_xlabel('Vkick: Apre: uni Mhe: power Vkick: maxwell', fontdict=font)

ax7.hist(Apre_dist_log, bins=bins)
ax8.hist(Apre_dist_uniform, bins=bins)
ax9.hist(Apre_dist_uniform, bins=bins)
ax7.hist(Apre_dist_log[successful_binaries_log_ben2_ben2], bins=bins)
ax8.hist(Apre_dist_uniform[successful_binaries_uni_uni_uni], bins=bins)
ax9.hist(Apre_dist_uniform[successful_binaries_uni_power_maxwell], bins=bins)
ax7.set_xlabel('Apre: Apre: Log Mhe: ben2 Vkick: ben2', fontdict=font)
ax8.set_xlabel('Apre: Apre: uni Mhe: uni Vkick: uni', fontdict=font)
ax9.set_xlabel('Apre: Apre: uni Mhe: power Vkick: maxwell', fontdict=font)
plot.suptitle('Number of Succesful Binaries Post Supernova for a Variety of Init Cond')
plot.show()

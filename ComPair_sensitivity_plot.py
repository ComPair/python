import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii

a=ascii.read("digitized_alex_sensitivities.dat",names=['eng','sensit'])
l=ascii.read("differential_flux_sensitivity_p8r2_source_v6_all_10yr_zmax100_n10.0_e1.50_ts25_000_090.txt",names=['emin','emax','e2diff','tmp'])

#print a['energy']
energy=a['eng']
sens=a['sensit']

erg2mev=624151.
lateng=(l["emin"]+l["emax"])/2.

plt.figure()
#LAT
plt.plot(lateng,l["e2diff"]*erg2mev,color='magenta',lw=2)
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlim([1e-2,1e6])
plt.gca().set_ylim([1e-8,1e-2])
plt.gca().set_xlabel('Energy (MeV)')
plt.gca().set_ylabel(r'Sensitivity $\times\ E^2$ [$\gamma$ MeV s$^{-1}$ cm$^{-2}$]')
plt.annotate('Fermi-LAT', xy=(5e2,2e-6),xycoords='data',fontsize=12,color='magenta')

#EGRET
ind=np.arange(69,74,1)
plt.plot(energy[ind],sens[ind],color='blue',lw=2)
plt.annotate('EGRET', xy=(1e2,1e-4),xycoords='data',fontsize=12,color='blue')

#SPI
ind=np.arange(20,46)
plt.plot(energy[ind],sens[ind],color='green',lw=2)
plt.annotate('SPI', xy=(6e-2,1e-4),xycoords='data',fontsize=12,color='green')

#COMPTEL
comptel_energy=[0.73295844,0.8483429,1.617075,5.057877,16.895761,29.717747]
comptel_sens=[6.566103E-4,3.6115389E-4,1.4393721E-4,1.6548172E-4,2.36875E-4,3.390693E-4]
plt.plot(comptel_energy,comptel_sens,color='orange',lw=2)
plt.annotate('COMPTEL', xy=(5,5e-4),xycoords='data',fontsize=12,color='orange')

#NuSTAR
ind=np.arange(84,147)
plt.plot(energy[ind]*1e-3,sens[ind]*(energy[ind]/1e3)**2*1e3,color='purple',lw=2)
plt.annotate('NuSTAR', xy=(0.1,3e-8),xycoords='data',fontsize=12,color='purple')

#ComPair
compair_eng=np.array([0.316,1,3.16,10,31.6,100,316.])
tracked=np.array([  1.54877155e-05,   4.84546681e-06,   5.28735667e-06, 6.53265846e-05, 0, 0, 0])
untracked=np.array([  2.49626245e-06,   1.82264874e-06,   1.54100276e-05, 9.59603201e-05, 0, 0, 0])
pair=np.array([ 0, 0,   5.62236032e-05, 3.19254897e-05,   1.71183233e-05,   1.61203804e-05, 2.19339e-05])
w=tracked > 0
plt.plot(compair_eng[w],tracked[w],color='black',lw=3)
w=pair > 0
plt.plot(compair_eng[w],pair[w],'r--',color='black',lw=3)
plt.annotate('ComPair', xy=(1,1e-6),xycoords='data',fontsize=18)

plt.savefig('new_sensitivity.pdf')

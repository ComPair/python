#!/usr/bin/env python
"""
------------------------------------------------------------------------

Scripts to make plots for the ComPair MidEx Proposal Section D

------------------------------------------------------------------------
"""

import numpy as np
import matplotlib.pylab as plot
from astropy.io import ascii
from astropy.io import fits
from scipy import interpolate


def loginterpol(x,y,x1):

	f=interpolate.interp1d(np.log10(x),np.log10(y),bounds_error=False,fill_value="extrapolate",kind='linear')
	y1=10**f(np.log10(x1))

	return y1

def agn_plot():

	fig = plot.figure()
	yrange = [1e-14,1e-8]#[1e-8,1e-3]#[1e-13, 1e2]
	xrange = [1e-11,1e6]
	plot.fill_between([0.2,10e3],[yrange[1],yrange[1]],[yrange[0],yrange[0]],facecolor='yellow',interpolate=True,color='yellow',alpha=0.5)
	plot.annotate('AMEGO',xy=(0.7,2e-9),xycoords='data',fontsize=26,color='black')
	plot.annotate('PMN J0641-0320',xy=(1e-10,3e-9),xycoords='data',fontsize=16,color='black')
	plot.annotate('z=1.196',xy=(1e-10,1e-9),xycoords='data',fontsize=16,color='black')

	## AGN
	a=ascii.read("data/marco_AGN_data_points_fig_6.txt",names=['logfreq','logflux'])
	logfreq=a['logfreq']
	logflux=a['logflux']
	h=6.6261e-27 #erg s
	erg2mev=624151.
	agn_energy=10**logfreq*h*erg2mev #Hz * erg s
	agn_flux=10**logflux#*erg2mev #erg cm-2 s-1
	agn_energy=np.append(agn_energy,2.6e4)
	agn_flux=np.append(agn_flux,1e-14)

	b=ascii.read("data/0641_0320nustar_torus.dat",names=['x','y','z','tot','f','fbump','ftot','syn','ssc','ext1c'])
	logfreq=b['x']
	lognufnu=b['ftot']
	agn_energy2=10**logfreq*h*erg2mev 
	agn_flux2=10**lognufnu#*erg2mev
	ssc=10**b['ssc']#*erg2mev
	syn=10**b['syn']#*erg2mev
	fbump=10**b['fbump']#*erg2mev
	f=10**b['f']#*erg2mev

	plot.plot(agn_energy,agn_flux,color='darkblue',lw=2)
	plot.plot(agn_energy2,agn_flux2,color='cornflowerblue',lw=2)
	#plot.plot(agn_energy2,fbump,'r--',color='lightblue',lw=2)
	#plot.plot(agn_energy2,ssc,'r--',color='cyan',lw=2)

	xrt=ascii.read("data/0641_0320xrt.dat",names=['logfreq','logflux'])
	xrt_energy=(10**xrt['logfreq'])*h*erg2mev
	xrt_flux=(10**xrt['logflux'])#*erg2mev
	#plot.scatter(xrt_energy,xrt_flux,color='blue')

	opt=ascii.read("data/0641_0320radio_optical.dat",names=['logfreq','logflux'])
	opt_energy=(10**opt['logfreq'])*h*erg2mev
	opt_flux=(10**opt['logflux'])#*erg2mev
	plot.scatter(opt_energy,opt_flux,color='blue')

	nustar=ascii.read("data/0641_0320nustar_ajello_fab.dat",names=['logfreq','logflux','logflux_yerr0','logflux_yerr1'])
	ns_energy=(10**nustar['logfreq'])*h*erg2mev
	ns_flux=(10**nustar['logflux'])#*erg2mev
	plot.scatter(ns_energy,ns_flux,color='blue')

	lat=ascii.read("data/LAT_spec_NuSTARobs2.txt",names=['ener','ed_ener','eu_ener','flux','ed_flux','eu_flux','ulim_flux','TS','Npred'])
	plot.scatter(lat['ener'][0:3],lat['flux'][0:3]/erg2mev)
	plot.errorbar(lat['ener'][0:3],lat['flux'][0:3]/erg2mev,xerr=[lat['ed_ener'][0:3],lat['eu_ener'][0:3]],yerr=[lat['ed_flux'][0:3]/erg2mev,lat['eu_flux'][0:3]/erg2mev],capsize=0,fmt="none")

	plot.xscale('log')
	plot.yscale('log')
	plot.ylim(yrange)
	plot.xlim(xrange)
	plot.xlabel(r'Energy (MeV)')
	plot.ylabel(r'$\nu$ F$_{\nu}$ (erg cm$^{-2}$ s$^{-1}$)')

	plot.savefig('../plots/agn_plot.eps', bbox_inches='tight')
	plot.savefig('../plots/agn_plot.png', bbox_inches='tight')
	plot.show()
	plot.close()

def magentar_plot():

	xmm=ascii.read('data/xmm_SED_paper.dat',names=['energy','err_energy','flux','err_flux'])
	integral=ascii.read('data/integral_SED_paper.dat',names=['energy','err_energy','flux','err_flux'])
	comptel=ascii.read('data/comptel_SED_paper.dat',names=['energy','err_energy','flux','err_flux'])
	lat={'energy_min':[1e5,1e5,1e6],'energy_max':[1e7,1e6,1e7],'flux':[1.6e-4,2.1e-4,8.13e-5]}
	#magic={'energy_min':[200e3,775e3],'energy_max':[774e3,3000e3],'flux':[6.29e-5,8.99e-5]}
	models=ascii.read('data/magnetar_4U0142+614.dat',names=['energy','flux'])

	fig=plot.figure()
	xrange=[1e-4,1e5]
	yrange=[5e-13,3e-10]
	s=1e-8

	plot.fill_between([0.2,10e3],[yrange[1],yrange[1]],[yrange[0],yrange[0]],facecolor='yellow',interpolate=True,color='yellow',alpha=0.5)
	plot.annotate('AMEGO',xy=(1.5,1.7e-10),xycoords='data',fontsize=26,color='black')

	plot.errorbar(xmm['energy']*1e-3,xmm['flux']*s,xerr=xmm['err_energy']*1e-3,yerr=xmm['err_flux']*s,capsize=0, fmt='.',color='blue')
	plot.errorbar(integral['energy']*1e-3,integral['flux']*s,xerr=integral['err_energy']*1e-3,yerr=integral['err_flux']*s,capsize=0, fmt='.',color='blue')
	plot.errorbar(comptel['energy']*1e-3,comptel['flux']*s,xerr=comptel['err_energy']*1e-3,yerr=0.5*comptel['flux']*s,capsize=0, fmt='.',color='blue',uplims=True)
	
	m=np.array(10**((np.log10(lat['energy_min'])+np.log10(lat['energy_max']))/2))*1e-3
	lowe=m-np.array(lat['energy_min'])*1e-3
	highe=np.array(lat['energy_max'])*1e-3-m
	plot.errorbar(m,np.array(lat['flux'])*s,xerr=[lowe,highe],yerr=0.3*np.array(lat['flux'])*s,capsize=0, fmt='.',color='blue',uplims=True)

	# m=np.array(10**((np.log10(magic['energy_min'])+np.log10(magic['energy_max']))/2))
	# lowe=m-np.array(magic['energy_min'])
	# highe=np.array(magic['energy_max'])-m
	# #plot.errorbar(m,np.array(magic['flux'])*s,xerr=[lowe,highe],yerr=0.2*np.array(magic['flux'])*s,capsize=0, fmt='.',color='blue',uplims=True)

	erg2kev=624151.*1e3
	color=['red','magenta','black','black']
	lines=['','','','r--']

	num=[0,18,38,67,len(models['energy'])]
	for i in range(0,4):
		n0=num[i]
		n1=num[i+1]
		ind=np.arange(n0,n1)
		x=np.array(models['energy'])*1e-3
		y=np.array(models['flux'])/erg2kev
#		n2=(n1-n0)/2+n0
		w=np.where(y[ind]-max(y[ind]) == 0)
		n2=ind[w[0][0]]
		x0=loginterpol(y[n0:n2],x[n0:n2],yrange[0])
		x1=loginterpol(y[n2:n1],x[n2:n1],yrange[0])
		y=np.append(np.append(yrange[0],y[ind]),yrange[0])
		x=np.append(np.append(x0,x[ind]),x1)
#		y=np.append(yrange[0],y)
#		x=np.append(x0,x)
		plot.plot(x,y,lines[i],color=color[i])


	#y=np.append(np.append(yrange[0],models['flux'][ind]/erg2kev),yrange[0])
	#x=loginterpol(models['flux'][ind[0:9]]/erg2kev,models['energy'][ind[0:9]*1e-3,y)
	#plot.plot(x,y,color='red')
	#print x[0],x[19]
	#print y[0],y[19]

	# ind=np.arange(18,38)
	# plot.plot(models['energy'][ind]*1e-3,models['flux'][ind]/erg2kev,color='magenta')

	# ind=np.arange(38,67)
	# plot.plot(models['energy'][ind]*1e-3,models['flux'][ind]/erg2kev,color='black')

	# ind=np.arange(67,len(models['flux']))
	# plot.plot(models['energy'][ind]*1e-3,models['flux'][ind]/erg2kev,'r--',color='black')


	plot.xscale('log')
	plot.yscale('log')
	plot.ylim(yrange)
	plot.xlim(xrange)
	plot.xlabel(r'Energy (MeV)')
	plot.ylabel(r'Energy$^2 \times $ Flux (Energy) (erg cm$^{-2}$ s$^{-1}$)')
	plot.show()

def nova_plot():

	erg2mev=624151.

	fig=plot.figure()
	yrange = [1e-6,2e-4]
	xrange = [1e-1,1e5]
	plot.fill_between([0.2,10e3],[yrange[1],yrange[1]],[yrange[0],yrange[0]],facecolor='yellow',interpolate=True,color='yellow',alpha=0.5)
	plot.annotate('AMEGO',xy=(3,9e-5),xycoords='data',fontsize=26,color='black')

	lat=ascii.read("data/NMon2012.LAT.dat",names=['energy','en_low','en_high','flux','flux_err','tmp'])
	plot.scatter(lat['energy'],lat['flux']*erg2mev,color='red')
	plot.errorbar(lat['energy'],lat['flux']*erg2mev,xerr=[lat['en_low'],lat['en_high']],yerr=lat['flux_err']*erg2mev,ecolor='red',capsize=0,fmt='none')
	latul=ascii.read("data/NMon2012.LAT.limits.dat",names=['energy','en_low','en_high','flux','tmp1','tmp2','tmp3','tmp4'])
	plot.errorbar(latul['energy'],latul['flux']*erg2mev,xerr=[latul['en_low'],latul['en_high']],yerr=0.5*latul['flux']*erg2mev,uplims=True,ecolor='red',capsize=0,fmt='none')
	plot.scatter(latul['energy'],latul['flux']*erg2mev,color='red')

	leptonic=ascii.read("data/sp-NMon12-IC-best-fit-1MeV-30GeV.txt",names=['energy','flux'],data_start=1)
	hadronic=ascii.read("data/sp-NMon12-pi0-and-secondaries.txt",names=['energy','flux1','flux2'],data_start=1)	

	plot.plot(leptonic['energy'],leptonic['flux']*erg2mev,'r--',color='black',lw=2,label='Leptonic')
	plot.plot(hadronic['energy'],hadronic['flux2']*erg2mev,color='black',lw=2,label='Hadronic+Secondary Leptons')

	plot.legend(loc='upper right',fontsize='small',frameon=False,framealpha=0.5)
	plot.xscale('log')
	plot.yscale('log')
	plot.ylim(yrange)
	plot.xlim(xrange)
	plot.xlabel(r'Energy (MeV)')
	plot.ylabel(r'Energy$^2 \times $ Flux (Energy) (erg cm$^{-2}$ s$^{-1}$)')
	plot.title('Nova V339 Del 2013')
	plot.savefig('Nova_SED.png', bbox_inches='tight')
	plot.savefig('Nova_SED.eps', bbox_inches='tight')
	plot.show()
	plot.close()

def SMBH_mass(save=False):

	fig=plot.figure()

	#a=ascii.read('data/BHmass_dist.dat',names=['mass','N'],data_start=1)
	#mass=(np.round(a['mass']*100.)/100.)
	#N=np.array(np.round(a['N']*100),dtype=np.int64)

	high=ascii.read('data/BH_mass_High_z.txt',names=['mass'])
	low=ascii.read('data/BH_mass_Low_z.txt',names=['mass'])

	loghigh=np.log10(high['mass'])
	loglow=np.log10(low['mass'])
	#ind1=np.arange(0,13,1)
	#ind2=np.arange(13,len(mass),1)

	#m1=np.repeat(mass[ind1],N[ind1])
	#w1=np.repeat(np.repeat(1./max(N[ind1]),len(ind1)),N[ind1])
	#m2=np.repeat(mass[ind2],N[ind2])
	#w2=np.repeat(np.repeat(1./max(N[ind2]),len(ind2)),N[ind2])

	low_bin=np.logspace(np.min(loglow),np.max(loglow),num=14)
	plot.hist(low['mass'],bins=low_bin,color='blue',weights=np.repeat(1./28,len(low)))
	high_bin=np.logspace(np.min(loghigh),np.max(loghigh),num=10)
	plot.hist(high['mass'],bins=high_bin,color='red',alpha=0.7,weights=np.repeat(1./28,len(high)))

	plot.annotate('Low Redshift (z < 3) Blazars',xy=(1.5e9,0.8),xycoords='data',fontsize=14,color='blue')
	plot.annotate('High Redshift (z > 3) Blazars',xy=(1.5e9,0.5),xycoords='data',fontsize=14,color='red')

	plot.xlim([5e7,5e10])
	plot.ylim([0,1.05])
	plot.xscale('log')
#	plot.yscale('log')
	plot.xlabel(r'Black Hole Mass (M$_{\odot}$)')
	plot.ylabel('Fraction of Known Blazars')
#	plot.title('Supermassive Black Hole Mass Evolution')

	if save:
		plot.savefig('SMBH_mass.png', bbox_inches='tight')
		plot.savefig('SMBH_mass.eps', bbox_inches='tight')
	else:
		plot.show()

	return

def UNIDplot(save=False):

	import FigureOfMeritPlotter
	from scipy import interpolate

	data=FigureOfMeritPlotter.parseEventAnalysisLogs(directory='../Simulations/PerformancePlotTraFiles/',silent=True)
	ComPair=FigureOfMeritPlotter.plotAllSourceSensitivities(data,angleSelection=1.0,doplot=False)

	latcat=fits.open('data/gll_psc_v14.fit')
	data=latcat[1].data
	#print latcat.info()
	#print latcat[1].columns
	wunid=np.where((data.field('ASSOC1') == ' ') & (data.field('ASSOC2') == ' '))
	Tot=len(wunid[-1])
	print Tot
	w=np.where((data.field('ASSOC1') == ' ') & (data.field('ASSOC2') == ' ') & (data.field('SpectrumType')=='PowerLaw'))
	tot=len(w[-1])
	print tot
	flux=data[w].field('Flux1000')
	emin=1
	emax=1000.
	erg2mev=624151.

	index=data[w].field('Spectral_Index')

	fig=plot.figure()

	yrange = [1e-15,4e-15]
	plot.fill_between([0.2,10e3],[yrange[1],yrange[1]],[yrange[0],yrange[0]],facecolor='yellow',interpolate=True,color='yellow',alpha=0.5)
	plot.annotate('AMEGO',xy=(4,1.5e-15),xycoords='data',fontsize=26,color='black')

	yrange = [5e-15,2.5e-14]
	plot.fill_between([100,300e3],[yrange[1],yrange[1]],[yrange[0],yrange[0]],facecolor='magenta',interpolate=True,color='magenta',alpha=0.5)
	plot.annotate('LAT',xy=(400,8e-15),xycoords='data',fontsize=26,color='black')


	logenergy=np.arange(-1,6,0.1)
	energy=10**logenergy
  	tck=interpolate.splrep(np.log10(ComPair[0][1:]),np.log10(ComPair[1][1:]/erg2mev),s=0)
  	w=np.where((energy > 0.5) & (energy < 500))
	compflux=10**interpolate.splev(np.log10(energy[w]),tck,der=0)

	N_nodet=0.
	N_det=0.
	N_highpeak=0.
	for i in range(len(index)):
		norm=flux[i]*(1.-index[i])/(emax**(1-index[i])-emin**(1-index[i]))
		f=norm*(energy/1e3)**(-index[i])*1e3#(energy**(2.-index[i])/(2-index[i]))
		e2f=f*energy**2/erg2mev**2
		if (index[i] < 2):
			color='grey'
			N_nodet+=1
		if (index[i] > 2) & (max(e2f[w]-compflux) > 0):
			color='red'
			N_det+=1
		if (max(e2f[w]-compflux) < 0) & (index[i] > 2):
			color='blue'
			N_highpeak+=1
		plot.plot(energy,e2f,'r:',color=color)

	print 'N_nodet: ',N_nodet, N_nodet/tot, N_nodet/Tot
	print 'N_det: ',N_det, N_det/tot, N_det/Tot
	print 'N_highpeak: ',N_highpeak, N_highpeak/tot, N_highpeak/Tot

	plot.annotate(r"AMEGO Detectable (Peak < LAT)",xy=(21,5e-9),xycoords='data',fontsize=14,color='red')
	plot.annotate(r'AMEGO Non-Detectable (Peak < LAT)',xy=(9,2e-9),xycoords='data',fontsize=14,color='blue')
	plot.annotate(r'AMEGO Non-Detectable (Peak $\geq$ LAT)',xy=(8,8e-10),xycoords='data',fontsize=14,color='grey')


	#print data[w].field('Spectral_Index')

	#plot.plot(ComPair[0][1:],ComPair[1][1:]/erg2mev,color='blue',lw=3)

	plot.xlim([1e-1,1e5])
	plot.ylim([1e-15,1e-8])
	plot.xscale('log')
	plot.yscale('log')
	plot.xlabel(r'Energy (MeV)')
	plot.ylabel(r'Energy$^2 \times $ Flux (Energy) (erg cm$^{-2}$ s$^{-1}$)')
	plot.title('Fermi-LAT Unidentified Sources in the MeV Band')
	if save:
		plot.savefig('UNID_SED.png', bbox_inches='tight')
		plot.savefig('UNID_SED.eps', bbox_inches='tight')
	else:
		plot.show()
	plot.close()
	latcat.close()
	return

def FillingTheGap(save=False):

	## ComPair
	fig = plot.figure()
	yrange = [1e-6,4e-3]#[1e-13, 1e2]
	xrange = [1e-3,1e7]
	plot.fill_between([0.2,10e3],[yrange[1],yrange[1]],[yrange[0],yrange[0]],facecolor='yellow',interpolate=True,color='yellow',alpha=0.5)
	plot.annotate('AMEGO',xy=(4,1.5e-3),xycoords='data',fontsize=26,color='black')

	## AGN
	a=ascii.read("data/marco_AGN_data_points_fig_6.txt",names=['logfreq','logflux'])
	logfreq=a['logfreq']
	logflux=a['logflux']
	h=6.6261e-27 #erg s
	erg2mev=624151.
	agn_energy=10**logfreq*h*erg2mev #Hz * erg s
	arbfact=5
	agn_flux=10**logflux*erg2mev*arbfact #erg cm-2 s-1

	i = np.where(agn_energy < 0.1)
	plot.plot(agn_energy[i],agn_flux[i],color='navy',lw=2)
	i = np.where((agn_energy > 0.1) & (agn_energy < 200))
	plot.plot(agn_energy[i],agn_flux[i],'r--',color='navy',lw=2)	
	i=np.where(agn_energy > 200)
	plot.plot(agn_energy[i],agn_flux[i],color='navy',lw=2)

	b=ascii.read("data/0641_0320nustar_torus.dat",names=['x','y','z','tot','f','fbump','ftot','syn','ssc','ext1c'])
	logfreq=b['x']
	lognufnu=b['ftot']
	agn_energy2=10**logfreq*h*erg2mev 
	agn_flux2=10**lognufnu*erg2mev*arbfact

	i = np.where(agn_energy2 < 0.1)
	plot.plot(agn_energy2[i],agn_flux2[i],color='cornflowerblue',lw=2)
	i = np.where((agn_energy2 > 0.1) & (agn_energy2 < 200))
	plot.plot(agn_energy2[i],agn_flux2[i],'r--',color='cornflowerblue',lw=2)	
	i=np.where(agn_energy2 > 200)
	plot.plot(agn_energy2[i],agn_flux2[i],color='cornflowerblue',lw=2)

	xrt=ascii.read("data/0641_0320xrt.dat",names=['logfreq','logflux'])
	xrt_energy=(10**xrt['logfreq'])*h*erg2mev
	xrt_flux=(10**xrt['logflux'])*erg2mev*arbfact
	plot.scatter(xrt_energy,xrt_flux,color='blue')

	nustar=ascii.read("data/0641_0320nustar_ajello_fab.dat",names=['logfreq','logflux','logflux_yerr0','logflux_yerr1'])
	ns_energy=(10**nustar['logfreq'])*h*erg2mev
	ns_flux=(10**nustar['logflux'])*erg2mev*arbfact
	plot.scatter(ns_energy,ns_flux,color='blue')

	lat=ascii.read("data/LAT_spec_NuSTARobs2.txt",names=['ener','ed_ener','eu_ener','flux','ed_flux','eu_flux','ulim_flux','TS','Npred'])
	plot.scatter(lat['ener'],lat['flux']*arbfact)
	plot.errorbar(lat['ener'][0:3],lat['flux'][0:3]*arbfact,xerr=[lat['ed_ener'][0:3],lat['eu_ener'][0:3]],yerr=[lat['ed_flux'][0:3]*arbfact,lat['eu_flux'][0:3]*arbfact],capsize=0,fmt="none")

	## Pulsar example

	pulsar_eng=np.array([0.012943256,0.018285165,0.031053474,0.05153211,0.08552302,0.21973862,137.03448,237.55414])
	pulsar_flux=np.array([1.7420283E-5,2.2255874E-5,3.0082629E-5,3.842357E-5,5.0966246E-5,7.149577E-5,1.4489453E-5,6.674534E-6])

	pulsar_eng_ul=np.array([421.64273,748.32324])
	pulsar_flux_ul=np.array([4.0049385E-6,2.314023E-6])

	xs = np.linspace(np.min(np.log10(1e-3)), np.max(np.log10(1e3)), 300)

	pulsar_Energy = 10**xs

	e0=[0.100518,0.100518*0.4,0.100518] # norm energy
	k=[1.574e-2,1.574e-2*2,1.574e-2*5] # normalization
	gamma=[-1.233,-1.18,-1.2] # 
	ec=[0.078,0.8,5e-4] #cutoff energy
	beta=[0.286,0.4,0.18] # cutoff slope
	arbfact=1

	color=['darkgreen','sage']

	for j in range(0,2):

		flux=k[j]*(pulsar_Energy/e0[j])**gamma[j]*np.exp(-(pulsar_Energy/ec[j])**beta[j])
		pulsar_Flux=flux*pulsar_Energy**2

		i = np.where(pulsar_Energy < 0.2)
		plot.plot(pulsar_Energy[i],pulsar_Flux[i]*arbfact,color=color[j],lw=2)
		i = np.where((pulsar_Energy > 0.2) & (pulsar_Energy < 100))
		plot.plot(pulsar_Energy[i],pulsar_Flux[i]*arbfact,'r--',color=color[j],lw=2)	
		i=np.where(pulsar_Energy > 100)
		plot.plot(pulsar_Energy[i],pulsar_Flux[i]*arbfact,color=color[j],lw=2)

	plot.scatter(pulsar_eng,pulsar_flux*arbfact,color='green')

	errfrac=np.concatenate((np.repeat(0.1,6),(0.4,0.6)),axis=0)
	plot.errorbar(pulsar_eng,pulsar_flux*arbfact,yerr=errfrac*pulsar_flux*arbfact,color='green',ecolor='green',capsize=0,fmt="none")
	plot.errorbar(pulsar_eng_ul,pulsar_flux_ul*arbfact,yerr=0.5*pulsar_flux_ul*arbfact,color='green',ecolor='green',capsize=0,uplims=True,fmt="none")
	plot.scatter(pulsar_eng_ul,pulsar_flux_ul*arbfact,color='green')

	## Nova

	arbfact=0.5
	#osse=ascii.read("data/NVel1999.OSSE.dat",names=['energy','en_low','en_high','flux','flux_err'])
	#plot.scatter(osse['energy'],osse['flux']*erg2mev*arbfact,color='red')
	#plot.errorbar(osse['energy'],osse['flux']*erg2mev*arbfact,xerr=[osse['en_low'],osse['en_high']],yerr=osse['flux_err']*erg2mev*arbfact,ecolor='red',capsize=0,fmt='none')

	lat=ascii.read("data/NMon2012.LAT.dat",names=['energy','en_low','en_high','flux','flux_err','tmp'])
	plot.scatter(lat['energy'],lat['flux']*erg2mev*arbfact,color='red')
	plot.errorbar(lat['energy'],lat['flux']*erg2mev*arbfact,xerr=[lat['en_low'],lat['en_high']],yerr=lat['flux_err']*erg2mev*arbfact,ecolor='red',capsize=0,fmt='none')
	latul=ascii.read("data/NMon2012.LAT.limits.dat",names=['energy','en_low','en_high','flux','tmp1','tmp2','tmp3','tmp4'])
	plot.errorbar(latul['energy'],latul['flux']*erg2mev*arbfact,xerr=[latul['en_low'],latul['en_high']],yerr=0.5*latul['flux']*erg2mev*arbfact,uplims=True,ecolor='red',capsize=0,fmt='none')
	plot.scatter(latul['energy'],latul['flux']*erg2mev*arbfact,color='red')

	#models=ascii.read("data/data-NovaMon2012.txt",names=['energy','leptonic','hadronic'],data_start=1)
	#mo=['leptonic','hadronic']
	colors=['orangered','coral','darkred']
	leptonic=ascii.read("data/sp-NMon12-IC-best-fit-1MeV-30GeV.txt",names=['energy','flux'],data_start=1)
	hadronic=ascii.read("data/sp-NMon12-pi0-and-secondaries.txt",names=['energy','flux1','flux2'],data_start=1)	
	for j in range(0,2):
		if (j == 0):
			energy=leptonic['energy']
			flux=leptonic['flux']
		if (j == 1):
			energy=hadronic['energy']			
			flux=hadronic['flux2']
		if (j == 2):
			flux=hadronic['flux1']
		i=np.where(energy < 0.2)
		plot.plot(energy[i],flux[i]*erg2mev*arbfact,color=colors[j],lw=2)
		i=np.where((energy > 0.2) & (energy <100 ))
		plot.plot(energy[i],flux[i]*erg2mev*arbfact,'r--',color=colors[j],lw=2)
		i=np.where(energy > 100)
		plot.plot(energy[i],flux[i]*erg2mev*arbfact,color=colors[j],lw=2)


	# PWNe
	# From Torres et al. JHEAP, 2014, 1, 31
	# G54.1+0.3

	# arbfact=1#e-3
	# pwne_model_eng=np.array([860.1604,1866.5002,3879.583,7087.079,10898.913,17498.09,28093.0,34838.23,41383.043,299664.66,6644004.5,3.5596056E7,1.29464152E8,4.91573856E8,3.71615181E9,2.36500337E10,8.2392531E10,2.52272067E11,5.47416343E11,8.7887118E11])*1e-6
	# pwne_model_flux=np.array([1.2608298E-11,9.389613E-12,5.4490882E-12,2.8233574E-12,1.1928516E-12,4.0173412E-13,1.7362255E-13,1.3226478E-13,1.4482099E-13,3.1305714E-13,7.409742E-13,8.684221E-13,8.2992185E-13,6.3223027E-13,2.9917822E-13,8.0318205E-14,1.925147E-14,4.119833E-15,8.816484E-16,1.0E-16])*erg2mev*arbfact

	# xs = np.linspace(np.min(np.log10(1e-3)), np.max(np.log10(5e3)), 300)
	# pwne_Energy=10**xs

	# tck=interpolate.splrep(np.log10(pwne_model_eng),np.log10(pwne_model_flux),s=0)
	# pwne_Flux=10**interpolate.splev(np.log10(pwne_Energy),tck,der=0)

	# pwne_data_eng=np.array([1.7125564E7,5.4741636E7,1.74980896E8,2.46901296E8,3.9639744E8,6.3641197E8,1.11359936E9,1.64041331E9,2.52272051E9,4.0502016E9])*1e-6
	# pwne_data_upper_flux=np.array([4.2462813E-12,4.7560116E-12,1.5817023E-11,3.7911818E-12,1.2202063E-12,2.2001422E-12,1.8772538E-12,2.6377E-12,1.1144125E-12,2.1508195E-12])*erg2mev*arbfact
	# pwne_data_flux=np.array([0,0,0,1.4628772E-12,3.2757992E-13,8.489537E-13,7.579663E-13,1.2202063E-12,2.3313897E-13,6.9224937E-13])*erg2mev*arbfact
	# pwne_data_lower_flux=np.array([8.883369E-13,1.1144125E-12,3.3089764E-12,3.5867788E-13,3.242955E-14,3.5867788E-13,2.4395433E-13,5.039727E-13,1.7188176E-14,1.2356735E-13])*erg2mev*arbfact

	# i = np.where(pwne_Energy < 0.5)
	# plot.plot(pwne_Energy[i],pwne_Flux[i],color='red',lw=2)
	# i = np.where((pwne_Energy > 0.5) & (pwne_Energy < 200))
	# plot.plot(pwne_Energy[i],pwne_Flux[i],'r--',color='red',lw=2)	
	# i=np.where(pwne_Energy > 200)
	# plot.plot(pwne_Energy[i],pwne_Flux[i],color='red',lw=2)

	# plot.errorbar(pwne_data_eng[3:],pwne_data_flux[3:],yerr=[pwne_data_flux[3:]-pwne_data_lower_flux[3:],pwne_data_upper_flux[3:]-pwne_data_flux[3:]],color='tomato',fmt='o',ecolor='tomato',capsize=0,lw=2)
	# plot.errorbar(pwne_data_eng[0:3],pwne_data_upper_flux[0:3],yerr=pwne_data_upper_flux[0:3]-pwne_data_lower_flux[0:3],color='tomato',fmt='o',ecolor='tomato',uplims=True,lw=2)

	# plot stuff
	plot.xscale('log')
	plot.yscale('log')
	plot.ylim(yrange)
	plot.xlim(xrange)
	plot.xlabel(r'Energy (MeV)')
	plot.ylabel(r'Energy$^2 \times $ Flux (Energy) (arbitrarily scaled)')

	im1=plot.imread('data/AGN_UnifiedModel.jpg')
	newax1=fig.add_axes([0.73,0.65,0.15,0.15],anchor='NE')
	newax1.imshow(im1)
	for axis in ['top','bottom','right','left']:
			newax1.spines[axis].set_linewidth(4)
			newax1.spines[axis].set_color('blue')
	newax1.set_xticks([])
	newax1.set_yticks([])
#	newax1.axis('off')
	plot.title('Jets',color='blue',fontsize=12)

	im2=plot.imread('data/pulsar.jpg')
	newax2=fig.add_axes([0.73,0.4,0.15,0.16],anchor='NE')
	newax2.imshow(im2)
	for axis in ['top','bottom','right','left']:
			newax2.spines[axis].set_linewidth(4)
			newax2.spines[axis].set_color('green')
	newax2.set_xticks([])
	newax2.set_yticks([])
#	newax2.axis('off')
	plot.title('Compact Objects',color='green',fontsize=12)

	im3=plot.imread('data/Classical_Nova_Final.jpg')
	newax3=fig.add_axes([0.73,0.18,0.15,0.15],anchor='NE')
	newax3.imshow(im3)
	for axis in ['top','bottom','right','left']:
			newax3.spines[axis].set_linewidth(4)
			newax3.spines[axis].set_color('red')
	newax3.set_xticks([])
	newax3.set_yticks([])

#	newax3.axis('off')
	plot.title('Shocks',color='red',fontsize=12)

	if save:
		#plot.savefig('SED_science_themes.eps', bbox_inches='tight')
		plot.savefig('SED_science_themes.pdf', bbox_inches='tight')
		plot.savefig('SED_science_themes.png', bbox_inches='tight')
	plot.show()
	plot.close()

	## to do
	### add colored box around each image
	### add labels for jets, compact objects, shocks + maybe source type

	return
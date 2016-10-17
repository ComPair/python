#!/usr/bin/env python
"""
------------------------------------------------------------------------

Scripts to make plots for the ComPair MidEx Proposal Section D

------------------------------------------------------------------------
"""

import numpy
import matplotlib.pylab as plot
from astropy.io import ascii

def FillingTheGap(save=False):

	from scipy import interpolate

	## ComPair
	fig = plot.figure()
	yrange = [1e-7,1e-3]#[1e-13, 1e2]
	xrange = [1e-3,1e7]
	plot.fill_between([0.5,500],[yrange[1],yrange[1]],[yrange[0],yrange[0]],facecolor='yellow',interpolate=True,color='yellow',alpha=0.5)
	plot.annotate('ComPair',xy=(2,1.5e-7),xycoords='data',fontsize=18,color='black')

	## AGN
	a=ascii.read("data/marco_AGN_data_points_fig_6.txt",names=['logfreq','logflux'])
	logfreq=a['logfreq']
	logflux=a['logflux']
	h=6.6261e-27 #erg s
	erg2mev=624151.
	agn_energy=10**logfreq*h*erg2mev #Hz * erg s
	arbfact=4#e3
	agn_flux=10**logflux*erg2mev*arbfact #erg cm-2 s-1

	i = numpy.where(agn_energy < 0.1)
	plot.plot(agn_energy[i],agn_flux[i],color='blue',lw=2)
	i = numpy.where((agn_energy > 0.1) & (agn_energy < 200))
	plot.plot(agn_energy[i],agn_flux[i],'r--',color='blue',lw=2)	
	i=numpy.where(agn_energy > 200)
	plot.plot(agn_energy[i],agn_flux[i],color='blue',lw=2)

	b=ascii.read("data/0641_0320nustar_torus.dat",names=['x','y','z','tot','f','fbump','ftot','syn','ssc','ext1c'])
	logfreq=b['x']
	lognufnu=b['ftot']
	agn_energy2=10**logfreq*h*erg2mev 
	agn_flux2=10**lognufnu*erg2mev*arbfact

	i = numpy.where(agn_energy2 < 0.1)
	plot.plot(agn_energy2[i],agn_flux2[i],color='cornflowerblue',lw=2)
	i = numpy.where((agn_energy2 > 0.1) & (agn_energy2 < 200))
	plot.plot(agn_energy2[i],agn_flux2[i],'r--',color='cornflowerblue',lw=2)	
	i=numpy.where(agn_energy2 > 200)
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

	pulsar_eng=numpy.array([0.012943256,0.018285165,0.031053474,0.05153211,0.08552302,0.21973862,137.03448,237.55414])
	pulsar_flux=numpy.array([1.7420283E-5,2.2255874E-5,3.0082629E-5,3.842357E-5,5.0966246E-5,7.149577E-5,1.4489453E-5,6.674534E-6])

	pulsar_eng_ul=numpy.array([421.64273,748.32324])
	pulsar_flux_ul=numpy.array([4.0049385E-6,2.314023E-6])

	xs = numpy.linspace(numpy.min(numpy.log10(1e-3)), numpy.max(numpy.log10(1e3)), 300)

	pulsar_Energy = 10**xs

	e0=[0.100518,0.100518*0.4,0.100518] # norm energy
	k=[1.574e-2,1.574e-2*2,1.574e-2*5] # normalization
	gamma=[-1.233,-1.18,-1.2] # 
	ec=[0.078,0.8,5e-4] #cutoff energy
	beta=[0.286,0.4,0.18] # cutoff slope
	arbfact=0.5

	color=['purple','darkorchid','mediumpurple']

	for j in range(0,3):

		flux=k[j]*(pulsar_Energy/e0[j])**gamma[j]*numpy.exp(-(pulsar_Energy/ec[j])**beta[j])
		pulsar_Flux=flux*pulsar_Energy**2

		i = numpy.where(pulsar_Energy < 0.2)
		plot.plot(pulsar_Energy[i],pulsar_Flux[i]*arbfact,color=color[j],lw=2)
		i = numpy.where((pulsar_Energy > 0.2) & (pulsar_Energy < 100))
		plot.plot(pulsar_Energy[i],pulsar_Flux[i]*arbfact,'r--',color=color[j],lw=2)	
		i=numpy.where(pulsar_Energy > 100)
		plot.plot(pulsar_Energy[i],pulsar_Flux[i]*arbfact,color=color[j],lw=2)

	plot.scatter(pulsar_eng,pulsar_flux*arbfact,color='purple')

	errfrac=numpy.concatenate((numpy.repeat(0.1,6),(0.4,0.6)),axis=0)
	plot.errorbar(pulsar_eng,pulsar_flux*arbfact,yerr=errfrac*pulsar_flux*arbfact,color='purple',ecolor='purple',capsize=0,fmt="none")
	plot.errorbar(pulsar_eng_ul,pulsar_flux_ul*arbfact,yerr=0.8*pulsar_flux_ul*arbfact,color='purple',ecolor='purple',capsize=0,uplims=True,fmt="none")

	## Nova

	arbfact=0.1
	osse=ascii.read("data/NVel1999.OSSE.dat",names=['energy','en_low','en_high','flux','flux_err'])
	plot.scatter(osse['energy'],osse['flux']*erg2mev*arbfact,color='red')
	plot.errorbar(osse['energy'],osse['flux']*erg2mev*arbfact,xerr=[osse['en_low'],osse['en_high']],yerr=osse['flux_err']*erg2mev*arbfact,ecolor='red',capsize=0,fmt='none')

	lat=ascii.read("data/NMon2012.LAT.dat",names=['energy','en_low','en_high','flux','flux_err','tmp'])
	plot.scatter(lat['energy'],lat['flux']*erg2mev*arbfact,color='red')
	plot.errorbar(lat['energy'],lat['flux']*erg2mev*arbfact,xerr=[lat['en_low'],lat['en_high']],yerr=lat['flux_err']*erg2mev*arbfact,ecolor='red',capsize=0,fmt='none')
	latul=ascii.read("data/NMon2012.LAT.limits.dat",names=['energy','en_low','en_high','flux','tmp1','tmp2','tmp3','tmp4'])
	plot.errorbar(latul['energy'],latul['flux']*erg2mev*arbfact,xerr=[latul['en_low'],latul['en_high']],yerr=0.8*latul['flux']*erg2mev*arbfact,uplims=True,ecolor='red',capsize=0,fmt='none')
	plot.scatter(latul['energy'],latul['flux']*erg2mev*arbfact,color='red')


	# PWNe
	# From Torres et al. JHEAP, 2014, 1, 31
	# G54.1+0.3

	# arbfact=1#e-3
	# pwne_model_eng=numpy.array([860.1604,1866.5002,3879.583,7087.079,10898.913,17498.09,28093.0,34838.23,41383.043,299664.66,6644004.5,3.5596056E7,1.29464152E8,4.91573856E8,3.71615181E9,2.36500337E10,8.2392531E10,2.52272067E11,5.47416343E11,8.7887118E11])*1e-6
	# pwne_model_flux=numpy.array([1.2608298E-11,9.389613E-12,5.4490882E-12,2.8233574E-12,1.1928516E-12,4.0173412E-13,1.7362255E-13,1.3226478E-13,1.4482099E-13,3.1305714E-13,7.409742E-13,8.684221E-13,8.2992185E-13,6.3223027E-13,2.9917822E-13,8.0318205E-14,1.925147E-14,4.119833E-15,8.816484E-16,1.0E-16])*erg2mev*arbfact

	# xs = numpy.linspace(numpy.min(numpy.log10(1e-3)), numpy.max(numpy.log10(5e3)), 300)
	# pwne_Energy=10**xs

	# tck=interpolate.splrep(numpy.log10(pwne_model_eng),numpy.log10(pwne_model_flux),s=0)
	# pwne_Flux=10**interpolate.splev(numpy.log10(pwne_Energy),tck,der=0)

	# pwne_data_eng=numpy.array([1.7125564E7,5.4741636E7,1.74980896E8,2.46901296E8,3.9639744E8,6.3641197E8,1.11359936E9,1.64041331E9,2.52272051E9,4.0502016E9])*1e-6
	# pwne_data_upper_flux=numpy.array([4.2462813E-12,4.7560116E-12,1.5817023E-11,3.7911818E-12,1.2202063E-12,2.2001422E-12,1.8772538E-12,2.6377E-12,1.1144125E-12,2.1508195E-12])*erg2mev*arbfact
	# pwne_data_flux=numpy.array([0,0,0,1.4628772E-12,3.2757992E-13,8.489537E-13,7.579663E-13,1.2202063E-12,2.3313897E-13,6.9224937E-13])*erg2mev*arbfact
	# pwne_data_lower_flux=numpy.array([8.883369E-13,1.1144125E-12,3.3089764E-12,3.5867788E-13,3.242955E-14,3.5867788E-13,2.4395433E-13,5.039727E-13,1.7188176E-14,1.2356735E-13])*erg2mev*arbfact

	# i = numpy.where(pwne_Energy < 0.5)
	# plot.plot(pwne_Energy[i],pwne_Flux[i],color='red',lw=2)
	# i = numpy.where((pwne_Energy > 0.5) & (pwne_Energy < 200))
	# plot.plot(pwne_Energy[i],pwne_Flux[i],'r--',color='red',lw=2)	
	# i=numpy.where(pwne_Energy > 200)
	# plot.plot(pwne_Energy[i],pwne_Flux[i],color='red',lw=2)

	# plot.errorbar(pwne_data_eng[3:],pwne_data_flux[3:],yerr=[pwne_data_flux[3:]-pwne_data_lower_flux[3:],pwne_data_upper_flux[3:]-pwne_data_flux[3:]],color='tomato',fmt='o',ecolor='tomato',capsize=0,lw=2)
	# plot.errorbar(pwne_data_eng[0:3],pwne_data_upper_flux[0:3],yerr=pwne_data_upper_flux[0:3]-pwne_data_lower_flux[0:3],color='tomato',fmt='o',ecolor='tomato',uplims=True,lw=2)

	# plot stuff
	plot.xscale('log')
	plot.yscale('log')
	plot.ylim(yrange)
	plot.xlim(xrange)
	plot.xlabel(r'Energy (MeV)')
	plot.ylabel(r'$E^2 \times $ Flux (arbitrarily scaled)')

	im1=plot.imread('data/AGN_UnifiedModel.jpg')
	newax1=fig.add_axes([0.73,0.65,0.15,0.15],anchor='NE')
	newax1.imshow(im1)
	newax1.axis('off')

	im2=plot.imread('data/pulsar.jpg')
	newax2=fig.add_axes([0.73,0.4,0.15,0.15],anchor='NE')
	newax2.imshow(im2)
	newax2.axis('off')

	im3=plot.imread('data/Classical_Nova_Final.jpg')
	newax3=fig.add_axes([0.73,0.2,0.15,0.15],anchor='NE')
	newax3.imshow(im3)
	newax3.axis('off')

	plot.savefig('SED_science_themes.png', bbox_inches='tight')
	plot.show()
	plot.close()

	return
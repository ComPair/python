#!/usr/bin/env python
"""
------------------------------------------------------------------------

Scripts to make simulated observations

------------------------------------------------------------------------
"""

import numpy as np
import matplotlib.pylab as plot
from astropy.io import ascii
#from astropy.io import fits
from scipy import interpolate,integrate
import FigureOfMeritPlotter 
import EnergyRes
import os
import re
import importlib
import glob
import EventAnalysis

def plot_AMEGO_background_sim(dir=None,exposure=100,doplot=True,silent=False,events=None,data=None):

	if not data:
		psdir='../Simulations/AMEGO4x4PerformancePlotTraFiles/'
		data=FigureOfMeritPlotter.parseEventAnalysisLogs(psdir)

	trafiles=glob.glob(dir+'*tra')

	ntra=len(trafiles)

	energy_ComptonEvents=[]
	event_type_Compton=[]

	energy_PairEvents=[]
	event_type_Pair=[]

	nunknown=0
	ncompton=0
	npair=0
#	unknownInARM=0
	ncomptonInARM=0
	npairInARM=0
	ComptonInARM=[]
	PairInARM=[]

	for tra in trafiles:
		print tra
		events=EventAnalysis.parse(tra,sourceTheta=1)

		if events:
			energy_ComptonEvents=np.append([energy_ComptonEvents],[events['energy_ComptonEvents']])
			if 'Photon' in tra:
				event_type_Compton=np.append([event_type_Compton],np.repeat('photon',len(events['energy_ComptonEvents'])))
			else:
				event_type_Compton=np.append([event_type_Compton],np.repeat('particle',len(events['energy_ComptonEvents'])))

			energy_PairEvents=np.append([energy_PairEvents],[events['energy_pairElectron']+events['energy_pairPositron']])
			if 'Photon' in tra:
				event_type_Pair=np.append([event_type_Pair],np.repeat('photon',len(events['energy_pairElectron'])))
			else:
				event_type_Pair=np.append([event_type_Pair],np.repeat('particle',len(events['energy_pairElectron'])))

			nunknown=nunknown+events['numberOfUnknownEventTypes']
			ncompton=ncompton+events['numberOfComptonEvents']
			npair=npair+events['numberOfPairEvents']

			wCompinARM,wPairinARM=FigureOfMeritPlotter.applyARM(data,events,angleSelection=1.0)

			if len(wCompinARM)>0:
				mask=np.zeros(events['numberOfComptonEvents'],dtype=bool)
				mask[wCompinARM]=True
				ComptonInARM=np.append(ComptonInARM,mask)
				ncomptonInARM=ncomptonInARM+len(wCompinARM)
			
			if len(wPairinARM)>0:
				mask=np.zeros(events['numberOfPairEvents'],dtype=bool)
				mask[wPairinARM]=True
				PairInARM=np.append(PairInARM,mask)
				npairInARM=npairInARM+len(wPairinARM)


	if not silent:
		print 'Unknown = ',nunknown
		print 'Compton = ',ncompton
		print 'Pair = ',npair
#		print 'Unknown in ARM = ',unknownInARM
		print 'Compton in ARM = ',ncomptonInARM
		print 'Pair in ARM = ',npairInARM


	if doplot:
		fig=plot.figure()
		bins = 10**(np.arange(-2,3,0.1))
		# plot.hist(energy_ComptonEvents[event_type_Compton=='photon']/exposure,bins=bins,log=True,label='photon - Compton')
		# plot.hist(energy_ComptonEvents[event_type_Compton=='particle']/exposure,bins=bins,log=True,label='particle - Compton')#,alpha=0.5)
		# plot.hist(energy_ComptonEvents[np.where((event_type_Compton=='photon') & (ComptonInARM==1))]/exposure,bins=bins,log=True,label='photon - Compton in ARM')#,alpha=0.5)
		# plot.hist(energy_ComptonEvents[np.where((event_type_Compton=='particle') & (ComptonInARM==1))]/exposure,bins=bins,log=True,label='particle - Compton in ARM')#,alpha=0.5)

		# plot.hist(energy_PairEvents[event_type_Pair=='photon']/exposure,bins=bins,log=True,label='photon - Pair')#,alpha=0.5)
		# plot.hist(energy_PairEvents[event_type_Pair=='particle']/exposure,bins=bins,log=True,label='particle - Pair')#,alpha=0.5)
		# plot.hist(energy_PairEvents[np.where((event_type_Pair=='photon') & (PairInARM==1))]/exposure,bins=bins,log=True,label='photon - Pair in ARM')#,alpha=0.5)
		# plot.hist(energy_PairEvents[np.where((event_type_Pair=='particle') & (PairInARM==1))]/exposure,bins=bins,log=True,label='particle - Pair in ARM')#,alpha=0.5)

		plot.hist(energy_ComptonEvents/exposure,bins=bins,log=True,label='Compton')
		plot.hist(energy_ComptonEvents[np.where(ComptonInARM==1)]/exposure,bins=bins,log=True,label='Compton in ARM')#,alpha=0.5)

		plot.hist(energy_PairEvents/exposure,bins=bins,log=True,label='Pair')#,alpha=0.5)
		plot.hist(energy_PairEvents[np.where(PairInARM==1)]/exposure,bins=bins,log=True,label='Pair in ARM')#,alpha=0.5)


		plot.legend(fontsize=8)
		plot.xscale('log')
		plot.xlabel('Energy (MeV)')
		plot.xlim([0.1,1e3])
		plot.ylabel(r'Background Count Rate (ph s$^{-1}$)')
		plot.savefig('Sim_all_Background.pdf', bbox_inches='tight')
		plot.savefig('Sim_all_Background.png', bbox_inches='tight')
		plot.show()

	# bins = 10**(np.arange(-2,3,0.1))
	# hist=plot.hist(energy/exposure,bins=bins,log=True)
	# w=np.where(hist[1] != 0)
	# background_eng=hist[1][w]
	# background_rate=hist[0][w]

	### bkg/Aeff/Angres = ph/cm2/s/sr
	### * E**2

	return #background_eng,background_rate,event_type


def plot_AMEGO_background(data,showbackground=False,plotTotals=False,plothandwavy=False,\
	varyalbedo=False,varybkg=False,plotCosmicPhotons=False,plotall=False,plotCosmicRays=False,\
	plotAlbedo=False):


	Energy=[]
	for key in data.keys():
		energy = float(key.split('_')[1].replace('MeV',''))
		Energy.append(energy)
	Energy = np.array(Energy)
	Energy=np.unique(Energy)

	# digitized from Stong, Moskalenko & Reimer 2000, Figure 8, top right
  	# high latitude |b|>5 deg
  	# multiply by 10 to vaguely account for the albedo background

	#oldeng2=np.array([0.10355561,0.3534914,1.2920963,4.659387,8.969312,18.735151,38.081676,69.40132,144.98259,227.4451,342.42523,462.24567,725.01324,939.413,1908.1061,28725.793])
	#olde2int2=np.array([2.7943178E-4,3.57757E-4,4.8821748E-4,6.806025E-4,8.0072926E-4,9.1560354E-4,0.0010469892,0.0011638523,0.0013691497,0.0015439879,0.0016334692,0.0017039803,0.0018284274,0.0018672496,0.0017879958,0.0014717471])

	# Alex's new background numbers from Gruber et al. (1999) and Weidenspointer et al. (2000) and >100 MeV from Ackermann et al. (2015)
	#alex_eng2=np.array([0.5,0.8,1.0,2.0,3.0,5.0,8.0,10.0,50.0,100.,200,500])
	#alex_e2int2=np.array([2e-2,1e-2,7e-3,3e-3,2e-3,8e-4,4e-4,3e-4,3e-5,3.2e-6,2e-6,6e-7])*alex_eng2

	# From Acero et al. (2016) - arxiv:1602.07246 |b| > 10 deg Galactic Diffuse
	#lateng=np.array([59.03219,85.70306,130.86836,206.94946,279.0981,377.88,515.6217,751.6096,1202.7435,1746.5424,2439.6077,3421.0806,5185.3374,7920.732,12578.761,19214.623,39437.465,88850.945,493146.28])
  	#late2int=np.array([8.6141995E-4,0.0011453591,0.0015459998,0.002040656,0.0022842947,0.002390658,0.0024465262,0.0025040577,0.002326049,0.0020965973,0.0018339512,0.0015110897,0.0011299471,8.449544E-4,6.041775E-4,4.417733E-4,2.6613174E-4,1.7280209E-4,7.286812E-5])
  	lateng=np.array([58.665302,83.7944,127.701385,212.20918,296.02475,493.0605,740.34045,1265.6293,2019.2109,3006.2268,4828.027,8546.594,18742.852,42185.098,152450.55,496614.97])
  	late2int_galactic=np.array([8.653016E-4,0.0011343559,0.0015828605,0.0020333533,0.0022578337,0.002416496,0.0023796277,0.002305653,0.0019558307,0.0016045898,0.0011626304,7.918069E-4,4.5331568E-4,2.5003447E-4,1.3304557E-4,7.2556504E-5])
  	# From Ackermann et al. (2015) - ApJ 799 86 isotropic EGB
  	lateng_igrb=np.array([120,170,240,340,490,690,900,1350,1950,2750,3850,5450,7750,11050,15500,22200,31000,43500,61500,86000,120000,170000,245000,350000,495000,700000.])
  	lat_igrb=np.array([3.7e-6,2.3e-6,1.5e-6,9.7e-7,6.7e-7,4.9e-7,3e-7,1.8e-7,1.1e-7,6.9e-8,4.2e-8,2.6e-8,1.7e-8,1.2e-8,6.8e-9,4.4e-9,2.7e-9,1.8e-9,1.1e-9,6.2e-10,3.1e-10,1.9e-10,8.9e-11,6.3e-11,2.1e-11,9.7e-12])
  	late2int_igrb0=lat_igrb*lateng_igrb
  	tck=interpolate.splrep(np.log10(lateng_igrb),np.log10(late2int_igrb0),s=0)
	late2int_igrb=10**interpolate.splev(np.log10(lateng),tck,der=0)
	late2int=late2int_galactic+late2int_igrb
  	
  	# COMPTEL * EGRET Galactic Diffuse from Gruber et al. (1999)
  	gruber_eng=np.array([2.978623,5.1983213,9.07216,13.32116,19.94295,32.241817,44.707794,72.33151,136.47008,278.06522,545.20044,1132.5265,3079.0847,6774.5522,17384.865,41301.277,105963.19,317014.44,1315024.2,6868901.5,2.2191038E7,8.6879376E7])*1e-3
  	gruber_e2int=np.array([5.278219,4.1341214,3.2380166,2.592378,1.8563249,1.2433603,0.8325035,0.36483333,0.13988705,0.049062684,0.01799201,0.00590178,0.0014491071,4.9711246E-4,1.3645743E-4,4.6819663E-5,1.4694342E-5,4.219293E-6,8.482257E-7,1.4922263E-7,4.098341E-8,9.63159E-9])*gruber_eng

  	# COMPTEL extragalactic background from Weidenspointner (2001)
  	#wp_eng=np.array([0.10258658,0.15835446,0.34505856,0.56054914,0.9957529,1.9590781,3.4359305,10.434804,57.75009,135.85852,377.31998])
  	#wp_e2int=np.array([0.024865912,0.016798664,0.010788648,0.0067556994,0.00433873,0.002821951,0.0022472343,0.0021910856,0.002363885,0.0015768929,0.0014071865])

  	## Combining things
  	eng2=np.append(gruber_eng[0:16],lateng[2:])
  	e2int2=np.append(gruber_e2int[0:16],late2int[2:])
	# interpolate background at our energies
	tck=interpolate.splrep(np.log10(eng2),np.log10(e2int2),s=0)
	logbackground=interpolate.splev(np.log10(Energy),tck,der=0)
	background=10.**logbackground

	### Andreas files
	dir="../Simulations/BackgroundFiles/"
	files=[]
	bkgs=[]
	bkgnames=[]
	for file in os.listdir(dir):
		if file.endswith(".spectrum.dat"):
			files=np.append(files,os.path.join(dir,file))
			d=ascii.read(dir+file,data_start=1,names=['thing','energy','flux'])
			bkgs.append(d)
			file.split('.')[0]
			c=re.findall('([A-Z][a-z]*)',file)
			bkgnames.append(c[0]+' '+c[1])

	plot.figure()
	i=0
	fact=1e3 # keV to MeV
	totbkg=bkgs[0]['flux']
	totphotbkg=np.zeros(len(totbkg))
	for file in files:
		if plotall:
			plot.plot(bkgs[i]['energy']*1e-3,bkgs[i]['flux']*fact*(bkgs[i]['energy']*1e-3)**2,label=bkgnames[i])
		if plotCosmicRays:
			if ("Photons" not in file) & (("Cosmic" in file) or ("Trapped" in file)):
				plot.plot(bkgs[i]['energy']*1e-3,bkgs[i]['flux']*fact*(bkgs[i]['energy']*1e-3)**2,label=bkgnames[i])
		if plotAlbedo:
			if ("Albedo" in file) & ("Photons" not in file):
				plot.plot(bkgs[i]['energy']*1e-3,bkgs[i]['flux']*fact*(bkgs[i]['energy']*1e-3)**2,label=bkgnames[i])				
		if (i > 0) & ("Photons" not in file):
			totbkg=totbkg+bkgs[i]['flux']
		# if (i > 0) & ("Photons" in file):
		# 	totphotbkg=totphotbkg+bkgs[i]['flux']
		i=i+1

	if plotall or plothandwavy:
		plot.plot(Energy,background*10.,color='red',linestyle='--',linewidth=2,label="AMEGO Preliminary (x10)")

	eng=bkgs[0]['energy']*1e-3
	EnergykeV=np.array(eng*1e3)
	w=np.where((eng>=100) & (eng <=300e3))
	w=w[0]
	tck=interpolate.splrep(np.log10(lateng),np.log10(late2int_galactic),s=0)
	interp_late2int_galactic=10**interpolate.splev(np.log10(eng[w]),tck,der=0)
	tck=interpolate.splrep(np.log10(lateng),np.log10(late2int_igrb),s=0)
	interp_late2int_igrb=10**interpolate.splev(np.log10(eng[w]),tck,der=0)

	if plotall or plotCosmicPhotons:
		i=7
		if not plotall: 
			plot.plot(gruber_eng,gruber_e2int,color='cyan',label='HEAO, COMPTEL, EGRET - Gruber et al. (1999)')
		plot.plot(eng[w],interp_late2int_galactic,'r--',color='magenta',label='LAT Galactic Diffuse - Acero et al. (2016)')
		plot.plot(eng[w],interp_late2int_igrb,'r:',color='magenta',lw=2,label='LAT Isotropic Diffuse - Ackermann et al. (2015)')

	if varybkg:
		backeff=[1e-5]
	else:
		backeff=[1]

	if plotTotals:
		plot.plot(eng,totbkg*fact*backeff*eng**2,label='Total Particle Background',color='black')

	#### Recaculate albedo photons adding LAT
	Ajello_albedo=0.0148/(pow(EnergykeV/33.7, [1,5]) + pow(EnergykeV/33.7, [1,-1.72]))
	Mizuno_albedo=np.concatenate((\
		1010.*pow(EnergykeV[EnergykeV<20000]/1000, [1,1.34]) / 1000.0 / 10000.0,\
		7290.*pow(EnergykeV[EnergykeV>=20000]/1000, [1,2]) / 1000.0 / 10000.0))#,\
#		29000*pow(EnergykeV[EnergykeV>1000000]/1000, [1.,2.79]) / 1000.0 / 10000.0))
#	LAT_albedo=pow(EnergykeV/1000,[7.4e1,2.79])*1e-3 
	LAT_albedo=np.concatenate((\
		pow(EnergykeV[EnergykeV<2e3]/1000, [0.6,2]) /1000,\
		pow(EnergykeV[(EnergykeV>=2e3) & (EnergykeV<1e6)]/1000, [5,2.4]) / 1000.0,\
		pow(EnergykeV[(EnergykeV>=1e6) & (EnergykeV<1e8)]/1000, [74,2.79]) / 1000.0,\
		pow(EnergykeV[EnergykeV>=1e8]/1000, [2.5e5,3.5]) / 1000.0))
#	LAT_albedo=(EnergykeV/1000)**(-6.9-1.9*np.log10(EnergykeV/1e3))*1e-3
#		29000*pow(EnergykeV[EnergykeV>1000000]/1000, [1.,2.79]) / 1000.0 / 10000.0))

	Rcut_desired = 12.6
	Rcut_Mizuno = 4.5
	ScalerMizuno = 1#pow(Rcut_desired/Rcut_Mizuno, [1.,1.13])

	MizunoValue=Mizuno_albedo[EnergykeV==1000]*ScalerMizuno
	AjelloValue=Ajello_albedo[EnergykeV==1000]
	ScalerAjello = 1#MizunoValue/AjelloValue
	ScalerLAT=1.#Mizuno_albedo[EnergykeV==1e5]/LAT_albedo[EnergykeV==1e5]*ScalerMizuno
#	print ScalerLAT

	albedo_bkg=np.concatenate((Ajello_albedo[EnergykeV<1000]*ScalerAjello,\
		Mizuno_albedo[(EnergykeV>=1000) & (EnergykeV<100000)]*ScalerMizuno,\
		LAT_albedo[(EnergykeV>=100000)]*ScalerLAT))
	if plotAlbedo:
		plot.plot(eng,albedo_bkg*1e3*eng**2,label='Albedo Photons')

	photbackeff=[1,1e-1,1e-2,1e-3]
	if not varyalbedo: photbackeff=[1]
	ls=['-','--',':','-.']
	i=0
	for eff in photbackeff:
	#	totphotbkg=(bkgs[2]['flux']*eff+bkgs[7]['flux'])
		totphotbkg=(albedo_bkg+bkgs[7]['flux'])
		totphotbkg[w]=totphotbkg[w]+(interp_late2int_galactic+interp_late2int_igrb)/fact/eng[w]**2
		if eff==1: 
			streff=''
		else: streff=' '+str(eff)
		if plotTotals: plot.plot(eng,totphotbkg*fact*eng**2,label='Total Photon Background',color='blue',linestyle=ls[i])
		i=i+1

	# if varybkg:
	# 	backeff=[1e-4,1e-5,1e-6]
	# else:
	# 	backeff=[1]
	# ls=['--',':','-.']
	# i=0
	
	# for eff in backeff:
	#  	tot=totbkg*eff+totphotbkg#*peff
	# # 	if plotTotals: plot.plot(eng,tot*fact*eng**2,label='Total Photon+Particle Background ('+str(eff)+')',color='black',linestyle=ls[i])
	# # 	i=i+1

	add=''
	if plotCosmicPhotons: add='CosmicPhotons_'
	if plotCosmicRays: add=add+'CosmicRays_'
	if plotAlbedo: add=add+'Albedo_'
	if plotall: add='All_'
	if plotTotals: add=add+'Totals_'
	if varybkg: add=add+'varybkg_'

	plot.legend(loc='upper left',fontsize=6)
	plot.xscale('log')
	plot.yscale('log')
	plot.xlabel(r'Energy (MeV)')
#	plot.ylabel(r'Flux (ph cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$)')
	plot.ylabel(r'$E^2 \times$ dN/dE (MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
	plot.xlim([1e-3,1e6])
	plot.ylim([1e-5,1e3])
	plot.savefig('plots/'+add+'Background.pdf', bbox_inches='tight')
	plot.savefig('plots/'+add+'Background.png', bbox_inches='tight')
	plot.show()

# 	if showbackground:
# 		plot.figure()
# 		plot.plot(oldeng2,olde2int2,color='green',label='Strong, Moskalenko, Reimer (2000)')
# #		plot.scatter(oldeng2,olde2int2,color='red')
# #		plot.annotate('Strong, Moskalenko, Reimer (2000)',xy=(1e-2,1e-4),xycoords='data',fontsize=12,color='red')
# 		plot.plot(alex_eng2,alex_e2int2,color='blue',label='Alex (new)')
# 		# plot.scatter(alex_eng2,alex_e2int2,color='blue')
# 		# plot.annotate('Alex (new)',xy=(1e2,2e-4),xycoords='data',fontsize=12,color='blue')
# #		plot.plot(lateng,late2int,color='magenta',label='Ackermann et al. (2016) - LAT Extragalactic')
# 		# plot.scatter(lateng,late2int,color='magenta')
# 		plot.plot(lateng,late2int_galactic,'r--',color='magenta',label='Acero et al. (2016) - LAT Galactic Diffuse')
# #		plot.scatter(lateng,late2int_galactic,color='magenta')
# 		plot.plot(lateng,late2int_igrb,'r:',color='magenta',lw=2,label='Ackermann et al. (2015) - LAT Isotropic')
# #		plot.scatter(lateng,late2int_igrb,color='magenta')
# #		plot.annotate('Ackermann et al. (2016) - LAT Extragalactic',xy=(5,1e-5),xycoords='data',fontsize=12,color='magenta')
# #		plot.annotate('Acero et al. (2016) - LAT Galactic',xy=(1e2,3e-3),xycoords='data',fontsize=12,color='magenta')
# 		plot.plot(gruber_eng,gruber_e2int,color='cyan',label='Gruber et al. (1999) - HEAO, COMPTEL, EGRET')
# #		plot.scatter(gruber_eng,gruber_e2int,color='cyan')
# #		plot.annotate('Gruber et al. (1999) - HEAO, COMPTEL, EGRET',xy=(2e-3,5e-2),xycoords='data',fontsize=12,color='cyan')
# 		plot.plot(wp_eng,wp_e2int,'r:',color='orange',label='Weidenspointner (2001)')
# #		plot.scatter(wp_eng,wp_e2int,color='orange')
# #		plot.annotate('Weidenspointner (2001)',xy=(2e-3,1e-3),xycoords='data',fontsize=12,color='orange')
# #		plot.plot(eng2,e2int2,color='purple')
# #		plot.scatter(eng2,e2int2,color='purple')	
# 		plot.plot(Energy,background,color='red',label='Combined Background')
# 		plot.plot(Energy,background*10,color='red',linestyle='--',label='Combined Background x10')
# #		plot.scatter(Energy,background,color='green')
# #		plot.annotate('Interpolated Used Bkg',xy=(1,1e-2),xycoords='data',fontsize=12,color='green')
# 		plot.xscale('log')
# 		plot.yscale('log')
# 		plot.xlabel(r'Energy (MeV)')
# 		plot.ylabel(r'$E^2 \times$ Flux (MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
# 		plot.xlim([1e-2,1e6])
# 		plot.ylim([1e-5,10])
# 		plot.legend(loc='upper right',fontsize=8)
# #		plot.title('Diffuse Background')
# 		plot.savefig('plots/Background.pdf', bbox_inches='tight')
# 		plot.savefig('plots/Background.png', bbox_inches='tight')
# 		plot.show()

	return eng,totbkg*fact,LAT_albedo

def AMEGO_background(particle_rejection=None,albedo_rejection=None):

	### Andreas files
	dir="../Simulations/BackgroundFiles/"
	files=[]
	bkgs=[]
	bkgnames=[]
	for file in os.listdir(dir):
		if file.endswith(".spectrum.dat"):
			files=np.append(files,os.path.join(dir,file))
			d=ascii.read(dir+file,data_start=1,names=['thing','energy','flux'])
			bkgs.append(d)
			file.split('.')[0]
			c=re.findall('([A-Z][a-z]*)',file)
			bkgnames.append(c[0]+' '+c[1])
	i=0
	keV2MeV=1e3 # keV to MeV
	totbkg=bkgs[0]['flux']
	totphotbkg=np.zeros(len(totbkg))
	for file in files:
		if (i > 0) & ("Photons" not in file):
			totbkg=totbkg+bkgs[i]['flux']
		i=i+1

	eng=bkgs[0]['energy']*1e-3
	w=np.where((eng>=100) & (eng <=300e3))
	w=w[0]

	## adding LAT background
	# From Acero et al. (2016) - arxiv:1602.07246 |b| > 10 deg Galactic Diffuse
  	lateng=np.array([58.665302,83.7944,127.701385,212.20918,296.02475,493.0605,740.34045,1265.6293,2019.2109,3006.2268,4828.027,8546.594,18742.852,42185.098,152450.55,496614.97])
  	late2int_galactic=np.array([8.653016E-4,0.0011343559,0.0015828605,0.0020333533,0.0022578337,0.002416496,0.0023796277,0.002305653,0.0019558307,0.0016045898,0.0011626304,7.918069E-4,4.5331568E-4,2.5003447E-4,1.3304557E-4,7.2556504E-5])
  	# From Ackermann et al. (2015) - ApJ 799 86 isotropic EGB
  	lateng_igrb=np.array([120,170,240,340,490,690,900,1350,1950,2750,3850,5450,7750,11050,15500,22200,31000,43500,61500,86000,120000,170000,245000,350000,495000,700000.])
  	lat_igrb=np.array([3.7e-6,2.3e-6,1.5e-6,9.7e-7,6.7e-7,4.9e-7,3e-7,1.8e-7,1.1e-7,6.9e-8,4.2e-8,2.6e-8,1.7e-8,1.2e-8,6.8e-9,4.4e-9,2.7e-9,1.8e-9,1.1e-9,6.2e-10,3.1e-10,1.9e-10,8.9e-11,6.3e-11,2.1e-11,9.7e-12])
  	late2int_igrb0=lat_igrb*lateng_igrb
  	tck=interpolate.splrep(np.log10(lateng_igrb),np.log10(late2int_igrb0),s=0)
	late2int_igrb=10**interpolate.splev(np.log10(lateng),tck,der=0)
	late2int=late2int_galactic+late2int_igrb

	tck=interpolate.splrep(np.log10(lateng),np.log10(late2int_galactic),s=0)
	interp_late2int_galactic=10**interpolate.splev(np.log10(eng[w]),tck,der=0)
	tck=interpolate.splrep(np.log10(lateng),np.log10(late2int_igrb),s=0)
	interp_late2int_igrb=10**interpolate.splev(np.log10(eng[w]),tck,der=0)

	## summing photon background, scaling Albedo by 0.1
	albedo_rejection=[1]
	totphotbkg=(bkgs[2]['flux']*albedo_rejection+bkgs[7]['flux'])
	totphotbkg[w]=totphotbkg[w]+(interp_late2int_galactic+interp_late2int_igrb)/keV2MeV/eng[w]**2

	particle_rejection=[1e-4]
	tot=totbkg*particle_rejection+totphotbkg


	return Energy,Background

def pow(x,p,xnorm=1):

	norm=p[0]
	pow1=p[1]

	f=norm*(x/xnorm)**(-pow1)

	return f

def bknpow(x,p):

	norm=p[0]
	pow1=p[1]
	break1=p[2]
	pow2=p[3]

	f=np.zeros(len(x))
	f[x<break1]=norm*x[x<break1]**(-pow1)
	f[x>=break1]=norm*break1**(pow2-pow1)*x[x>=break1]**(-pow2)

	return f

def get_norm(model,flux,emin,emax,params):

	## flux = integral over model
	## get out normalization for model spectral in AMEGO_obs

	int1=integrate.quad(model,emin,emax,args=[1,params])[0]

	return norm

def AMEGO_obs(data=None,model=None,params=None,energy=None,flux=None,angle=0.8):

	# What Would AMEGO see?
	# simdir = '../Simulations/AMEGO4x4PerformancePlotTraFiles/'
	# data=FigureOfMeritPlotter.parseEventAnalysisLogs(simdir)
	# need to input source spectrum with some normalization and shape
	# 	convolve with energy resolution & effective area 
	# 	see what is above background
	#   eventually expand to light curve given spectrum and flux in some bins

	# read energy resolution for a given angle from a sim
	# combine tracked & Pair, split at 10 MeV
	# smear with EnergyRes

	models=['pow','bknpow','band','logparabola']

	emin=np.log10(0.3)	# 300 keV
	emax=4 				# 10 GeV
	if (model != None) & (energy == None):
		energy=np.logspace(emin,emax,num=50)
		mod=importlib.import_module('ScienceSims')
		flux=getattr(mod,model)(energy,params)

	ER_Energy,st,sp=FigureOfMeritPlotter.plotEnergyResolution(data,angleSelections=angle)
	EngRes=np.append(st,sp)
	smeared_energy=np.array(EnergyRes.EnergyResolutionFromFile(input_energy=energy, eres_Energy=ER_Energy, eres_Res=EngRes))

	Aeff_Untracked, Aeff_Tracked, Aeff_Pair, \
		Aeff_Energy=FigureOfMeritPlotter.plotEffectiveArea(data,angleSelections=angle)
	Aeff=np.append(Aeff_Tracked[Aeff_Energy<10],Aeff_Pair[Aeff_Energy>=10])

	f=interpolate.interp1d(np.log10(Aeff_Energy),np.log10(Aeff),bounds_error=False,fill_value="extrapolate",kind='linear')
	Aeff2=10**f(np.log10(smeared_energy))

	erg2mev=624151.  # 1 erg 
	amego_ctr=flux*Aeff2*erg2mev/smeared_energy

	bkg_eng,bkg_flux=AMEGO_background(data,varybkg=True,plottots=True,varyalbedo=True)
	Aeff3=10**f(np.log10(bkg_eng))
	bkg_ctr=bkg_flux*Aeff3/bkg_eng

	fig=plot.figure()
	m=np.invert(np.isnan(amego_ctr))
	plot.plot(smeared_energy[m],amego_ctr[m])#*smeared_energy[m]**2)
	plot.plot(bkg_eng,bkg_ctr,color='red')
	plot.xscale('log')
	plot.yscale('log')
	plot.xlabel('Energy (MeV)')
	plot.ylabel(r'dN/dE (ph cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
	plot.show()

	return smeared_energy,amego_ctr

def setup_source_files(dir=None,geometryfile=None,outdir=None,exposure=1):

	if outdir == None: 
		outdir=dir
	if geometryfile == None:
		geometryfile="/data/slag2/ComPair/src/Geometry/AMEGO_4x4TowerModel/AmegoBase.geo.setup"

	for file in os.listdir(dir):
		if file.endswith(".partial.source"):
			print(file)
			f=open(dir+file,'r')
			d=f.read()
			f.close()
			comp=file.split('.')
			comp=comp[0]
			outf=open(outdir+comp+'.source','w')

			s=['# Compute Background Run for Cosima',\
			'# A background source at the zenith',\
			'Version          1',\
			'Geometry         '+geometryfile,\
			'CheckForOverlaps 1000 0.01',\
			'# Physics list',\
			'PhysicsListEM                        Livermore',\
			'PhysicsListEMActivateFluorescence    false',\
			'# Output formats','StoreCalibrated                      true',\
			'StoreSimulationInfo                  true',\
			'StoreSimulationInfoIonization        false',\
			'DiscretizeHits                       true',
			' ',\
			'Run SpaceSim',\
			'SpaceSim.FileName         '+comp,\
			'SpaceSim.Time             '+str(exposure),d]


			[outf.write(st+'\n') for st in s]

			outf.close()








	
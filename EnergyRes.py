#!/usr/bin/env python
"""
------------------------------------------------------------------------

Script to add energy resolution:

Author: Regina Caputo (regina.caputo@nasa.gov)
Date: April 14, 2017

Usage Examples: Returns Energy smeared with a Gaussian
import EnergyRes
EnergyRes.smear(10.,0.5)
EnergyRes.ResolutionRange([5,10,25],[0.1,0.05,0.01])
EnergyRes.EnergyResolutionFromFile([5,10,25],filename="file.txt",doPlot=True)
EnergyRes.AngularResolutionFromFile([5.],[100.],filename="file.txt",doPlot=True)

------------------------------------------------------------------------
"""

#import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy import interpolate

def smear(val, sigma):
	"""
	This is a function to smear a single given value energy by a gaussian with width sigma
	"""
	s=np.random.normal(val, sigma, 1)

	return s[0]


def ResolutionRange(val, sigma):
	"""
	This is a function to smear an array of vals (energies or angles) by a gaussian with an array of widths
	"""
		
	if(len(val) != len(sigma)):
		if len(sigma) !=1:
			print "provide either one sigma for all values or sigmas for every value"
			return

	sVals=[]

	for i in range(len(val)):
		if len(sigma)==1:
			sVals.append(smear(val[i],sigma[0]))
		else:
			sVals.append(smear(val[i],sigma[i]))


	return sVals


def lineFun(x, a, b):
	#takes input angle in degrees
	return a*x+b

def AngularResolutionFromFile(input_energy=None, input_angle=None, ares_energy=None, ares_ares=None, energyanglefile=None, aresfile="AResTest.txt", doPlot=False, listResults=False): 
	"""
	This is a function that takes desired angles and energies and interpolates the angular resolution 
	from the inputs from a text file of values (energy or angle) vs. sigmas
	Val is the values that you would like the interpolated resolution for
	filename is the list of energy vs. resolutions
	there is a plot function that allows you to see how the interpolation went. 
	Note: This is a simple linear interpolation... 
	"""

	energy=[]
	angle=[]

	if input_energy==None and input_angle==None:
		energyangle=open(energyanglefile,'r')
		for line in energyangle:
			data=line.split()
			energy.append(float(data[0]))
			angle.append(float(data[1]))
	elif energyanglefile==None:
		if input_energy is not None and input_angle is not None:
			energy=input_energy
			angle=input_angle
		else:
			print "please enter energy and angle OR a file with that information"
			return


	if ares_energy != None:
		x=ares_energy
		y=ares_ares
	else:
		file = open(aresfile, 'r')
		x=[]
		y=[]
	
		for line in file:
			vals=line.split()
			if len(vals)==2:
				x.append(float(vals[0]))
				y.append(float(vals[1]))
			else:
				print "check file format"

	interpRes=[]
	tE=[]

	#print x, y 

	f = interp1d(x, y)
	#f = interp1d(x, y, kind='cubic')
	#f=interpolate.splrep(x, y, s=0)
	xend=len(x)-1
		
	interpRes=f(energy)

	if doPlot:
		#xnew=np.linspace(x[0],x[xend],num=41,endpoint=True)
		plt.plot(x,y,'o', energy, interpRes, 'x')
		plt.legend(['Measured Eres','Linear Interp', 'Interp. Eres'],loc='best')
		plt.xscale('log')
		plt.yscale('log')
		plt.show()

	smearedValues = ResolutionRange(angle, interpRes)

	for j in range(len(smearedValues)):
		tE.append(round(smearedValues[j],2))

	if listResults:
		print "angle: ", angle
		print "Smeared Values: ", tE
		print "Ares: ", interpRes


	return smearedValues


def EnergyResolutionFromFile(input_energy=None, energyfile=None, eres_Energy=None, eres_Res=None, filename="EResTest.txt", doPlot=False, listResults=False):
	"""
	This is a function that takes desired energies and interpolates the energy resolution 
	from the inputs from a text file of values (energy or angle) vs. sigmas
	Val is the values that you would like the interpolated resolution for
	filename is the list of energy vs. resolutions
	there is a plot function that allows you to see how the interpolation went. 
	Note: This is a simple linear interpolation... 
	"""

	if input_energy==None:
		f_energy=open(energyfile,'r')
		for line in f_energy:
			data=line.split()
			energy.append(float(data[0]))
	elif energyfile==None:
		if input_energy is not None:
			energy=input_energy
		else:
			print "please enter energy and angle OR a file with that information"
			return

	if eres_Energy != None:
		x=eres_Energy
		y=eres_Res
	else:	
		file = open(filename, 'r')
		x=[]
		y=[]

		for line in file:
			vals=line.split()
			if len(vals)==2:
				x.append(float(vals[0]))
				y.append(float(vals[1]))
			else:
				print "check file format"

	interpRes=[]
	tE=[]

	f = interp1d(x, y,bounds_error=False,fill_value='extrapolate')
	#f = interp1d(x, y, kind='cubic')
	#f=interpolate.splrep(x, y, s=0)
	xend=len(x)-1
		
	interpRes=f(energy)

	if doPlot:
		xnew=np.linspace(x[0],x[xend],num=41,endpoint=True)
		plt.plot(x,y,'o',xnew,f(xnew),'--', energy, interpRes, 'x')
		plt.legend(['Measured Eres','Linear Interp', 'Interp. Eres'],loc='best')
		plt.xscale('log')
		plt.yscale('log')
		plt.show()

	smearedValues = ResolutionRange(energy, interpRes)

	for j in range(len(smearedValues)):
		tE.append(round(smearedValues[j],2))

	if listResults:
		print "Energies: ", energy
		print "Smeared Values: ", tE
		print "Eres: ", interpRes


	return smearedValues











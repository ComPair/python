#!/usr/bin/env python
"""
------------------------------------------------------------------------

Script to add energy resolution:

Author: Regina Caputo (regina.caputo@nasa.gov)
Date: April 14, 2017

Usage Examples: Returns Energy smeared with a Gaussian
import EnergyRes
EnergyRes.Resolution(E)

------------------------------------------------------------------------
"""

#import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy import interpolate

def Resolution(Energy, sigma):
	"""
	This is a function to smear a single given Energy by a gaussian with width sigma
	"""
	s=np.random.normal(Energy, sigma, 1)

	return s[0]

def ResolutionRange(Energy, sigma):
	"""
	This is a function to smear an array of Energies by a gaussian with an array of widths
	"""
		
	if(len(Energy) != len(sigma)):
		if len(sigma) !=1:
			print "provide either one sigma for all E or sigmas for every E"
			return

	sEnergy=[]

	for i in range(len(Energy)):
		if len(sigma)==1:
			sEnergy.append(Resolution(Energy[i],sigma[0]))
		else:
			sEnergy.append(Resolution(Energy[i],sigma[i]))

	return sEnergy


def lineFun(x, a, b):
	#takes input angle in degrees
	return a*x+b


def ResolutionFromFile(Energy, filename="EResTest.txt", doPlot=False):
	"""
	This is a function that inputs a text file of Energy vs. sigmas
	"""

	file = open(filename, 'r')
	x=[]
	y=[]

	for line in file:
		vals=line.split()
		x.append(float(vals[0]))
		y.append(float(vals[1]))

	f = interp1d(x, y)
	#f = interp1d(x, y, kind='cubic')
	#f=interpolate.splrep(x, y, s=0)
	xend=len(x)-1

	interpRes=f(Energy)

	if doPlot:
		xnew=np.linspace(x[0],x[xend],num=41,endpoint=True)
		plt.plot(x,y,'o',xnew,f(xnew),'--', Energy, interpRes, 'x')
		plt.legend(['Measured Eres','Linear Interp', 'Interp. Eres'],loc='best')
		plt.xscale('log')
		plt.yscale('log')
		plt.show()

	return interpRes










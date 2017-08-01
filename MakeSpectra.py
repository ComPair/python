#!/usr/bin/env python
"""
------------------------------------------------------------------------

A script to so something. 

Following: https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/source_models.html

Author: 
Date: 

Usage Examples:


------------------------------------------------------------------------
"""

import os
import time
import sys
import fileinput
import numpy
from itertools import product, combinations
from collections import OrderedDict
from scipy.stats import norm
from scipy.optimize import leastsq
import scipy.optimize
import math
import glob


#These functions represent the dN/dE of the source

def powerlaw(x,N0,E0,g):

	return N0*(x/E0)**(g)

def powerlaw2(x,N0,Emin,Emax,g):

	return N0*(g+1)*x**g/(Emax**(g+1)-Emin**(g+1))

def brokenpowerlaw(x,N0,Eb,g1,g2):
	result=-99999999.
	if x<Eb:
		result=powerlaw(x,N0,Eb,g1)
	else:
		result=powerlaw(x,N0,Eb,g2)
	return result

#def brokenpowerlaw2():

def logparabola(x,N0,Eb,a,b):

	return N0*(x/Eb)**(-1.0*(a+b*numpy.log10(x/Eb)))

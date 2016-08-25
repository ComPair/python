import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.integrate import quad
import math
import ROOT

from astropy.io import ascii
from operator import truediv

E0=np.array([0.1, 0.316, 1.0, 3.16, 10., 31.6, 100., 316.])
A_eff=np.array([10., 10., 10., 10., 10., 10., 10., 10.])
PSF=np.array([1., 1., 1., 1., 1., 1., 1., 1.]) # use the one at Cos(th)=0.8
Bbkg=np.array([0.00278318, 0.00346008, 0.00447618, 0.00594937, 0.00812853, 0.0100297, 0.0124697, 0.0161290])
time = 6.3*10**6

def omega(PSF):
    return 2*math.pi*(1-np.cos(2*PSF*math.pi/180.))
    
def Isrc(E,time, Aeff, nsig, domega, bkg):
    
    arg = np.sqrt((nsig**4/4.)+(nsig**2*Bbkg*Aeff*time*domega/E))
    num = E/(Aeff*time)*(nsig**2/2+arg)
    return num

test = Isrc(E0, time, A_eff, 3., omega(PSF), Bbkg)

print test


#def plsrc(x, A, index, E0):
#    return A*(x/E0)**(-1.*index)

#def intplsrc(Emin, Emax):
#
#    nbEbins = 10000
#    de      = (math.log10(Emax)-math.log10(Emin))/nbEbins
#    logEmin = math.log10(Emin)
#    Eedges  = []
#    for iE in range(0,nbEbins+1):
#        Eedges.append(10**(logEmin+iE*de))
#
#    phFlux=0.0
#    for iE in range(0,nbEbins):
#        deltaE = Eedges[iE+1]-Eedges[iE]
#        Ecenter = Eedges[iE]+deltaE/2.
#        phFlux += plsrc(Ecenter, 1, 2, 1)*deltaE
#
#    return phFlux

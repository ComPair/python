#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Henrike Fleischhack, CUA/NASA-GSFC.                                   #
#                                                                              #
#------------------------------------------------------------------------------#

"""
Calculate integral flux sensitivities for bursts at a given time scale and energy range, given effective area and background rates. No ARM cuts!

This is similar to running MEGALib's SensitivityOptimizer with very open cuts
 (and no actual optimization), however:
    * Here, we get the background rates from (binned) histograms rather than counting events,
      so the resolution is limited by the width of the energy bins.
    * Here, we find the effective area by interpolating the effective area vs true energy
      curves from mono-energetic simulations and finding the weighted average
      given a certain energy spectrum.
      The SensitivityOptimizer counts events falling in a given *reconstructed* energy range.
    
Note that no cuts are applied here, including no ARM cuts. 
This is appropriate for very short time scales (bursts, giant magnetar flares)
but not for the survey!
    
usage: BurstSensitivities.py [-h] -f [EAFILES [EAFILES ...]] -b BGFILE
                             [-em EMIN] [-eM EMAX] [-t EXPOSURE]
                             [-s SIGNIFICANCE] [-x INDEX] [-c CUTOFF] [-o]

optional arguments:
  -h, --help            show this help message and exit
  -f [EAFILES [EAFILES ...]], --eafiles [EAFILES [EAFILES ...]]
                        input .txt file(s) with effective areas (from
                        plotFigureOfMwerit.py) (default: None)
  -b BGFILE, --bgfile BGFILE
                        input .root file(s) with background rate histograms
                        (already normalized by exposure time) (default: None)
  -em EMIN, --emin EMIN
                        Minimum energy (MeV) (default: 0.1)
  -eM EMAX, --emax EMAX
                        Maximum energy (MeV) (default: 10)
  -t EXPOSURE, --exposure EXPOSURE
                        Desired exposure time (s) (default: 1)
  -s SIGNIFICANCE, --significance SIGNIFICANCE
                        Desired signal/noise threshold (sigma) (default: 6)
  -x INDEX, --index INDEX
                        Spectral index (default: -2)
  -c CUTOFF, --cutoff CUTOFF
                        Cutoff energy (MeV) (default: inf)

"""

import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from scipy.interpolate import interp1d
from scipy.optimize import minimize

__description__ = 'Integral flux sensitivities for bursts'


"""Command-line switches.
"""
formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-f', '--eafiles', type=str, required=True, nargs='*',
                    help='input .txt file(s) with effective areas (from plotFigureOfMwerit.py)')
PARSER.add_argument('-b', '--bgfile', type=str, required=True,
                    help='input .root file(s) with background rate histograms (already normalized by exposure time)')
PARSER.add_argument('-em', '--emin', type=float, default=0.1,
                    help='Minimum energy (MeV)')
PARSER.add_argument('-eM', '--emax', type=float, default=10,
                    help='Maximum energy (MeV)')
PARSER.add_argument('-t', '--exposure', type=float, default=1, help='Desired exposure time (s)')
PARSER.add_argument('-s', '--significance', type=float, default=6, help='Desired signal/noise threshold (sigma)')
PARSER.add_argument('-x', '--index', type=float, default=-2, help='Spectral index')
PARSER.add_argument('-c', '--cutoff', type=float, default=np.inf, help='Cutoff energy (MeV)')

FLAG_STR = '(?<=\_)\w+(?=\.txt)'


def get_bg_counts(root_file, emin, emax, type):
    """read in ROOT file with BG rate histograms (dN/dE/dT) and return
    the rate dN/dT for a given energy range and event type (UC, TC, P)"""
    
    bgRate = {}
    _file0 = ROOT.TFile(root_file)
    
    
    try:
    
            totalBG = _file0.Get(f"totalBG_{type}")
            bin1 = totalBG.FindBin( emin*1e3 ) #energies are in keV!
            bin2 = totalBG.FindBin( emax*1e3 )
            rate = totalBG.Integral(bin1, bin2, "width" )
    
    except:
            rate = None

    return rate


def read_effective_areas(ea_files):
    """read in txt files with effective area vs true energy (as in compareFigureOfMerit.py)
    return a dictionary with the effective areas (as interp1d objects) and event type (UC, TC, P)"""

    theEAs = {}

    for ea_file in ea_files:

        print('Parsing %s ...'%ea_file)

        file_basename = os.path.basename(ea_file)
        m = re.search(r'%s'%FLAG_STR, file_basename)
        type = m.group(0)
    
        en, ea = np.loadtxt(ea_file).T
        
        ea[np.isnan(ea)] = 0

        f = interp1d(en, ea, kind='cubic')
        
        theEAs[type] = f
        
    return theEAs


def average_effective_area( EA, emin, emax, index, cutoff):
    """return effective area averaged over a given energy range, weighted by a given energy spectrum (cutoff powerlaw)"""

    E = np.linspace( emin, emax, 1000)
    dNdE = E**args.index * np.exp(-E/args.cutoff)
    
    x1 = 0
    x2 = 0
    
    for e, n in zip(E, dNdE):
        x1 += n * EA(e)
        x2 += n
    
    if x2 == 0:
        return 0
    else:
        return x1/x2
    
    
def get_sensitivity(time, significance, Aeff, BGrate):
    """Calculate minimum detectable flux (in photons/cm2/s) for a given exposure time, detection threshold,
        effective area, and background rate"""
        
    arg = np.sqrt((significance**4/4.)+(significance**2*BGrate*time))
    num = (significance**2/2+arg)/(Aeff*time)
    return num



if __name__ == '__main__':

    args = PARSER.parse_args()

    theEAs = read_effective_areas(args.eafiles)
    
    print(f"{args.significance}σ sensitvities for T = {args.exposure:.3g}s")
    
    Emin = args.emin
    Emax = args.emax

    for type in theEAs.keys():

        Aeff = average_effective_area( theEAs[type], Emin, Emax, args.index, args.cutoff )
        BGrate = get_bg_counts(args.bgfile, Emin, Emax, type )

        #only try to calculate the sensitivity if we have events
        if BGrate is not None and (Aeff > 0):
            sensi = get_sensitivity(args.exposure, args.significance, Aeff, BGrate)
            
            print( f"{type:2s} Events: {Emin:5.2f} - {Emax:5.2f} MeV, "
                    f" Aeff = {Aeff:6.4g} cm2, "
                    f"BG = {BGrate:5.1f} Hz, min. Flux = {sensi:7.4g} ph/cm2/s" )

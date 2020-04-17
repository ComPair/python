
# Author: Donggeun Tak (takdg123@gmail.com)
# Date: April 3rd, 2020

# Before using this code, check geofile, Log_E, and angles parameters.
# This code will generate source files and runCosima.sh

# Modifications: 
# Jan (lommler@uni-mainz.de)
# Date: April 15th, 2020
# removed fixed seed (we don't want to relive lhcb...)
# sims are now gziped, BEWARE: the automated analysis pipeline needs some tweaks to works with gziped files, for event analysis either modify the parts using the sim file and extract the tra-files manually or build a wrapper -> still to do
# sim files are now stored in a subfolder 'data/', you need to keep the concatinating sim file (it is the one without '.incXX.' in the name) in the parent directory of data/ (mcosima assumes the same path as for simulation when creating the concatination file). 
# changed PhysicsList to LivermorePol (you never know when some one asks you to look into polarization in the least expected moment...
# Globals:
# int threads:  lets you define the number of threads you want to allocate to simulating, 
# int targetTriggers: number of total triggered events for each energy bin
# int triggers:  gives the number of triggered events that are simulated derived from threads and targetTriggers

from numpy import *
from math import *


#here put your geometry file and source type
geofile= 'Geometry/AMEGO_Midex/TradeStudies/Calorimeter/BaseSideLogs/AmegoBase.geo.setup'
OneBeam= 'FarFieldPointSource'
threads = 20
targetTriggers = 100000
triggers = int(math.ceil(targetTriggers/threads))

#define your energies and angles 
#Log_E=[2.2,2.5,2.7,3,3.2,3.5,3.7,4,4.2,4.5,4.7,5,5.2,5.5,5.7,6,6.2,6.5,6.7]
#angles  =[0,25.8,36.9,45.6,53.1,60]
Log_E=[2.2,2.5,2.7,3,3.2,3.5,3.7,4,4.2,4.5,4.7,5,5.2,5.5,5.7,6]
angles=[0]
#gives the cosTheta array
def ang2cos(allAng):
   ang =[]
   for i in allAng:
      a= round(cos(math.radians(i)),1)
      ang.append(a) 
   return ang

#in keV [316,501,1000,1585, ... ]
def logE2ene(allEne):
   ene =[]
   for ee in allEne:
      a=int(10**ee)
      ene.append(a) 
   return ene

energies=logE2ene(Log_E)
cos_ang =ang2cos(angles)

with open("./runMcosima.sh", mode='w') as f:
   for myene in energies:
      for cosTh,ang in zip(cos_ang,angles):
         
         # this is to print all the parameters combinations
         #print (geofile, OneBeam, myene/1000., cosTh, OneBeam, ang, myene)
         
         #this is just a long string, with all the raws of the .source file, and the energies/angles values
         string= "# An example run for Cosima \n# This was created with the python wrapper --> create_source_file.py <--\n\nVersion          1 \nGeometry         %s // Update this to your path \nCheckForOverlaps 1000 0.01 \nPhysicsListEM    LivermorePol \n\nStoreCalibrate                 true\nStoreSimulationInfo            true\nStoreOnlyEventsWithEnergyLoss  true  // Only relevant if no trigger criteria is given! \nDiscretizeHits                 true \n\nRun FFPS \nFFPS.FileName              data/%s_%.3fMeV_Cos%.1f \nFFPS.NTriggers             %s \n\n\nFFPS.Source One \nOne.ParticleType        1 \nOne.Beam                %s  %.1f 0 \nOne.Spectrum            Mono  %i\nOne.Flux                1000.0 "%(geofile, OneBeam, myene/1000., cosTh, triggers, OneBeam, ang, myene)
         source_file='%s_%.3fMeV_Cos%.1f.source'%(OneBeam,myene/1000.,cosTh)
         sf=open(source_file,'w')
         sf.write(string)
         sf.close()

         runCode='%s_%.3fMeV_Cos%.1f.source'%(OneBeam,myene/1000.,cosTh)
         # -t gives the number of simulations to start, 
         # -m gives the number of simultanous threads to start (leave it at least as high as the value of t, otherwise the console gets spammed with waiting messages),
         # -w blocks the console till the simulations are finished (advantage: if you have to kill the simulations you can ctrl+c-kill all the sims and don't have to use the killAllCosimas script or hunt the threads down manually)
         # -z gzips the sim file, beware: you have to modify the automated analysis scripts to make everything work again
         f.write("mcosima -t {} -m {} -z -w {}\n".format(threads, threads, runCode))
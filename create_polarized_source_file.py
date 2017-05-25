"""
------------------------------------------------------------------------
This script aims being part of the ComPair multicore analysis chain

It is a wrapper through energies and angles to create cosima source files
One only needs to define the path to the geometry file and the type of source (e.g. FarFieldPointSource)
Energies and angles may be adjusted according to the use preference

Author: Sara Buson (sara.buson@gmail.com)
Author: Fabian Kislat (fkislat@wustl.edu)
------------------------------------------------------------------------
"""

from numpy import *
from math import *

#here put your geometry file and source type
geofile= '$COMPAIRPATH/Geometry/AMEGO_4x4TowerModel/AmegoBase.geo.setup'
OneBeam= 'FarFieldPointSource'
Output=  'PolarizedFarFieldPointSource'

#define your energies and angles
#unlike the unpolarized simulations, focus on the Compton regime
Log_E=[2.2,2.5,2.7,3,3.2,3.5,3.7,4,4.2]
angles  =[0,25.8,36.9,45.6,53.1,60]



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

for myene in energies:
   for cosTh,ang in zip(cos_ang,angles):
      
      # this is to print all the parameters combinations
      #print (geofile, OneBeam, myene/1000., cosTh, OneBeam, ang, myene)
      
      #this is just a long string, with all the raws of the .source file, and the energies/angles values
      string="""# An example run for Cosima
# This was created with the python wrapper --> create_source_file.py <--

Version          1
Geometry         %s // Update this to your path
CheckForOverlaps 1000 0.01
PhysicsListEM    Livermore-Pol

StoreCalibrate                 true
StoreSimulationInfo            true
StoreOnlyEventsWithEnergyLoss  true  // Only relevant if no trigger criteria is given!
DiscretizeHits                 true

Run FFPS
FFPS.FileName              %s_%.3fMeV_Cos%.1f
FFPS.NTriggers             100000


FFPS.Source One
One.ParticleType        1
One.Beam                %s  %.1f 0
One.Spectrum            Mono  %i
One.Flux                1000.0
One.Polarization        RelativeX 1.0 0.
""" % (geofile, Output, myene/1000., cosTh, OneBeam, ang, myene)
      
      source_file='%s_%.3fMeV_Cos%.1f.source'%(Output,myene/1000.,cosTh)
      sf=open(source_file,'w')
      sf.write(string)
      sf.close()

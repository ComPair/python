
from numpy import *
from math import *

#here put your geometry file and source type
geofile= '/net/slag2/dtak/Geometry/AMEGO_Midex/AmegoBase.geo.setup'
OneBeam= 'FarFieldPointSource'

#define your energies and angles 
Log_E=[2.2,2.5,2.7,3,3.2,3.5,3.7,4,4.2,4.5,4.7,5,5.2,5.5,5.7,6,6.2,6.5,6.7]
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
      string= "# An example run for Cosima \n# This was created with the python wrapper --> create_source_file.py <--\n\nVersion          1 \nGeometry         %s // Update this to your path \nCheckForOverlaps 1000 0.01 \nPhysicsListEM    Livermore \n\nStoreCalibrate                 true\nStoreSimulationInfo            true\nStoreOnlyEventsWithEnergyLoss  true  // Only relevant if no trigger criteria is given! \nDiscretizeHits                 true \n\nRun FFPS \nFFPS.FileName              %s_%.3fMeV_Cos%.1f \nFFPS.NTriggers             100000 \n\n\nFFPS.Source One \nOne.ParticleType        1 \nOne.Beam                %s  %.1f 0 \nOne.Spectrum            Mono  %i\nOne.Flux                1000.0 "%(geofile, OneBeam, myene/1000., cosTh, OneBeam, ang, myene)
      source_file='%s_%.3fMeV_Cos%.1f.source'%(OneBeam,myene/1000.,cosTh)
      sf=open(source_file,'w')
      sf.write(string)
      sf.close()
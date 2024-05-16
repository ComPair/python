
# Author: Donggeun Tak (takdg123@gmail.com)
# Date: April 3rd, 2020

# Before using this code, check geofile, Log_E, and angles parameters.
# This code will generate source files and runCosima.sh


from numpy import *
from math import *
import math

# here put your geometry file and source type
# geofile = '/Users/njmille2/ComPair/Geometry/AMEGO_Midex/AmegoXBase.geo.setup'
# geofile = '/Users/njmille2/ComPair/Geometry/AMEGO_Midex/TradeStudies/Tracker/BasePixelTracker_05mm_noSubThr/AmegoBase.geo.setup'
# geofile = '/Users/njmille2/ComPair/Geometry/AMEGO_Midex/TradeStudies/Tracker/BasePixelTracker_05mm/AmegoBase.geo.setup'
geofile = '/Users/njmille2/ComPair/Geometry/ComPair_23/ComPair23.geo.setup'
OneBeam = 'FarFieldPointSource'

# define your energies and angles
# Log_E = [2.2, 2.5, 2.7, 3, 3.2, 3.5, 3.7, 4, 4.2, 4.5, 4.7, 5, 5.2, 5.5, 5.7, 6, 6.2, 6.5, 6.7]
# angles= [0, 25.8, 36.9, 45.6, 53.1, 60]
# Log_E = [2.2, 2.7, 4.2, 5, 6.2]
# Log_E = [2.5, 3, 3.2, 3.5, 3.7, 5]
# Log_E = [2.2, 2.5, 2.7, 3, 3.2, 3.5, 3.7, 4.2]
# Log_E = [1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.7]
Log_E = [1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.7, 3, 3.2, 3.5, 3.7, 4.2, 4.4]
#Log_E = [1.4, 1.9, 3.5, 4.4]  # 1.4 = 25 keV, 4.4 = 25 MeV
#Log_E = [3.5, 4.4]

#angles = [0, 36.9]
angles = [0]


# gives the cosTheta array
def ang2cos(allAng):
    ang = []
    for i in allAng:
        a = round(cos(math.radians(i)), 1)
        ang.append(a)
    return ang


# in keV [316,501,1000,1585, ... ]
def logE2ene(allEne):
    ene = []
    for ee in allEne:
        a = int(10**ee)
        ene.append(a)
    return ene


energies = logE2ene(Log_E)
cos_ang = ang2cos(angles)

with open("./runCosima.sh", mode='w') as f:
   for myene in energies:
      for cosTh, ang in zip(cos_ang, angles):
         
         # this is to print all the parameters combinations
         # print (geofile, OneBeam, myene/1000., cosTh, OneBeam, ang, myene)
         
         # this is just a long string, with all the raws of the .source file, and the energies/angles values
         string = "# An example run for Cosima \n# This was created with the python wrapper --> create_source_file.py <--\n\nVersion          1 \nGeometry         %s // Update this to your path \nCheckForOverlaps 1000 0.01 \nPhysicsListEM    Livermore \n\nStoreCalibrate                 true\nStoreSimulationInfo            true\nStoreOnlyEventsWithEnergyLoss  true  // Only relevant if no trigger criteria is given! \nDiscretizeHits                 true \n\nRun FFPS \nFFPS.FileName              %s_%.3fMeV_Cos%.1f \nFFPS.NTriggers             400000 \n\n\nFFPS.Source One \nOne.ParticleType        1 \nOne.Beam                %s  %.1f 0 \nOne.Spectrum            Mono  %i\nOne.Flux                1000.0 "%(geofile, OneBeam, myene/1000., cosTh, OneBeam, ang, myene)
         source_file = '%s_%.3fMeV_Cos%.1f.source'%(OneBeam,myene/1000.,cosTh)
         sf = open(source_file, 'w')
         sf.write(string)
         sf.close()

         runCode = '%s_%.3fMeV_Cos%.1f.source'%(OneBeam,myene/1000.,cosTh)
         f.write("cosima -s 120 {}\n".format(runCode))

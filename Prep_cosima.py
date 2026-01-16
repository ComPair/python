
# Author: Donggeun Tak (takdg123@gmail.com)
# Date: April 3rd, 2020

# Before using this code, check geofile, Log_E, and angles parameters.
# This code will generate source files and runCosima.sh


from numpy import *
from math import *
import math
import os

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

def main(Log_E, angles, N_Triggers, geofile, OneBeam, working_directory):
    energies = logE2ene(Log_E)
    cos_ang = ang2cos(angles)

    with open(f"{working_directory}/runCosima.sh", mode='w') as f:
        for myene in energies:
            for cosTh, ang in zip(cos_ang, angles):
                
                # this is to print all the parameters combinations
                # print (geofile, OneBeam, myene/1000., cosTh, OneBeam, ang, myene)
                
                # this is just a long string, with all the raws of the .source file, and the energies/angles values
                string = f"# An example run for Cosima \n# This was created with the python wrapper --> create_source_file.py <--\n\nVersion          1 \nGeometry         %s // Update this to your path \nCheckForOverlaps 1000 0.01 \nPhysicsListEM    Livermore \n\nStoreCalibrate                 true\nStoreSimulationInfo            true\nStoreOnlyEventsWithEnergyLoss  true  // Only relevant if no trigger criteria is given! \nDiscretizeHits                 true \n\nRun FFPS \nFFPS.FileName              %s_%.3fMeV_Cos%.1f \nFFPS.NTriggers             {N_Triggers} \n\n\nFFPS.Source One \nOne.ParticleType        1 \nOne.Beam                %s  %.1f 0 \nOne.Spectrum            Mono  %i\nOne.Flux                1000.0 "%(geofile, OneBeam, myene/1000., cosTh, OneBeam, ang, myene)
                source_file = '%s_%.3fMeV_Cos%.1f.source'%(OneBeam,myene/1000.,cosTh)
                sf = open(f'{working_directory}/{source_file}', 'w')
                sf.write(string)
                sf.close()

                runCode = '%s_%.3fMeV_Cos%.1f.source'%(OneBeam,myene/1000.,cosTh)
                f.write("cosima -s 120 {}\n".format(f'{working_directory}/{runCode}'))


if __name__=='__main__':
    # here put your geometry file and source type
    geofile = '/Users/gsommer1/Software/ComPair2/Geometry/ComPair_23/ComPair23.geo.setup'
    OneBeam = 'FarFieldPointSource'

    # define your energies and angles

    Log_E = [2.0, 2.59897, 3.]

    angles = [0,30]

    N_Triggers=50000

    working_directory=os.getcwd()

    main(Log_E, angles, N_Triggers, geofile, OneBeam, working_directory)
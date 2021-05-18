# Adapted from https://github.com/ComPair/python/blob/master/Prep_cosima.py

# This code will generate source files and runCosima.sh


import numpy as np
import os

import argparse

Log_E=[1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.2,3.5,3.7,4,4.2,4.5,4.7,5,5.2,5.5,5.7,6]
#angles = [0, 25.8, 36.9, 45.6, 53.1, 60] #cos
angles = [0, 36.9] #cos

width = 1.0 # degrees

Log_E=np.array(Log_E)
angles=np.array(angles)

OneBeam= 'FarFieldAreaSource'


#gives the cosTheta array
def ang2cos(allAng):
    return np.round(np.cos(np.deg2rad(allAng)),1)

def logE2ene(allEne):
    return (10**allEne).astype(int)

energies=logE2ene(Log_E)
cos_ang =ang2cos(angles)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g",  help="Geometry file")
    parser.add_argument("-o", default=os.getcwd(), help="Output file dir (Default: current directory)")
    args = parser.parse_args()

    output_dir= args.o
    geometry_file = args.g
    
    os.makedirs(output_dir,exist_ok=True)


    with open("{0}/runCosima.sh".format(output_dir), mode='w') as f:
        for myene in energies:
            for cosTh,ang in zip(cos_ang,angles):
         
                # this is to print all the parameters combinations
                #print (geofile, OneBeam, myene/1000., cosTh, OneBeam, ang, myene)
         
                minTheta = max( ang - width, 0 )
                maxTheta = min( ang + width, 180 )
                minPhi = 0.0
                maxPhi = 360.0
            
                #this is just a long string, with all the raws of the .source file, and the energies/angles values
                string= "# An example run for Cosima \n# This was created with the python wrapper --> create_source_file.py <--\n\nVersion          1 \nGeometry         %s // Update this to your path \nCheckForOverlaps 1000 0.01 \nPhysicsListEM    LivermorePol \n\nStoreCalibrate                 true\nStoreSimulationInfo            true\nStoreOnlyEventsWithEnergyLoss  true  // Only relevant if no trigger criteria is given! \nDiscretizeHits                 true \nStoreSimulationInfoIonization  false \n\nRun FFPS \nFFPS.FileName              %s/%s_%.3fMeV_Cos%.1f \nFFPS.NTriggers             100000 \n\n\nFFPS.Source One \nOne.ParticleType        1 \nOne.Beam                %s  %.1f %.1f %.1f %.1f \nOne.Spectrum            Mono  %i\nOne.Flux                10000.0 "%(geometry_file, output_dir, OneBeam, myene/1000., cosTh, OneBeam, minTheta, maxTheta, minPhi, maxPhi, myene)
                source_file='%s/%s_%.3fMeV_Cos%.1f.source'%(output_dir, OneBeam,myene/1000.,cosTh)
                sf=open(source_file,'w')
                sf.write(string)
                sf.close()

                f.write("cosima -z -v 0 -s 120 {}\n".format(source_file))
        
    print("Run this script: {0}/runCosima.sh".format(output_dir) )



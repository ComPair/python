import numpy as np
import math

Log_E=[2.2,2.5,2.7,3,3.2,3.5,3.7,4,4.2,4.5,4.7,5,5.2,5.5,5.7,6,6.2,6.5,6.7]
#angles  =[0,25.8,36.9,45.6,53.1,60]
angles = [0, 36.9]
OneBeam= 'FarFieldPointSource'

#gives the cosTheta array
def ang2cos(allAng):
    ang =[]
    for i in allAng:
        a= round(np.cos(math.radians(i)),1)
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
with open("./runCosima.sh", mode='w') as f:
    for myene in energies:
        for cosTh,ang in zip(cos_ang,angles):
            source_file='%s_%.3fMeV_Cos%.1f.source'%(OneBeam,myene/1000.,cosTh)
            f.write("cosima -s 120 {}\n".format(source_file))
        
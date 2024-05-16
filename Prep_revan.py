
# Author: Donggeun Tak (takdg123@gmail.com)
# Date: April 3rd, 2020

# Before using this code, check geofile, Log_E, and angles parameters.
# This code will generate runRevan.sh file.

import numpy as np
import math
import argparse
import os

# Log_E=[2.2,2.5,2.7,3,3.2,3.5,3.7,4,4.2,4.5,4.7,5,5.2,5.5,5.7,6,6.2,6.5,6.7]
# Log_E = [2.2, 2.5, 2.7, 3, 3.2, 3.5, 3.7, 4.2] # 5
#Log_E = [1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.7]
Log_E = [1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.7, 3, 3.2, 3.5, 3.7, 4.2, 4.4]
# angles  =[0,25.8,36.9,45.6,53.1,60]
#Log_E = [1.4, 1.9, 3.5, 4.4]

angles = [0]
OneBeam = 'FarFieldPointSource'
geometry_file = 'ComPair23.geo.setup'


# gives the cosTheta array
def ang2cos(allAng):
    ang = []
    for i in allAng:
        a = round(np.cos(math.radians(i)), 1)
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


def create_file(source_file):
    with open('./{}.revan.cfg'.format(source_file), mode='w') as cfg:
        with open('./{}'.format(base_file)) as base:
            for line in base.readlines():
                if line.find('<GeometryFileName>') != -1:
                    cfg.write('<GeometryFileName>{}/{}</GeometryFileName>\n'.format(geometry_dir,geometry_file))
                elif line.find('<DataFileName>') != -1:
                    cfg.write('<DataFileName>{}/{}.inc1.id1.sim</DataFileName>\n'.format(source_dir, source_file))
                else:
                    cfg.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("g",  help="Geometry file dir")
    parser.add_argument("-f", default=os.getcwd(), help="Source file dir (Default: current directory)")
    parser.add_argument("-b", default='revan_AMEGO_X.cfg', help="Base file")
    args = parser.parse_args()

    source_dir = args.f
    geometry_dir = args.g
    base_file = args.b

    args = parser.parse_args()
    with open("./runRevan.sh", mode='w') as f:
        for myene in energies:
            for cosTh, ang in zip(cos_ang, angles):
                source_file = '%s_%.3fMeV_Cos%.1f'%(OneBeam,myene/1000.,cosTh)
                create_file(source_file)
                f.write("revan -a -n -f {}.inc1.id1.sim -g {}/{} -c {}.revan.cfg\n".format(source_file, geometry_dir,geometry_file, source_file))

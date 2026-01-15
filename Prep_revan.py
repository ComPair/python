
# Author: Donggeun Tak (takdg123@gmail.com)
# Date: April 3rd, 2020

# Before using this code, check geofile, Log_E, and angles parameters.
# This code will generate runRevan.sh file.

import numpy as np
import math
import argparse
import os


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



def create_file(source_dir, source_file, base_file, geofile):
    with open(f'{source_dir}/{source_file}.revan.cfg', mode='w') as cfg:
        with open(base_file) as base:
            for line in base.readlines():
                if line.find('<GeometryFileName>') != -1:
                    cfg.write('<GeometryFileName>{}</GeometryFileName>\n'.format(geofile))
                elif line.find('<DataFileName>') != -1:
                    cfg.write('<DataFileName>{}/{}.inc1.id1.sim</DataFileName>\n'.format(source_dir, source_file))
                else:
                    cfg.write(line)

def main(Log_E, angles, geofile, OneBeam, revan_config_file, source_dir):
    energies = logE2ene(Log_E)
    cos_ang = ang2cos(angles)
    with open(f"{source_dir}/runRevan.sh", mode='w') as f:
            for myene in energies:
                for cosTh, ang in zip(cos_ang, angles):
                    source_file = '%s_%.3fMeV_Cos%.1f'%(OneBeam,myene/1000.,cosTh)
                    create_file(source_dir, source_file, revan_config_file, geofile)
                    f.write("revan -a -n -f {}.inc1.id1.sim -g {} -c {}.revan.cfg\n".format(f'{source_dir}/{source_file}', geofile, f'{source_dir}/{source_file}'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", default=os.getcwd(), help="Source file dir (Default: current directory)")
    parser.add_argument("-b", default='revan_AMEGO_X.cfg', help="Base file")
    args = parser.parse_args()

    source_dir = args.f
    base_file = args.b

    args = parser.parse_args()

    Log_E = [2.0, 2.59897, 3.]

    angles = [0,30]
    OneBeam = 'FarFieldPointSource'
    geofile = '/Users/gsommer1/Software/ComPair2/Geometry/ComPair_23/ComPair23.geo.setup'

    main(Log_E=Log_E, angles=angles, geofile=geofile, OneBeam=OneBeam, revan_config_file=base_file, source_dir=source_dir)

    

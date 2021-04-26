# Adapted from https://github.com/ComPair/python/blob/master/Prep_revan.py
# This code will generate soft links of the source files
# into the output directory and generate runRevan.sh file.
# Please provide a geometry file & the full paths to the source and output directories.
# (The symlink is needed to be able to run revan on the same input files for multiple
# revan config files/geometry files without overwriting the results, since revan
# always puts the output files in the same directory as the input, and there does not
# seem to be a way to specify the file name.

import numpy as np
import math
import argparse
import os


#to be adjusted before running
Log_E=[1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.2,3.5,3.7,4,4.2,4.5,4.7,5,5.2,5.5,5.7,6]
theta_angles  =[0, 25.8, 36.9, 45.6, 53.1, 60]
phi_angles=[0.0]


OneBeam= 'FarFieldPointSource'

Log_E=np.array(Log_E)
angles=np.array(angles)
phi_angles=np.array(phi_angles)


#gives the cosTheta array
def ang2cos(allAng):
    return np.round(np.cos(np.deg2rad(allAng)),1)

#in keV [316,501,1000,1585, ... ]
def logE2ene(allEne):
    return (10**allEne).astype(int)


energies=logE2ene(Log_E)
cos_ang =ang2cos(angles)

def create_config_file(source_dir, source_file, output_dir, geometry_file):
    with open('{}/{}.revan.cfg'.format(output_dir, source_file), mode='w') as cfg:
        with open('./{}'.format(base_file)) as base:
            for line in base.readlines():
                if line.find('<GeometryFileName>') != -1:
                    cfg.write('<GeometryFileName>{}</GeometryFileName>\n'.format(geometry_file))
                elif line.find('<DataFileName>') != -1:
                    cfg.write('<DataFileName>{}/{}.inc1.id1.sim.gz</DataFileName>\n'.format(source_dir, source_file))
                else:
                    cfg.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g",  help="Geometry file")
    parser.add_argument("-f", default=os.getcwd(), help="Source file dir (Default: current directory)")
    parser.add_argument("-o", default=os.getcwd(), help="Output file dir (Default: current directory); will create symbolic links to input files if output directory is different from input")
    parser.add_argument("-b", default='revan_AMEGO_X.cfg', help="Base revan config file, default: revan_AMEGO_X.cfg")
    args = parser.parse_args()

    source_dir= args.f
    output_dir= args.o
    geometry_file = args.g
    base_file= args.b
    
    if not os.path.exists( output_dir):
        os.makedirs( output_dir )

    args = parser.parse_args()
    with open("{}/runRevan.sh".format(output_dir), mode='w') as f:
        for myene in energies:
            for cosTh,ang in zip(cos_ang,theta_angles):
                for phi in phi_angles:
                    source_file='%s_%.3fMeV_Cos%.1f_Phi%.1f'%(OneBeam,myene/1000.,cosTh, phi)
                    
                    create_config_file(source_dir, source_file, output_dir, geometry_file)
                    
                    simfile_in="{0}/{1}.inc1.id1.sim.gz".format(source_dir, source_file)
                    simfile_out="{0}/{1}.inc1.id1.sim.gz".format(output_dir, source_file)

                    if not os.path.exists( simfile_out ):
                        if not os.path.exists( simfile_in ):
                            print (simfile_in)
                        assert os.path.exists( simfile_in )
                        os.symlink(simfile_in, simfile_out)
                    
                    f.write("revan -a --oi -n -f {3} -g {2} -c {0}/{1}.revan.cfg >& {0}/log_{1}.txt \n".format(output_dir, source_file, geometry_file, simfile_out) )
                    
    print("Run this script: {}/runRevan.sh".format(output_dir))
                

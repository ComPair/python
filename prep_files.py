
# Author: Grant Sommer (grant.sommer@gwu.edu)
# Date: January 15, 2026

# Before running get_performance_plots.sh, in this code update:
#   energies, angles, Ntriggers, geofile, OneBeam, revan config file, and the main directory for these simulations
#
#
# This code is a wrapper for Prep_cosima.py, Prep_revan.py
# This code will return runCosima.sh and the .source files as well as runRevan.sh

import os
import sys
import Prep_cosima
import Prep_revan

### User Inputs Below ###

Log_E = [2.0, 3.0]
angles = [0]
N_Triggers=50000
geofile='/data/slag2/gsommer1/ComPair2/Simulation_Work/Geometry/ComPair_23/ComPair23.geo.setup'
# geofile = '/Users/gsommer1/Software/ComPair2/Geometry/ComPair_23/ComPair23.geo.setup'
OneBeam = 'FarFieldPointSource'
revan_config_file='/data/slag2/gsommer1/ComPair2/Simulation_Work/Automated/python/revan_ComPair2_Grant_1-14-26.cfg'
# revan_config_file='/Users/gsommer1/Software/ComPair2/Coding_Work/Automated_Simulations/python/revan_ComPair2_Grant_1-14-26.cfg'
working_directory='/data/slag2/gsommer1/ComPair2/Simulation_Work/Automated/python/1-16-26'
# working_directory='/Users/gsommer1/Software/ComPair2/Coding_Work/Automated_Simulations/python/AutomatedAnalysis_1-15-26_part_4'

### End User Inputs ###

# printing to capture in the main shell script
print(working_directory)

os.makedirs(f'{working_directory}/Plots', exist_ok=True)
os.makedirs(f'{working_directory}/Plots/FiguresOfMerit', exist_ok=True)

Prep_cosima.main(Log_E=Log_E, angles=angles, N_Triggers=N_Triggers, geofile=geofile, OneBeam=OneBeam, working_directory=working_directory)

Prep_revan.main(Log_E=Log_E, angles=angles, geofile=geofile, OneBeam=OneBeam, revan_config_file=revan_config_file, source_dir=working_directory)


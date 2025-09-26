There are a suit of python analysis tools which have been developed by the AMEGO team. Many of these duplicate the functionality of MEGAlib's program mimrec, but were written to allow for more flexibility and a better understanding for the user.

This description was initially written on confluence, but it's being included here years later since confluence access isn't easy: https://confluence.slac.stanford.edu/display/COMPair/AMEGO+Python+Tools

There was a lot of development of these tools years ago, and many of the scripts have become out of date but can be found on github using the tag Python2_OldAMEGOModel.

## Step-by-step: Using the ComPair python scripts for AMEGO Performance Plots

Prerequisites: 
- MEGAlib (http://megalibtoolkit.com/home.html)
- Geometry files (https://github.com/ComPair/Geometry)
- the Python scripts (https://github.com/ComPair/python)
*The paths in each of the commands below will need to be changed to point towards the files/directories on your computer


If you're not wanting to run the full simulations, and just use the Log Files from GitHub, then start at Step 5
# Step 1

First of all, we need to generate the source files. This can be done automatically by running Prep_cosima.py:
  $ python Prep_cosima.py

You can/should change the following parameters in Prep_cosima.py:
  Geofile = '/data/slag2/dtak/Geometry/AMEGO_Midex/AmegoBase.geo.setup'
  Log_E=[2.2,2.5,2.7,3,3.2,3.5,3.7,4,4.2,4.5,4.7,5,5.2,5.5,5.7,6,6.2,6.5,6.7]
  angles=[0,25.8,36.9,45.6,53.1,60]
  where "Geofile" is the address to the mass model, "Log_E" is the energy of mono-energetic sources to be simulated on a log scale (from 0.158 keV to 5011 MeV), and "angles" is the incident angle in degrees. For all of the simulations done in the past year (2019 - 2020), we have only looked at on-axis sources with angle=0.
  
Once all of the .source files and runCosima.sh have been created, use Cosima to perform the simulations. We've collated the .source files to run them all one after another in cosima with runCosima.sh, but this is cumbersome if you're doing it on your personal computer. You might have better luck starting them one after the other, or we may provide (or you can write your own to share!) a more sophisticated script soon...
  $ source runCosima.sh

# Step 2
After the .sim files have finished, the next step is to perform the event reconstruction with revan. You can use Prep_revan.py to do this automatically, which generates a runRevan.sh.
   $ python Prep_revan.py /data/slag2/dtak/Geometry/AMEGO_Midex/AmegoBase.geo.setup -b revan_AMEGO_X.cfg
   $ source runRevan.sh

   The revan_AMEGO_X.cfg file contains the configuration setting for revan and is included in the GitHub python repo.

# Step 3
After new simulations have been run and the revan reconstruction has been performed, you’re ready to run the AMEGO python tools to create instrument performance plots! The first thing to scan the .sim files to create the TriggerEfficiency.txt file, which is used in subsequent programs to calculate the effective area.

Start a python session:
  $ python
  > import EventAnalysis
  > EventAnalysis.getTriggerEfficiency(directory='/data/slag2/rcaputo/AMEGO/AMEGO4x4PerformancePlotSimFiles/')
where the directory here is the location of the .sim files. This creates a TriggerEfficiency.txt file in the same directory as the .sim files.

Then, create the .log files for each .tra file with the performCompleteAnalysis function. This fits the energy spectra and angular resolution for each energy and angle. The option showPlots should be used if you’re planning to use these values, because often times the energy and ARM fits suck for the Compton range. To ensure convergence, at this point many of these fits have been manually performed in mimrec (and then we've manually edited the .log file). *We should implement an energy cut before the ARM plot in this python program.

Create .log files for all .tra files in a directory:
  > EventAnalysis.performCompleteAnalysis(directory=‘/data/slag2/rcaputo/AMEGO/AMEGO4x4PerformancePlotTraFiles/TraFiles_R1/‘, showPlots=True)

Or for a single file:
  > EventAnalysis.performCompleteAnalysis(filename=‘/data/slag2/rcaputo/AMEGO/AMEGO4x4PerformancePlotTraFiles/TraFiles_R5/FarFieldPointSource_3162.277MeV_Cos1.0.inc1.id1.tra ‘, showPlots=True)

At this point, we recommend you to run the python program for a single file at a time since you may need to change the range of the angle and energy band to make sure the Gaussian fit converged. (Another option is to perform these fits in MEGAlib, and in that case, I have manually edited each .log file in the Compton regime to include the numbers I achieved from mimrec fits.) You may also need to remove some unreasonable numbers in the .log files; e.g., fit from low photon statistics. *performCompleteAnalysis needs to be made more robust so this isn't an issue.

# Step 4
Then, to make the performance plots! (The plots below are from the simulations performed March 2019 for AMEGO-X)
  > import FigureOfMeritPlotter
  > triggerEfficiencyFilename ='/data/slag2/rcaputo/AMEGO/AMEGO4x4PerformancePlotSimFiles/TriggerEfficiency.txt'

Load the directory for log files and the TriggerEfficiency.txt file
  > directory='/data/slag2/rcaputo/AMEGO/AMEGO4x4PerformancePlotTraFiles/LogFilesFinalAnalyis/'
  > data = FigureOfMeritPlotter.parseEventAnalysisLogs(directory, triggerEfficiencyFilename=triggerEfficiencyFilename)
  > FigureOfMeritPlotter.plotAngularResolution(data, angleSelections=1.0)

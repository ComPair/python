#!/usr/bin/env python
"""
------------------------------------------------------------------------

Script to parse and visualize results created by mimrec:

Author: Daniel Kocevski (dankocevski@gmail.com)
Date: August 25st, 2016

Usage Examples:

import FigureOfMeritPlotter

simulationIDs = '../Simulations/PerformancePlotSourceFiles/FarFieldPointSourceIDs.txt'
data = FigureOfMeritPlotter.parse()

# Effecitive Area vs Energy
FigureOfMeritPlotter.plotEffectiveArea(data)

# Effecitive Area vs Angle
FigureOfMeritPlotter.plotEffectiveAreaVsAngle(data)

# Energy Resolution vs Energy
FigureOfMeritPlotter.plotEnergyResolution(data)

# Energy Resolution vs Angle
FigureOfMeritPlotter.plotEnergyResolutionVsAngle(data)

# Angular Resolution vs Energy
FigureOfMeritPlotter.plotAngularResolution(data)

# Angular Resolution vs Angle
FigureOfMeritPlotter.plotAngularResolutionVsAngle(data)

# Sensitivity Curves
FigureOfMeritPlotter.plotSourceSensitivity(data)


------------------------------------------------------------------------
"""

import fileinput
import glob
import sys
import os 
import numpy
import matplotlib.pylab as plot
import math
from astropy.io import ascii


##########################################################################################

def omega(PSF):

	omega_solidAngle = 2*math.pi*(1-numpy.cos(2*PSF*math.pi/180.))

    return omega_solidAngle

##########################################################################################

def Isrc(E,time, Aeff, nsig, domega, Bbkg):
    
    # print E
    # print time
    # print Aeff
    # print nsig
    # print domega
    # print Bbkg
    arg = numpy.sqrt((nsig**4/4.)+(nsig**2*Bbkg*Aeff*time*domega/E))
    num = E/(Aeff*time)*(nsig**2/2+arg)
    return num


##########################################################################################

def parseSimLog(filename):

	for line in fileinput.input([filename]):

		# Parse the lines
		if 'ID 100000 ' in line:
			NumberOfSimulatedEvents = float(line.split()[2])

	fileinput.close()

	return NumberOfSimulatedEvents


##########################################################################################

def parseMimrecLog(filename, interactionType, verbose=False):

	# Return nan values if the file doesn't exist
	if os.path.isfile(filename) == False:
		if 'Pair' in interactionType:
			return numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan
		else:
			return numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan

	RMS = []
	FWHM = []
	Containment68 = None

	# Read the file line by line
	for line in fileinput.input([filename]):

		# Parse the lines
		if '1  Constant' in line:
			Constant = float(line.split()[2])
			ConstantError = float(line.split()[3])	

		if '2  Mean' in line:
			Mean = float(line.split()[2])
			MeanError = float(line.split()[3])

		if '3  Sigma' in line:
			Sigma = float(line.split()[2])
			SigmaError = float(line.split()[3])

		if 'RMS: ' in line:
			RMS.append(float(line.split()[1]))

		if 'Pair' in interactionType:

			if '0.68% containment:' in line and Containment68 == None:
				Containment68 = float(line.split()[2])

			if 'Pair events in histogram:' in line:
				NumberOfReconstructedEvents = float(line.split()[4])

		else:

			if 'Total FWHM of fit (not of data!): ' in line:
				FWHM.append(float(line.split()[7]))

			if 'Compton and pair events in histogram: ' in line:
				NumberOfReconstructedEvents = float(line.split()[6])

	# Close the input
	fileinput.close()

	# Return the results
	if 'Pair' in interactionType:
		return Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, Containment68, NumberOfReconstructedEvents
	else:
		return Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS[-2], FWHM[-2], NumberOfReconstructedEvents


##########################################################################################

def parse(sumlationsIDs=None):

	# Define the path to the python repository
	pythonumpyath = os.path.dirname(os.path.realpath(__file__))

	# Define the root path containg the ComPair repositories
	rootPath = pythonumpyath.replace('python','')

	# Define the PerformancePlotSourceFiles path
	sourceDirectory = rootPath + 'Simulations/PerformancePlotSourceFiles'

	# Define the PerformancePlotMimrecFiles path
	mimrecDirectory = rootPath + 'Simulations/PerformancePlotMimrecFiles'

	# Define the PerformancePlotTraFiles path
	traDirectory = rootPath + 'Simulations/PerformancePlotTraFiles' 

	if sumlationsIDs == None:
		sumlationsIDs = sourceDirectory + '/FarFieldPointSourceIDs.txt'

	# Create an dictionary to store the data
	data = {}


	# Read the sim file
	with open(sumlationsIDs) as filehandle:
		lines = filehandle.readlines()

	currentNumber = 1
	print '\numpyarsing mimrec logs...'

	# Loop through each of the lines
	for line in lines:
		lineContents = line.split()
		simulationName = lineContents[0].split('.inc')[0]		
		simulationFilename = lineContents[0].replace(':ID','')
		NumberOfSimulatedEvents = lineContents[2]

		# Enter an empty list for the initial value to which the results will be appended
		data[simulationName] = []

		# Add the number of simulated to the results dictionary
		data[simulationName].append(NumberOfSimulatedEvents)

		# Generate the mimrec log names
		mimrecFilename_tracked = mimrecDirectory + '/'  + simulationFilename.replace('.sim', '.mimrec_tracked.log')
		mimrecFilename_untracked = mimrecDirectory + '/'  + simulationFilename.replace('.sim', '.mimrec_untracked.log')
		mimrecFilename_pair = mimrecDirectory + '/'  + simulationFilename.replace('.sim', '.mimrec_pair.log')

		# Get the mimrec figures of merit and add the dictionary
		# print "Parsing: %s" % mimrecFilename_tracked
		Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents = parseMimrecLog(mimrecFilename_tracked, interactionType='Compton')
		data[simulationName].append([Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents])

		# Get the mimrec figures of merit and add the dictionary
		# print "Parsing: %s" % mimrecFilename_untracked	
		Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents = parseMimrecLog(mimrecFilename_untracked, interactionType='Compton')
		data[simulationName].append([Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents])

		# Get the mimrec figures of merit and add the dictionary
		# print "Parsing: %s" % mimrecFsilename_pair		
		Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, Containment68, NumberOfReconstructedPairEvents = parseMimrecLog(mimrecFilename_pair, interactionType='Pair')
		data[simulationName].append([Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, Containment68, NumberOfReconstructedPairEvents])

		currentNumber = currentNumber + 1


	print 'Done.'
	return data

##########################################################################################

def plotAngularResolution(data, angleSelections=[1,0.9,0.8,0.7,0.6,0.5], xlog=True, ylog=False, save=False):
	
	if hasattr(angleSelections, '__iter__') == False:
		angleSelections = [angleSelections]


	plotNumber = 1
	plot.figure(figsize=(10,12))

	for angleSelection in angleSelections:

		Energy = []
		FWHM_tracked = []
		FWHM_untracked = []
		Containment68 = []

		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			angle = float(key.split('_')[2].replace('Cos',''))

			if angle == angleSelection:
				Energy.append(energy)
				fwhm_tracked = data[key][1][7]
				fwhm_untracked = data[key][2][7]
				containment68 = data[key][3][6]

				FWHM_tracked.append(fwhm_tracked)
				FWHM_untracked.append(fwhm_untracked)
				Containment68.append(containment68)

		# Convert everything to a numpy array
		Energy = numpy.array(Energy)
		FWHM_tracked = numpy.array(FWHM_tracked)
		FWHM_untracked = numpy.array(FWHM_untracked)
		Containment68 = numpy.array(Containment68)

		# Sort by energy
		i = [numpy.argsort(Energy)]
		Energy = Energy[i]
		FWHM_tracked = FWHM_tracked[i]
		FWHM_untracked = FWHM_untracked[i]
		Containment68 = Containment68[i]

		# Plot the data
		ax = plot.subplot( str(len(angleSelections)) + str(10 + plotNumber) )

		plot.scatter(Energy, FWHM_tracked, color='darkgreen')
		plot.plot(Energy, FWHM_tracked, color='darkgreen', alpha=0.5, label='Compton (tracked)')

		plot.scatter(Energy, FWHM_untracked, color='blue')
		plot.plot(Energy, FWHM_untracked, color='blue', alpha=0.5, label='Compton (untracked)')

		plot.scatter(Energy, Containment68, color='darkred')
		plot.plot(Energy, Containment68, color='darkred', alpha=0.5, label='Pair')

		if plotNumber == 1:
			plot.title('Angular Resolution')			
			plot.legend(numpoints=1, scatterpoints=1, fontsize='small', frameon=True, loc='upper right')

		if xlog == True:
			plot.xscale('log')

		if ylog == True:
			plot.yscale('log')

		plot.ylabel(u'fwhm (deg)')

		plot.text(0.015, 0.8, u'%i\N{DEGREE SIGN}' % round(numpy.degrees(numpy.arccos(angleSelection))),
		        verticalalignment='bottom', horizontalalignment='left',
		        transform=ax.transAxes,
		        color='black', fontsize=12)


		if plotNumber != len(angleSelections):
			ax.set_xticklabels([])

			# labels = ax.get_yticklabels()
			# print labels
			# labels[0] = ""
			# labels[-1] = ""
			# ax.set_yticklabels(labels)


		if plotNumber == len(angleSelections):
			plot.xlabel(r'Energy (MeV)')


		plotNumber = plotNumber + 1


	plot.subplots_adjust(wspace=0, hspace=.2)

	if save == True:
		plot.savefig('AngularResolution.png', bbox_inches='tight')

	plot.show()

	plot.close()

##########################################################################################

def plotAngularResolutionVsAngle(data, energySelections=[0.3, 1.0, 3.16, 10.0, 31.6, 100, 316.0], xlog=False, ylog=False, save=False):
	
	if hasattr(energySelections, '__iter__') == False:
		energySelections = [energySelections]

	plotNumber = 1
	plot.figure(figsize=(10,12))

	for energySelection in energySelections:

		Angle = []
		FWHM_tracked = []
		FWHM_untracked = []
		Containment68 = []

		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			angle = float(key.split('_')[2].replace('Cos',''))
			angle = round(numpy.degrees(numpy.arccos(angle)))

			if energy == energySelection:
				Angle.append(angle)
				fwhm_tracked = data[key][1][7]
				fwhm_untracked = data[key][2][7]
				containment68 = data[key][3][6]

				FWHM_tracked.append(fwhm_tracked)
				FWHM_untracked.append(fwhm_untracked)
				Containment68.append(containment68)

		# Convert everything to a numpy array
		Angle = numpy.array(Angle)
		FWHM_tracked = numpy.array(FWHM_tracked)
		FWHM_untracked = numpy.array(FWHM_untracked)
		Containment68 = numpy.array(Containment68)

		# Sort by energy
		i = [numpy.argsort(Angle)]
		Angle = Angle[i]
		FWHM_tracked = FWHM_tracked[i]
		FWHM_untracked = FWHM_untracked[i]
		Containment68 = Containment68[i]

		# Plot the data
		ax = plot.subplot( str(len(energySelections)) + str(10 + plotNumber) )

		plot.scatter(Angle, FWHM_tracked, color='darkgreen')
		plot.plot(Angle, FWHM_tracked, color='darkgreen', alpha=0.5, label='Compton (tracked)')

		plot.scatter(Angle, FWHM_untracked, color='blue')
		plot.plot(Angle, FWHM_untracked, color='blue', alpha=0.5, label='Compton (untracked)')

		plot.scatter(Angle, Containment68, color='darkred')
		plot.plot(Angle, Containment68, color='darkred', alpha=0.5, label='Pair')

		if plotNumber == 1:
			plot.title('Angular Resolution')
			plot.legend(numpoints=1, scatterpoints=1, fontsize='small', frameon=True, loc='upper right')

		if xlog == True:
			plot.xscale('log')

		if ylog == True:
			plot.yscale('log')

		plot.ylabel(u'fwhm (deg)')

		plot.text(0.015, 0.8, '%s MeV' % energySelection,
		        verticalalignment='bottom', horizontalalignment='left',
		        transform=ax.transAxes,
		        color='black', fontsize=12)


		if plotNumber != len(energySelections):
			ax.set_xticklabels([])

			# labels = ax.get_yticklabels()
			# print labels
			# labels[0] = ""
			# labels[-1] = ""
			# ax.set_yticklabels(labels)


		if plotNumber == len(energySelections):
			plot.xlabel(r'$\theta$')


		plotNumber = plotNumber + 1


	plot.subplots_adjust(wspace=0, hspace=.2)

	if save == True:
		plot.savefig('AngularResolutionVsAngle.png', bbox_inches='tight')

	plot.show()

	plot.close()


##########################################################################################

def plotEnergyResolution(data, angleSelections=[1,0.9,0.8,0.7,0.6,0.5], xlog=True, ylog=False, save=False):
	
	if hasattr(angleSelections, '__iter__') == False:
		angleSelections = [angleSelections]

	plotNumber = 1
	plot.figure(figsize=(10,12))

	for angleSelection in angleSelections:

		Energy = []
		Sigma_tracked = []
		Sigma_untracked = []
		Sigma_pair = []

		SigmaError_tracked = []
		SigmaError_untracked = []
		SigmaError_pair = []

		# print 'Name, fwhm_tracked, fwhm_untracked, containment'

		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			angle = float(key.split('_')[2].replace('Cos',''))

			if angle == angleSelection:
				Energy.append(energy)

				Sigma_tracked.append(data[key][1][4])
				Sigma_untracked.append(data[key][2][4])
				Sigma_pair.append(data[key][3][4])

				SigmaError_tracked.append(data[key][1][5])
				SigmaError_untracked.append(data[key][2][5])
				SigmaError_pair.append(data[key][3][5])


		# Convert everything to a numpy array
		Energy = numpy.array(Energy)
		Sigma_tracked = numpy.array(Sigma_tracked)
		Sigma_untracked = numpy.array(Sigma_untracked)
		Sigma_pair = numpy.array(Sigma_pair)

		SigmaError_tracked = numpy.array(SigmaError_tracked)
		SigmaError_untracked = numpy.array(SigmaError_untracked)
		SigmaError_pair = numpy.array(SigmaError_pair)

		# Sort by energy
		i = [numpy.argsort(Energy)]
		Energy = Energy[i]
		Sigma_tracked = Sigma_tracked[i]
		Sigma_untracked = Sigma_untracked[i]
		Sigma_pair = Sigma_pair[i]
		SigmaError_tracked = SigmaError_tracked[i]
		SigmaError_untracked = SigmaError_untracked[i]
		SigmaError_pair = SigmaError_pair[i]


		# Plot the data
		ax = plot.subplot( str(len(angleSelections)) + str(10 + plotNumber) )

		plot.scatter(Energy, Sigma_tracked, color='darkgreen')
		plot.errorbar(Energy, Sigma_tracked, yerr=SigmaError_tracked, color='darkgreen', fmt=None)		
		plot.plot(Energy, Sigma_tracked, color='darkgreen', alpha=0.5, label='Compton (tracked)')

		plot.scatter(Energy, Sigma_untracked, color='blue')
		plot.errorbar(Energy, Sigma_untracked, yerr=SigmaError_untracked, color='blue', fmt=None)				
		plot.plot(Energy, Sigma_untracked, color='blue', alpha=0.5, label='Compton (untracked)')

		plot.scatter(Energy, Sigma_pair, color='darkred')
		plot.errorbar(Energy, Sigma_pair, yerr=SigmaError_pair, color='darkred', fmt=None)		
		plot.plot(Energy, Sigma_pair, color='darkred', alpha=0.5, label='Pair')

		if plotNumber == 1:
			plot.title('Energy Resolution')			
			plot.legend(numpoints=1, scatterpoints=1, fontsize='small', frameon=True, loc='upper left')

		if xlog == True:
			plot.xscale('log')

		if ylog == True:
			plot.yscale('log')

		plot.ylabel(r'$\sigma$')

		plot.text(1-0.015, 0.8, u'%i\N{DEGREE SIGN}' % round(numpy.degrees(numpy.arccos(angleSelection))),
		        verticalalignment='bottom', horizontalalignment='right',
		        transform=ax.transAxes,
		        color='black', fontsize=12)


		if plotNumber != len(angleSelections):
			ax.set_xticklabels([])

		if plotNumber == len(angleSelections):
			plot.xlabel(r'Energy (MeV)')

		plotNumber = plotNumber + 1


	plot.subplots_adjust(wspace=0, hspace=.2)

	if save == True:
		plot.savefig('EnergyResolution.png', bbox_inches='tight')

	plot.show()

	plot.close()


##########################################################################################

def plotEnergyResolutionVsAngle(data, energySelections=[0.3, 1.0, 3.16, 10.0, 31.6, 100, 316.0], xlog=False, ylog=False, save=False):
	
	if hasattr(energySelections, '__iter__') == False:
		energySelections = [energySelections]

	plotNumber = 1
	plot.figure(figsize=(10,12))

	for energySelection in energySelections:

		Angle = []
		Sigma_tracked = []
		Sigma_untracked = []
		Sigma_pair = []

		SigmaError_tracked = []
		SigmaError_untracked = []
		SigmaError_pair = []

		# print 'Name, fwhm_tracked, fwhm_untracked, containment'

		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			angle = float(key.split('_')[2].replace('Cos',''))

			if energy == energySelection:
				Angle.append(angle)

				Sigma_tracked.append(data[key][1][4])
				Sigma_untracked.append(data[key][2][4])
				Sigma_pair.append(data[key][3][4])

				SigmaError_tracked.append(data[key][1][5])
				SigmaError_untracked.append(data[key][2][5])
				SigmaError_pair.append(data[key][3][5])


		# Convert everything to a numpy array
		Angle = numpy.array(Angle)
		Sigma_tracked = numpy.array(Sigma_tracked)
		Sigma_untracked = numpy.array(Sigma_untracked)
		Sigma_pair = numpy.array(Sigma_pair)

		SigmaError_tracked = numpy.array(SigmaError_tracked)
		SigmaError_untracked = numpy.array(SigmaError_untracked)
		SigmaError_pair = numpy.array(SigmaError_pair)

		# Sort by energy
		i = [numpy.argsort(Angle)]
		Angle = Angle[i]
		Sigma_tracked = Sigma_tracked[i]
		Sigma_untracked = Sigma_untracked[i]
		Sigma_pair = Sigma_pair[i]
		SigmaError_tracked = SigmaError_tracked[i]
		SigmaError_untracked = SigmaError_untracked[i]
		SigmaError_pair = SigmaError_pair[i]


		# Plot the data
		ax = plot.subplot( str(len(energySelections)) + str(10 + plotNumber) )

		plot.scatter(Angle, Sigma_tracked, color='darkgreen')
		plot.errorbar(Angle, Sigma_tracked, yerr=SigmaError_tracked, color='darkgreen', fmt=None)		
		plot.plot(Angle, Sigma_tracked, color='darkgreen', alpha=0.75, label='Compton (tracked)')

		plot.scatter(Angle, Sigma_untracked, color='blue')
		plot.errorbar(Angle, Sigma_untracked, yerr=SigmaError_untracked, color='blue', fmt=None)				
		plot.plot(Angle, Sigma_untracked, color='blue', alpha=0.75, label='Compton (untracked)')

		plot.scatter(Angle, Sigma_pair, color='darkred')
		plot.errorbar(Angle, Sigma_pair, yerr=SigmaError_pair, color='darkred', fmt=None)		
		plot.plot(Angle, Sigma_pair, color='darkred', alpha=0.75, label='Pair')

		if plotNumber == 1:
			plot.title('Energy Resolution')						
			plot.legend(numpoints=1, scatterpoints=1, fontsize='small', frameon=True, loc='upper left')

		if xlog == True:
			plot.xscale('log')

		if ylog == True:
			plot.yscale('log')

		plot.ylabel(r'$\sigma$')

		plot.text(1-0.015, 0.8, '%s MeV' % energySelection,
		        verticalalignment='bottom', horizontalalignment='right',
		        transform=ax.transAxes,
		        color='black', fontsize=12)

		if plotNumber != len(energySelections):
			ax.set_xticklabels([])

		if plotNumber == len(energySelections):
			plot.xlabel(r'$\theta$')

		plotNumber = plotNumber + 1


	plot.subplots_adjust(wspace=0, hspace=.2)

	if save == True:
		plot.savefig('EnergyResolutionVsAngle.png', bbox_inches='tight')

	plot.show()

	plot.close()

##########################################################################################

def plotEffectiveArea(data, angleSelections=[1,0.9,0.8,0.7,0.6,0.5], xlog=True, ylog=False, save=False):

	if hasattr(angleSelections, '__iter__') == False:
		angleSelections = [angleSelections]

	plotNumber = 1
	plot.figure(figsize=(10,12))

	for angleSelection in angleSelections:

		Energy = []
		EffectiveArea_Tracked = []
		EffectiveArea_Untracked  = []
		EffectiveArea_Pair  = []

		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			angle = float(key.split('_')[2].replace('Cos',''))

			if angle == angleSelection:
				numberOfSimulatedEvents = float(data[key][0])
				numberOfReconstructedEvents_tracked = float(data[key][1][-1])
				numberOfReconstructedEvents_untracked = float(data[key][2][-1])
				numberOfReconstructedEvents_pair = float(data[key][3][-1])

				# Calculate the effective area
				effectiveArea_tracked = (numberOfReconstructedEvents_tracked/numberOfSimulatedEvents) * math.pi * 300**2
				effectiveArea_untracked = (numberOfReconstructedEvents_untracked/numberOfSimulatedEvents) * math.pi * 300**2
				effectiveArea_pair = (numberOfReconstructedEvents_pair/numberOfSimulatedEvents) * math.pi * 300**2

				# Store the results
				Energy.append(energy)
				EffectiveArea_Tracked.append(effectiveArea_tracked)
				EffectiveArea_Untracked.append(effectiveArea_untracked)
				EffectiveArea_Pair.append(effectiveArea_pair)


		# Convert everything to a numpy array
		Energy = numpy.array(Energy)
		EffectiveArea_Tracked = numpy.array(EffectiveArea_Tracked)
		EffectiveArea_Untracked = numpy.array(EffectiveArea_Untracked)
		EffectiveArea_Pair = numpy.array(EffectiveArea_Pair)

		# Sort by energy
		i = [numpy.argsort(Energy)]
		Energy = Energy[i]
		EffectiveArea_Tracked = EffectiveArea_Tracked[i]
		EffectiveArea_Untracked = EffectiveArea_Untracked[i]
		EffectiveArea_Pair = EffectiveArea_Pair[i]


		# Create the new subplot
		ax = plot.subplot( str(len(angleSelections)) + str(10 + plotNumber) )

		# Plot the data
		plot.scatter(Energy, EffectiveArea_Tracked, color='darkgreen')
		plot.plot(Energy, EffectiveArea_Tracked, color='darkgreen', alpha=0.5, label='Compton (tracked)')

		plot.scatter(Energy, EffectiveArea_Untracked, color='blue')
		plot.plot(Energy, EffectiveArea_Untracked, color='blue', alpha=0.5, label='Compton (untracked)')

		plot.scatter(Energy, EffectiveArea_Pair, color='darkred')
		plot.plot(Energy, EffectiveArea_Pair, color='darkred', alpha=0.5, label='Pair')	


		plot.ylabel(r'A$_{\mathrm{eff}}$ (cm$^2$)')

		plot.text(1-0.015, 0.8, u'%i\N{DEGREE SIGN}' % round(numpy.degrees(numpy.arccos(angleSelection))),
		        verticalalignment='bottom', horizontalalignment='right',
		        transform=ax.transAxes,
		        color='black', fontsize=12)

		if plotNumber == 1:
			plot.title('Effective Area')			

		if xlog == True:
			plot.xscale('log')

		if ylog == True:
			plot.yscale('log')


		if plotNumber != len(angleSelections):
			ax.set_xticklabels([])

		if plotNumber == len(angleSelections):
			plot.xlabel(r'Energy (MeV)')

		plotNumber = plotNumber + 1


	plot.subplots_adjust(wspace=0, hspace=.2)

	if save == True:
		plot.savefig('EffectiveArea.png', bbox_inches='tight')

	plot.show()


 ##########################################################################################

def plotEffectiveAreaVsAngle(data, energySelections=[0.3, 1.0, 3.16, 10.0, 31.6, 100, 316.0], xlog=False, ylog=False, save=False, collapse=False):

	if hasattr(energySelections, '__iter__') == False:
		energySelections = [energySelections]

	plotNumber = 1

	# Create the new subplot
	if collapse == False:
		plot.figure(figsize=(10,12))
	else:
		plot.figure(figsize=(10, 6.39))
		ax = plot.subplot(111)
		colors = ['purple', 'blue', 'green', 'orange', 'brown', 'red', 'darkred']

	for energySelection in energySelections:

		Angle = []
		EffectiveArea_Tracked = []
		EffectiveArea_Untracked  = []
		EffectiveArea_Pair  = []

		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			angle = float(key.split('_')[2].replace('Cos',''))
			angle = round(numpy.degrees(numpy.arccos(angle)))

			if energy == energySelection:
				numberOfSimulatedEvents = float(data[key][0])
				numberOfReconstructedEvents_tracked = float(data[key][1][-1])
				numberOfReconstructedEvents_untracked = float(data[key][2][-1])
				numberOfReconstructedEvents_pair = float(data[key][3][-1])

				# Calculate the effective area
				effectiveArea_tracked = (numberOfReconstructedEvents_tracked/numberOfSimulatedEvents) * math.pi * 300**2
				effectiveArea_untracked = (numberOfReconstructedEvents_untracked/numberOfSimulatedEvents) * math.pi * 300**2
				effectiveArea_pair = (numberOfReconstructedEvents_pair/numberOfSimulatedEvents) * math.pi * 300**2

				# Store the results
				Angle.append(angle)
				EffectiveArea_Tracked.append(effectiveArea_tracked)
				EffectiveArea_Untracked.append(effectiveArea_untracked)
				EffectiveArea_Pair.append(effectiveArea_pair)


		# Convert everything to a numpy array
		Angle = numpy.array(Angle)
		EffectiveArea_Tracked = numpy.array(EffectiveArea_Tracked)
		EffectiveArea_Untracked = numpy.array(EffectiveArea_Untracked)
		EffectiveArea_Pair = numpy.array(EffectiveArea_Pair)

		# Sort by energy
		i = [numpy.argsort(Angle)]
		Angle = Angle[i]
		EffectiveArea_Tracked = EffectiveArea_Tracked[i]
		EffectiveArea_Untracked = EffectiveArea_Untracked[i]
		EffectiveArea_Pair = EffectiveArea_Pair[i]

		# Create the new subplot
		if collapse == False:
			ax = plot.subplot( str(len(energySelections)) + str(10 + plotNumber) )

		# Plot the data
		if collapse == False:
			plot.scatter(Angle, EffectiveArea_Tracked, color='darkgreen')
			plot.plot(Angle, EffectiveArea_Tracked, color='darkgreen', alpha=0.5, label='Compton (tracked)')

			plot.scatter(Angle, EffectiveArea_Untracked, color='blue')
			plot.plot(Angle, EffectiveArea_Untracked, color='blue', alpha=0.5, label='Compton (untracked)')

			plot.scatter(Angle, EffectiveArea_Pair, color='darkred')
			plot.plot(Angle, EffectiveArea_Pair, color='darkred', alpha=0.5, label='Pair')			

		else:
			plot.scatter(Angle, EffectiveArea_Tracked, color=colors[plotNumber-1])
			plot.plot(Angle, EffectiveArea_Tracked, color=colors[plotNumber-1], alpha=0.5)

			plot.scatter(Angle, EffectiveArea_Untracked, color=colors[plotNumber-1])
			plot.plot(Angle, EffectiveArea_Untracked, color=colors[plotNumber-1], alpha=0.5, linestyle='--')

			plot.scatter(Angle, EffectiveArea_Pair, color=colors[plotNumber-1])
			plot.plot(Angle, EffectiveArea_Pair, color=colors[plotNumber-1], alpha=0.5, linestyle='-.')



		if plotNumber == 1:
			plot.title('Effective Area')			
			plot.legend(numpoints=1, scatterpoints=1, fontsize='small', frameon=True, loc='upper left')

		plot.ylabel(r'A$_{\mathrm{eff}}$ (cm$^2$)')

		if collapse == False:
			plot.text(1-0.015, 0.8, '%s MeV' % energySelection,
			        verticalalignment='bottom', horizontalalignment='right',
		       		transform=ax.transAxes,
			        color='black', fontsize=12)
		# else:
		# 	plot.text(61, EffectiveArea_Tracked[-1], '%s MeV' % energySelection,
		# 	        verticalalignment='bottom', horizontalalignment='left',
		# 	        color='black', fontsize=12)

		if xlog == True:
			plot.xscale('log')

		if ylog == True:
			plot.yscale('log')


		if plotNumber != len(energySelections) and collapse == False:
			ax.set_xticklabels([])

		if plotNumber == len(energySelections) or collapse == True:
			plot.xlabel(r'$\theta$')

		plotNumber = plotNumber + 1


	plot.subplots_adjust(wspace=0, hspace=.2)

	if save == True:
		plot.savefig('EffectiveAreaVsAngle.png', bbox_inches='tight')

	plot.show()


########################################################################################## 

def plotSourceSensitivity(data, angleSelection=0.7, xlog=True, ylog=True, save=False, doplot=False):

	# Energy = numpy.array([0.316, 1.0, 3.16, 10., 31.6, 100., 316.])
	# energy = numpy.array([0.316, 1.0, 3.16, 10., 31.6, 100., 316.])

	A_eff = numpy.array([10., 10., 10., 10., 10., 10., 10.])
	background = numpy.array([0.00346008, 0.00447618, 0.00594937, 0.00812853, 0.0100297, 0.0124697, 0.0161290])
	exposure = 6.3*10**6

	Energy = []

	EffectiveArea_Tracked = []
	EffectiveArea_Untracked  = []
	EffectiveArea_Pair  = []
	EffectiveArea_Pair_ideal  = []

	FWHM_tracked = []
	FWHM_untracked = []
	Containment68 = []


	for key in data.keys():

		energy = float(key.split('_')[1].replace('MeV',''))
		angle = float(key.split('_')[2].replace('Cos',''))

		if angle == angleSelection:

			# Get the number of reconstructed events
			numberOfSimulatedEvents = float(data[key][0])
			numberOfReconstructedEvents_tracked = float(data[key][1][-1])
			numberOfReconstructedEvents_untracked = float(data[key][2][-1])
			numberOfReconstructedEvents_pair = float(data[key][3][-1])
			numberOfReconstructedEvents_pair_ideal = 100000.
            
			# Get the angular resolution
			fwhm_tracked = data[key][1][7]
			fwhm_untracked = data[key][2][7]
			containment68 = data[key][3][6]

			# Calculate the effective area
			effectiveArea_tracked = (numberOfReconstructedEvents_tracked/numberOfSimulatedEvents) * math.pi * 300**2
			effectiveArea_untracked = (numberOfReconstructedEvents_untracked/numberOfSimulatedEvents) * math.pi * 300**2
			effectiveArea_pair = (numberOfReconstructedEvents_pair/numberOfSimulatedEvents) * math.pi * 300**2
			effectiveArea_pair_ideal = (numberOfReconstructedEvents_pair_ideal/numberOfSimulatedEvents) * math.pi * 300**2

			# Store the effective area results
			Energy.append(energy)
			EffectiveArea_Tracked.append(effectiveArea_tracked)
			EffectiveArea_Untracked.append(effectiveArea_untracked)
			EffectiveArea_Pair.append(effectiveArea_pair)
			EffectiveArea_Pair_ideal.append(effectiveArea_pair_ideal)

			# Store the angular resolution results
			FWHM_tracked.append(fwhm_tracked)
			FWHM_untracked.append(fwhm_untracked)
			Containment68.append(containment68)

	# Convert everything to a numpy array
	Energy = numpy.array(Energy)
	EffectiveArea_Tracked = numpy.array(EffectiveArea_Tracked)
	EffectiveArea_Untracked = numpy.array(EffectiveArea_Untracked)
	EffectiveArea_Pair = numpy.array(EffectiveArea_Pair)
	EffectiveArea_Pair_ideal = numpy.array(EffectiveArea_Pair_ideal)

	FWHM_tracked = numpy.array(FWHM_tracked)
	FWHM_untracked = numpy.array(FWHM_untracked)
	Containment68 = numpy.array(Containment68)

	# Sort by energy
	i = [numpy.argsort(Energy)]
	Energy = Energy[i]
	EffectiveArea_Tracked = EffectiveArea_Tracked[i]
	EffectiveArea_Untracked = EffectiveArea_Untracked[i]
	EffectiveArea_Pair = EffectiveArea_Pair[i]
	EffectiveArea_Pair_ideal = EffectiveArea_Pair_ideal[i]

	FWHM_tracked = FWHM_tracked[i]
	FWHM_untracked = FWHM_untracked[i]
	Containment68 = Containment68[i]
	Containment68_idealang = 1.0
   
	Sensitivity_tracked = Isrc(Energy, exposure, EffectiveArea_Tracked, 3., omega(FWHM_tracked), background)
	Sensitivity_untracked = Isrc(Energy, exposure, EffectiveArea_Untracked, 3., omega(FWHM_untracked), background)
	Sensitivity_pair = Isrc(Energy, exposure, EffectiveArea_Pair, 3., omega(Containment68), background)
	Sensitivity_pair_ideal = Isrc(Energy, exposure, EffectiveArea_Pair_ideal, 3., omega(Containment68), background)
	Sensitivity_pair_idealang = Isrc(Energy, exposure, EffectiveArea_Pair_ideal, 3., omega(Containment68_idealang), background)

	plot.figure(figsize=(10, 6.39))
	ax = plot.subplot(111)

	plot.scatter(Energy, Sensitivity_tracked, color='darkgreen')
	plot.plot(Energy, Sensitivity_tracked, color='darkgreen', alpha=0.5, label='Compton (tracked)')

	plot.scatter(Energy, Sensitivity_untracked, color='blue')
	plot.plot(Energy, Sensitivity_untracked, color='blue', alpha=0.5, label='Compton (untracked)')

	plot.scatter(Energy, Sensitivity_pair, color='darkred')
	plot.plot(Energy, Sensitivity_pair, color='darkred', alpha=0.5, label='Pair')	

	if xlog == True:
		plot.xscale('log')

	if ylog == True:
		plot.yscale('log')

	plot.legend(numpoints=1, scatterpoints=1, fontsize='small', frameon=True, loc='upper left')

	if save == True:
		plot.savefig('SourceSensitivity.png', bbox_inches='tight')

	if doplot == True:
		plot.show()

	return Energy, Sensitivity_tracked, Sensitivity_untracked, Sensitivity_pair, Sensitivity_pair_ideal, Sensitivity_pair_idealang


##########################################################################################

def plotAllSourceSensitivities(data, angleSelection=0.8, xlog=True, ylog=True, save=False):

	ComPairSensitivity=plotSourceSensitivity(data,angleSelection=0.8,doplot=False)

	a=ascii.read("digitized_alex_sensitivities.dat",names=['eng','sensit'])
	l=ascii.read("differential_flux_sensitivity_p8r2_source_v6_all_10yr_zmax100_n10.0_e1.50_ts25_000_090.txt",names=['emin','emax','e2diff','tmp'])

	#print a['energy']
	energy=a['eng']
	sens=a['sensit']

	erg2mev=624151.
	lateng=(l["emin"]+l["emax"])/2.

	plot.clf()
#	plot.figure()
	#LAT
	plot.plot(lateng,l["e2diff"]*erg2mev,color='magenta',lw=2)
	plot.gca().set_xscale('log')
	plot.gca().set_yscale('log')
	plot.gca().set_xlim([1e-2,1e6])
	plot.gca().set_ylim([1e-8,1e-2])
	plot.gca().set_xlabel('Energy (MeV)')
	plot.gca().set_ylabel(r'Sensitivity $\times\ E^2$ [$\gamma$ MeV s$^{-1}$ cm$^{-2}$]')
	plot.annotate('Fermi-LAT', xy=(5e2,2e-6),xycoords='data',fontsize=14,color='magenta')

	#EGRET
	ind=numpy.arange(69,74,1)
	plot.plot(energy[ind],sens[ind],color='blue',lw=2)
	plot.annotate('EGRET', xy=(1e2,1e-4),xycoords='data',fontsize=14,color='blue')

	#SPI
	ind=numpy.arange(20,46)
	plot.plot(energy[ind],sens[ind],color='green',lw=2)
	plot.annotate('SPI', xy=(6e-2,1e-4),xycoords='data',fontsize=14,color='green')

	#COMPTEL
	comptel_energy=[0.73295844,0.8483429,1.617075,5.057877,16.895761,29.717747]
	comptel_sens=[6.566103E-4,3.6115389E-4,1.4393721E-4,1.6548172E-4,2.36875E-4,3.390693E-4]
	plot.plot(comptel_energy,comptel_sens,color='orange',lw=2)
	plot.annotate('COMPTEL', xy=(5,5e-4),xycoords='data',fontsize=14,color='orange')

	#NuSTAR
	ind=numpy.arange(84,147)
	plot.plot(energy[ind]*1e-3,sens[ind]*(energy[ind]/1e3)**2*1e3,color='purple',lw=2)
	plot.annotate('NuSTAR', xy=(0.1,3e-8),xycoords='data',fontsize=14,color='purple')

	#ComPair
	compair_eng=numpy.array([0.316,1,3.16,10,31.6,100,316.])
	tracked=ComPairSensitivity[1]#numpy.array([  1.54877155e-05,   4.84546681e-06,   5.28735667e-06, 6.53265846e-05, 0, 0, 0])
	untracked=ComPairSensitivity[2]#numpy.array([  2.49626245e-06,   1.82264874e-06,   1.54100276e-05, 9.59603201e-05, 0, 0, 0])
	pair=ComPairSensitivity[3]#numpy.array([ 0, 0,   5.62236032e-05, 3.19254897e-05,   1.71183233e-05,   1.61203804e-05, 2.19339e-05])
	pair_ideal=ComPairSensitivity[4]
	pair_idealang=ComPairSensitivity[5]
	w1=tracked > 0
	plot.plot(compair_eng[w1],tracked[w1],color='black',lw=3)
	w2=pair > 0
	plot.plot(compair_eng[w2],pair[w2],'r--',color='black',lw=3)
	w3=pair_ideal > 0
	plot.plot(compair_eng[w3],pair_ideal[w3],'r:',color='black',lw=3)
	plot.plot(compair_eng[2:],pair_idealang[2:],'r-.',color='black',lw=3)
	plot.annotate('ComPair', xy=(1,1e-6),xycoords='data',fontsize=20)

	if save == True:
		plot.savefig('new_sensitivity.png')
	plot.show()

	return

if __name__=='__main__':


	print dir_path
	if len(sys.argv) == 1:
		print 'Usage: FigureOfMeritParser.py filename'

	else:

		# Get the user supplied simulation output
		simulationFilename = sys.argv[1]

		results = run(simulationFilename)

		print ""
		print results


	







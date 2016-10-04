#!/usr/bin/env python
"""
------------------------------------------------------------------------

Script to parse and visualize results created by mimrec:

Author: Daniel Kocevski (dankocevski@gmail.com)
Date: August 25st, 2016


Usage Examples:

import FigureOfMeritPlotter

# Get the number of triggered and simulated events for a set of simulation files


# Import Mimrec data
simulationIDs = '../Simulations/PerformancePlotSourceFiles/FarFieldPointSourceIDs.txt'
data = FigureOfMeritPlotter.parseMimrecLogs(simulationIDs=simulationIDs)

# Import EventAnalysis data (Check the EventAnalysis doc string to see how to generate the TriggerEfficiency.txt file)
directory = '../Simulations/MyCustomSimulationRun'
triggerEfficiencyFilename = "../Simulations/MyCustomSimulationRun/TriggerEfficiency.txt"
data = FigureOfMeritPlotter.parseEventAnalysisLogs(directory, triggerEfficiencyFilename=triggerEfficiencyFilename)

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
from scipy import interpolate


##########################################################################################

def omega(PSF):

	psf=numpy.array(PSF).astype(float)
	omega_solidAngle = 2*math.pi*(1-numpy.cos(2*psf*math.pi/180.))

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
		if 'ID ' in line:
			numberOfSimulatedEvents = float(line.split()[2])

	fileinput.close()

	return numberOfSimulatedEvents


##########################################################################################

def getMimrecValues(filename, interactionType, verbose=False):

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

def parseMimrecLogs(sumlationsIDs=None):

	# Define the path to the python repository
	pythonPath = os.path.dirname(os.path.realpath(__file__))

	# Define the root path containg the ComPair repositories
	rootPath = pythonPath.replace('python','')

	# Define the PerformancePlotSourceFiles path
	sourceDirectory = rootPath + 'Simulations/PerformancePlotSourceFiles'

	# Define the PerformancePlotMimrecFiles path
	mimrecDirectory = rootPath + 'Simulations/PerformancePlotMimrecFiles'

	# Define the PerformancePlotTraFiles path
	traDirectory = rootPath + 'Simulations/PerformancePlotTraFiles' 

	if sumlationsIDs == None:
		sumlationsIDs = sourceDirectory + '/FarFieldPointSourceIDs.txt'

	# Create a dictionary to store the data
	data = {}


	# Read the sim file
	with open(sumlationsIDs) as filehandle:
		lines = filehandle.readlines()

	currentNumber = 1
	print 'parsing mimrec logs...'

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
		Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents = getMimrecValues(mimrecFilename_tracked, interactionType='Compton')
		data[simulationName].append([Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents])

		# Get the mimrec figures of merit and add the dictionary
		# print "Parsing: %s" % mimrecFilename_untracked	
		Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents = getMimrecValues(mimrecFilename_untracked, interactionType='Compton')
		data[simulationName].append([Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents])

		# Get the mimrec figures of merit and add the dictionary
		# print "Parsing: %s" % mimrecFsilename_pair		
		Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, Containment68, NumberOfReconstructedPairEvents = getMimrecValues(mimrecFilename_pair, interactionType='Pair')
		data[simulationName].append([Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, Containment68, NumberOfReconstructedPairEvents])

		currentNumber = currentNumber + 1


	print 'Done.'
	return data



##########################################################################################

def parseEventAnalysisLogs(directory, triggerEfficiencyFilename=None):

	if triggerEfficiencyFilename == None:
		triggerEfficiencyFilename = directory + '/TriggerEfficiency.txt'

	# Read the sim file
	with open(triggerEfficiencyFilename) as filehandle:
		triggerEfficiencyLines = filehandle.readlines()

	# Create a dictionary to store the data
	data = {}

	for triggerEfficiencyLine in triggerEfficiencyLines:
		lineContents = triggerEfficiencyLine.split()

		# Get the filename
		simulationName = lineContents[0].split('/')[-1]	

		# Create a key for this simulation name
		data[simulationName] = []

		# Get the number of triggered and simulated events
		numberOfTriggeredEvents = lineContents[1]		
		numberOfSimulatedEvents = lineContents[2]

		# Add the number of simulated to the results dictionary
		data[simulationName].append(numberOfSimulatedEvents)

		# Generate the log filename
		analysisLog = directory + '/' + simulationName.replace('.sim', '.log')

		try:
			print "Parsing: %s" % analysisLog
			with open(analysisLog) as filehandle:
				analysisLogLines = filehandle.readlines()

		except:
			print "*** Could not find log file: %s" % analysisLog

		# Loop through the analysis log file
		for analysisLogLine in analysisLogLines:

			if "Results for simulation:" in analysisLogLine:
				analysisLogLine = analysisLogLine.replace("Results for simulation: ",'')
				energy, energySearchUnit, cos, angle, filename = analysisLogLine.split()

			if "Compton Events Reconstructed: " in analysisLogLine:
				numberOfComptonEvents = analysisLogLine.split()[-1]

			if "Compton Energy Resolution (keV): " in analysisLogLine:
				FWHM_energyComptonEvents = analysisLogLine.split()[-1]

			if "Compton Angular Resolution (deg): " in analysisLogLine:
				FWHM_angleComptonEvents = analysisLogLine.split()[-1]

			if "Untracked Compton Events Reconstructed: " in analysisLogLine:
				numberOfUntrackedComptonEvents = analysisLogLine.split()[-1]

			if "Untracked Compton Energy Resolution (keV): " in analysisLogLine:
				FWHM_energyUntrackedComptonEvents = analysisLogLine.split()[-1]

			if "Untracked Compton Angular Resolution (deg): " in analysisLogLine:
				FWHM_angleUntrackedComptonEvents = analysisLogLine.split()[-1]

			if "Tracked Compton Events Reconstructed: " in analysisLogLine:
				numberOfTrackedComptonEvents = analysisLogLine.split()[-1]

			if "Tracked Compton Energy Resolution (keV): " in analysisLogLine:
				FWHM_energyTrackedComptonEvents = analysisLogLine.split()[-1]

			if "Tracked Compton Angular Resolution (deg): " in analysisLogLine:
				FWHM_angleTrackedComptonEvents = analysisLogLine.split()[-1]

			if "Pair Events Reconstructed: " in analysisLogLine:
				numberOfPairEvents = analysisLogLine.split()[-1]

			if "Pair Energy Resolution (keV): " in analysisLogLine:
				FWHM_pairEvents = analysisLogLine.split()[-1]

			if "Pair Angular Containment " in analysisLogLine:
				contaimentData_68 = analysisLogLine.split()[-1]

		# Add all the values to the results dictionary
		holder=-999. # This is required for historic reasons

		if float(energy) <=10.0:
			data[simulationName].append([holder, holder, holder, holder, FWHM_energyComptonEvents, holder, holder, FWHM_angleComptonEvents, numberOfComptonEvents])
			data[simulationName].append([holder, holder, holder, holder, FWHM_energyTrackedComptonEvents, holder, holder, FWHM_angleTrackedComptonEvents, numberOfTrackedComptonEvents])
			data[simulationName].append([holder, holder, holder, holder, FWHM_energyUntrackedComptonEvents, holder, holder, FWHM_angleUntrackedComptonEvents, numberOfUntrackedComptonEvents])
		else:
			data[simulationName].append([numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan])
			data[simulationName].append([numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan])
			data[simulationName].append([numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan])

		if float(energy) >=3.0:
			data[simulationName].append([holder, holder, holder, holder, FWHM_pairEvents, holder, contaimentData_68,numberOfPairEvents])
		else:
			data[simulationName].append([numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan])

	return data		


##########################################################################################

def plotAngularResolution(data, angleSelections=[1,0.9,0.8,0.7,0.6,0.5], xlog=True, ylog=False, save=False):
	
	if hasattr(angleSelections, '__iter__') == False:
		angleSelections = [angleSelections]


	plotNumber = 1
	if len(angleSelections)==1:
		plot.figure(figsize=(8,6))
	else:
		plot.figure(figsize=(10,12))

	for angleSelection in angleSelections:

		Energy = []
		FWHM_tracked = []
		FWHM_untracked = []
		Containment68 = []

		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			#angle = float(key.split('_')[2].replace('Cos',''))
			half = key.split('_')[2].replace('Cos','')
			angle = float(half.replace('.inc1.id1.sim',''))

			if angle == angleSelection:
				Energy.append(energy)
				fwhm_tracked = data[key][2][7]
				fwhm_untracked = data[key][3][7]
				containment68 = data[key][4][6]

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

	if save:
		plot.savefig('AngularResolution.png', bbox_inches='tight')

	plot.show()

	plot.close()

##########################################################################################

def plotAngularResolutionVsAngle(data, energySelections=None, xlog=False, ylog=False, save=False):

	plotNumber = 1

	if energySelections is None:
		if len(energySelections)==1:
			plot.figure(figsize=(8,6))
		else:
			plot.figure(figsize=(10,12))
		Energy = []	
		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			if energy not in Energy:
					Energy.append(energy)
		Energy=numpy.array(Energy)
		i = [numpy.argsort(Energy)]
		energySelections = Energy[i]
	else:	
		if type(energySelections) == float or type(energySelections) == int:
			energySelections=[energySelections]
		energySelections.sort(key=int)
		energySelections=numpy.array(energySelections,dtype=float)

	if len(energySelections)>6:
		energySelections=energySelections[[1,3,5,7,9,11]]
		print "plotting only every other energy: ", energySelections


	for energySelection in energySelections:

		Angle = []
		FWHM_tracked = []
		FWHM_untracked = []
		Containment68 = []

		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			#angle = float(key.split('_')[2].replace('Cos',''))
			half = key.split('_')[2].replace('Cos','')
			angle = numpy.array(float(half.replace('.inc1.id1.sim','')))
			angle = round(numpy.degrees(numpy.arccos(angle)))

			if energy == energySelection:
				Angle.append(angle)
				fwhm_tracked = data[key][2][7]
				fwhm_untracked = data[key][3][7]
				containment68 = data[key][4][6]

				FWHM_tracked.append(fwhm_tracked)
				FWHM_untracked.append(fwhm_untracked)
				Containment68.append(containment68)

		# Convert everything to a numpy array
		Angle = numpy.array(Angle)
		FWHM_tracked = numpy.array(FWHM_tracked)
		FWHM_untracked = numpy.array(FWHM_untracked)
		Containment68 = numpy.array(Containment68)

		# Sort by Angle
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

		if xlog:
			plot.xscale('log')

		if ylog:
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

	if save:
		plot.savefig('AngularResolutionVsAngle.png', bbox_inches='tight')

	plot.show()

	plot.close()


##########################################################################################

def plotEnergyResolution(data, angleSelections=[1,0.9,0.8,0.7,0.6,0.5], xlog=True, ylog=False, save=False):
	
	if hasattr(angleSelections, '__iter__') == False:
		angleSelections = [angleSelections]

	plotNumber = 1
	if len(angleSelections)==1:
		plot.figure(figsize=(8,6))
	else:
		plot.figure(figsize=(10,12))

	for angleSelection in angleSelections:

		Energy = []
		Sigma_tracked = []
		Sigma_untracked = []
		Sigma_pair = []

		# print 'Name, fwhm_tracked, fwhm_untracked, containment'

		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			#angle = float(key.split('_')[2].replace('Cos',''))
			half = key.split('_')[2].replace('Cos','')
			angle = numpy.array(float(half.replace('.inc1.id1.sim','')))

			if angle == angleSelection:
				Energy.append(energy)

				Sigma_tracked.append(data[key][2][4])
				Sigma_untracked.append(data[key][3][4])
				Sigma_pair.append(data[key][4][4])


		# Convert everything to a numpy array
		Energy = numpy.array(Energy)
		Sigma_tracked = numpy.array(Sigma_tracked)
		Sigma_untracked = numpy.array(Sigma_untracked)
		Sigma_pair = numpy.array(Sigma_pair)

		# Sort by energy
		i = [numpy.argsort(Energy)]
		Energy = Energy[i]
		Sigma_tracked = Sigma_tracked[i]
		Sigma_untracked = Sigma_untracked[i]
		Sigma_pair = Sigma_pair[i]

		# Plot the data
		ax = plot.subplot( str(len(angleSelections)) + str(10 + plotNumber) )

		i=Sigma_tracked != 'nan'
		st=numpy.double(Sigma_tracked[i])/numpy.double(Energy[i])*1e-3
		plot.scatter(Energy[i],st,color='darkgreen')
		plot.plot(Energy[i], st, color='darkgreen', alpha=0.5, label='Compton (tracked)')

		i=Sigma_untracked != 'nan'
		sut=numpy.double(Sigma_untracked[i])/numpy.double(Energy[i])*1e-3
		plot.scatter(Energy[i],sut,color='blue')
		plot.plot(Energy[i], sut, color='blue', alpha=0.5, label='Compton (untracked)')

		i=Sigma_pair != 'nan'
		sp=numpy.double(Sigma_pair[i])/numpy.double(Energy[i])*1e-3
		plot.scatter(Energy[i],sp,color='darkred')
		plot.plot(Energy[i],sp, color='darkred', alpha=0.5, label='Pair')

		if plotNumber == 1:
			plot.title('Energy Resolution')			
			plot.legend(numpoints=1, scatterpoints=1, fontsize='small', frameon=True, loc='upper left')

		if xlog:
			plot.xscale('log')

		if ylog:
			plot.yscale('log')

		plot.ylabel(r'$\sigma$ / Energy')

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

	if save:
		plot.savefig('EnergyResolution.png', bbox_inches='tight')

	plot.show()

	plot.close()


##########################################################################################

def plotEnergyResolutionVsAngle(data, energySelections=None, xlog=False, ylog=False, save=False):

	plotNumber = 1
	plot.figure(figsize=(10,12))

	if energySelections is None:
		Energy = []	
		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			if energy not in Energy:
					Energy.append(energy)
		Energy=numpy.array(Energy)
		i = [numpy.argsort(Energy)]
		energySelections = Energy[i]
	else:	
		if type(energySelections) == float or type(energySelections) == int:
			energySelections=[energySelections]
		energySelections.sort(key=int)
		energySelections=numpy.array(energySelections,dtype=float)

	if len(energySelections)>6:
		print "plotting only every other energy"
		energySelections=energySelections[[1,3,5,7,9,11]]

	for energySelection in energySelections:

		Angle = []
		Sigma_tracked = []
		Sigma_untracked = []
		Sigma_pair = []

		# print 'Name, fwhm_tracked, fwhm_untracked, containment'

		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			#angle = float(key.split('_')[2].replace('Cos',''))
			half = key.split('_')[2].replace('Cos','')
			angle = numpy.array(float(half.replace('.inc1.id1.sim','')))

			if energy == energySelection:
				Angle.append(angle)

				Sigma_tracked.append(data[key][2][4])
				Sigma_untracked.append(data[key][3][4])
				Sigma_pair.append(data[key][4][4])

		# Convert everything to a numpy array
		Angle = numpy.degrees(numpy.arccos(Angle))
		Angle = numpy.array(Angle)
		Sigma_tracked = numpy.array(Sigma_tracked)
		Sigma_untracked = numpy.array(Sigma_untracked)
		Sigma_pair = numpy.array(Sigma_pair)

		# Sort by energy
		i = [numpy.argsort(Angle)]
		Angle = Angle[i]
		Sigma_tracked = Sigma_tracked[i]
		Sigma_untracked = Sigma_untracked[i]
		Sigma_pair = Sigma_pair[i]

		# Plot the data
		ax = plot.subplot( str(len(energySelections)) + str(10 + plotNumber) )
		#ax = plot.subplots(len(energySelections)/2,2)
		#ax = plot.subplot(2,len(energySelections)/2+1,len(energySelections))

		
		i=Sigma_tracked =='nan'
		st=numpy.double(numpy.ma.array(Sigma_tracked,mask=i).compressed())/numpy.double(energySelection)*1.0e-3
		plot.plot(numpy.ma.array(Angle,mask=i).compressed(), st, color='darkgreen', alpha=0.75, label='Compton (tracked)', marker='o')
		
		i=Sigma_untracked == 'nan'
		sut=numpy.double(numpy.ma.array(Sigma_untracked,mask=i).compressed())/numpy.double(energySelection)*1e-3		
		plot.plot(numpy.ma.array(Angle,mask=i).compressed(), sut, color='blue', alpha=0.75, label='Compton (untracked)', marker='o')
		
		i=Sigma_pair == 'nan'
		sp=numpy.double(numpy.ma.array(Sigma_pair,mask=i).compressed())/numpy.double(energySelection)*1e-3		
		plot.plot(numpy.ma.array(Angle,mask=i).compressed(), sp, color='darkred', alpha=0.75, label='Pair', marker='o')

		if plotNumber == 1:
			plot.title('Energy Resolution')						
			plot.legend(numpoints=1, scatterpoints=1, fontsize='small', frameon=True, loc='upper left')

		if xlog:
			plot.xscale('log')

		if ylog:
			plot.yscale('log')

		plot.ylabel(r'$\sigma$ / Energy')
		plot.xlim([0,60])

		plot.text(1-0.015, 0.8, '%s MeV' % energySelection,
		        verticalalignment='top', horizontalalignment='right',
		        transform=ax.transAxes,
		        color='black', fontsize=12)

		#if plotNumber != len(energySelections):
		#	ax.set_xticklabels([])

		if plotNumber == len(energySelections):
			plot.xlabel(r'$\theta$ (deg)')

		plotNumber = plotNumber + 1


	plot.subplots_adjust(wspace=0, hspace=.2)

	if save:
		plot.savefig('EnergyResolutionVsAngle.png', bbox_inches='tight')

	plot.show()

	plot.close()

##########################################################################################

def plotEffectiveArea(data, angleSelections=[1,0.9,0.8,0.7,0.6,0.5], ideal=False, xlog=True, ylog=False, save=False):

	if hasattr(angleSelections, '__iter__') == False:
		angleSelections = [angleSelections]

	plotNumber = 1
	if len(angleSelections)==1:
		plot.figure(figsize=(8,6))
	else:
		plot.figure(figsize=(10,12))

	for angleSelection in angleSelections:

		Energy = []
		EffectiveArea_Tracked = []
		EffectiveArea_Untracked  = []
		EffectiveArea_Pair  = []

		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			#angle = float(key.split('_')[2].replace('Cos',''))
			half = key.split('_')[2].replace('Cos','')
			angle = float(half.replace('.inc1.id1.sim',''))
			#print energy, angle

			if angle == angleSelection:
				#print data
				numberOfSimulatedEvents = float(data[key][0])

				if ideal:
					#This removes the event selection on the final Aeff calculation
					#It does not change anything from the FWHM or the 68% containment
					#Compton events are multiplied by the ratio of tracked vs. untracked
					total_compton_events = float(data[key][1][-1])
					numberOfReconstructedEvents_tracked = 100000.*float(data[key][2][-1])/(total_compton_events)
					numberOfReconstructedEvents_untracked = 100000.*float(data[key][3][-1])/(total_compton_events)
					numberOfReconstructedEvents_pair = 100000.
				else:
					numberOfReconstructedEvents_tracked = float(data[key][2][-1])
					numberOfReconstructedEvents_untracked = float(data[key][3][-1])
					numberOfReconstructedEvents_pair = float(data[key][4][-1])

				#numberOfReconstructedEvents_tracked = float(data[key][1][-1])
				#numberOfReconstructedEvents_untracked = float(data[key][2][-1])
				#numberOfReconstructedEvents_pair = float(data[key][3][-1])

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

		if ideal:
			EffectiveArea_Pair[0]=numpy.nan
			EffectiveArea_Pair[1]=numpy.nan


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

		#plot.gca().set_ylim([0.001,10000.])

		if plotNumber == 1:
			plot.title('Effective Area')			

		if xlog:
			plot.xscale('log')

		if ylog:
			plot.yscale('log')


		if plotNumber != len(angleSelections):
			ax.set_xticklabels([])

		if plotNumber == len(angleSelections):
			plot.xlabel(r'Energy (MeV)')

		plotNumber = plotNumber + 1


	plot.subplots_adjust(wspace=0, hspace=.2)

	if save:
		plot.savefig('EffectiveArea.png', bbox_inches='tight')

	plot.show()


 ##########################################################################################

def plotEffectiveAreaVsAngle(data, energySelections=None, ideal=False, xlog=False, ylog=False, save=False, collapse=False):

	if energySelections is None:
		Energy = []	
		for key in data.keys():
			energy = float(key.split('_')[1].replace('MeV',''))
			if energy not in Energy:
					Energy.append(energy)
		Energy=numpy.array(Energy)
		i = [numpy.argsort(Energy)]
		energySelections = Energy[i]
	else:	
		if type(energySelections) == float or type(energySelections) == int:
			energySelections=[energySelections]
		energySelections.sort(key=int)
		energySelections=numpy.array(energySelections,dtype=float)

	if len(energySelections)>6:
		print "plotting only every other energy"
		energySelections=energySelections[[1,3,5,7,9,11]]

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
			#angle = float(key.split('_')[2].replace('Cos',''))
			half = key.split('_')[2].replace('Cos','')
			angle = numpy.array(float(half.replace('.inc1.id1.sim','')))
			angle = round(numpy.degrees(numpy.arccos(angle)))

			if energy == energySelection:

				if ideal:
					# This removes the event selection on the final Aeff calculation
					# It does not change anything from the FWHM or the 68% containment
					# Compton events are multiplied by the ratio of tracked vs. untracked
					total_compton_events = float(data[key][1][-1])
					numberOfReconstructedEvents_tracked = 100000.*float(data[key][2][-1])/(total_compton_events)
					numberOfReconstructedEvents_untracked = 100000.*float(data[key][3][-1])/(total_compton_events)
					numberOfReconstructedEvents_pair = 100000.
				else:
					numberOfReconstructedEvents_tracked = float(data[key][2][-1])
					numberOfReconstructedEvents_untracked = float(data[key][3][-1])
					numberOfReconstructedEvents_pair = float(data[key][4][-1])

				numberOfSimulatedEvents = float(data[key][0])

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

		if xlog:
			plot.xscale('log')

		if ylog:
			plot.yscale('log')


		if plotNumber != len(energySelections) and collapse == False:
			ax.set_xticklabels([])

		if plotNumber == len(energySelections) or collapse == True:
			plot.xlabel(r'$\theta$')

		plotNumber = plotNumber + 1


	plot.subplots_adjust(wspace=0, hspace=.2)

	if save:
		plot.savefig('EffectiveAreaVsAngle.png', bbox_inches='tight')

	plot.show()


########################################################################################## 

def plotSourceSensitivity(data, angleSelection=0.8, exposure = 6.3*10**6, ideal=False, doPSF=None, xlog=True, ylog=True, save=False, doplot=False):

	#background = numpy.array([0.00346008, 0.00447618, 0.00594937, 0.00812853, 0.0100297, 0.0124697, 0.0161290])
	#background=numpy.array([0.00346008,0.00378121,0.00447618,0.00504666,0.00594937,0.00712394,0.00812853,0.00881078,0.0100297,0.0109190,0.0124697,0.0139781,0.0161290])

	Energy = []
	EffectiveArea_Tracked = []
	EffectiveArea_Untracked  = []
	EffectiveArea_Pair  = []
	FWHM_tracked = []
	FWHM_untracked = []
	Containment68 = []

	for key in data.keys():

		energy = float(key.split('_')[1].replace('MeV',''))

		#angle = float(key.split('_')[2].replace('Cos',''))
		half = key.split('_')[2].replace('Cos','')
		angle = float(half.replace('.inc1.id1.sim',''))
		if angle == angleSelection:

			# Get the number of reconstructed events
			numberOfSimulatedEvents = float(data[key][0])

			if ideal:
				#This removes the event selection on the final Aeff calculation
				#It does not change anything from the FWHM or the 68% containment
				#Compton events are multiplied by the ratio of tracked vs. untracked
				total_compton_events = float(data[key][1][-1])
				numberOfReconstructedEvents_tracked = 100000.*float(data[key][2][-1])/(total_compton_events)
				numberOfReconstructedEvents_untracked = 100000.*float(data[key][3][-1])/(total_compton_events)
				numberOfReconstructedEvents_pair = 100000.
			else:
				numberOfReconstructedEvents_tracked = float(data[key][2][-1])
				numberOfReconstructedEvents_untracked = float(data[key][3][-1])
				numberOfReconstructedEvents_pair = float(data[key][4][-1])
            
			# Get the angular resolution
			fwhm_tracked = data[key][2][7]
			fwhm_untracked = data[key][3][7]
			containment68 = data[key][4][6]

			# Calculate the effective area
			effectiveArea_tracked = (numberOfReconstructedEvents_tracked/numberOfSimulatedEvents) * math.pi * 300**2
			effectiveArea_untracked = (numberOfReconstructedEvents_untracked/numberOfSimulatedEvents) * math.pi * 300**2
			effectiveArea_pair = (numberOfReconstructedEvents_pair/numberOfSimulatedEvents) * math.pi * 300**2
		
			# Store the effective area results
			Energy.append(energy)
			EffectiveArea_Tracked.append(effectiveArea_tracked)
			EffectiveArea_Untracked.append(effectiveArea_untracked)
			EffectiveArea_Pair.append(effectiveArea_pair)
		
			# Store the angular resolution results
			FWHM_tracked.append(fwhm_tracked)
			FWHM_untracked.append(fwhm_untracked)
			Containment68.append(containment68)

	# Convert everything to a numpy array
	Energy = numpy.array(Energy)
	EffectiveArea_Tracked = numpy.array(EffectiveArea_Tracked)
	EffectiveArea_Untracked = numpy.array(EffectiveArea_Untracked)
	EffectiveArea_Pair = numpy.array(EffectiveArea_Pair)

	FWHM_tracked = numpy.array(FWHM_tracked)
	FWHM_untracked = numpy.array(FWHM_untracked)
	Containment68 = numpy.array(Containment68)

	# Sort by energy
	i = [numpy.argsort(Energy)]
	Energy = Energy[i]
	EffectiveArea_Tracked = EffectiveArea_Tracked[i]
	EffectiveArea_Untracked = EffectiveArea_Untracked[i]
	EffectiveArea_Pair = EffectiveArea_Pair[i]

	FWHM_tracked = FWHM_tracked[i]
	FWHM_untracked = FWHM_untracked[i]
	Containment68 = Containment68[i]
	if doPSF is not None:
		#Note: this only changes the PSF for Pair events, not Compton!
		Containment68 = doPSF	

	# digitized from Stong, Moskalenko & Reimer 2000, Figure 8, top right
  	# high latitude |b|>5 deg
  	# multiply by 10 to vaguely account for the albedo background

	eng2=numpy.array([0.10355561,0.3534914,1.2920963,4.659387,8.969312,18.735151,38.081676,69.40132,144.98259,227.4451,342.42523,462.24567,725.01324,939.413,1908.1061,28725.793])
	e2int2=numpy.array([2.7943178E-4,3.57757E-4,4.8821748E-4,6.806025E-4,8.0072926E-4,9.1560354E-4,0.0010469892,0.0011638523,0.0013691497,0.0015439879,0.0016334692,0.0017039803,0.0018284274,0.0018672496,0.0017879958,0.0014717471])

	# interpolate background at our energies
	tck=interpolate.splrep(eng2,e2int2,s=0)
	background=interpolate.splev(Energy,tck,der=0)

	Sensitivity_tracked = Isrc(Energy, exposure, EffectiveArea_Tracked, 3., omega(FWHM_tracked), 10*background)
	Sensitivity_untracked = Isrc(Energy, exposure, EffectiveArea_Untracked, 3., omega(FWHM_untracked), 10*background)
	Sensitivity_pair = Isrc(Energy, exposure, EffectiveArea_Pair, 3., omega(Containment68), background)

	if doPSF:
		Sensitivity_pair[0]=numpy.nan
		Sensitivity_pair[1]=numpy.nan
	
	if plot.fignum_exists(1):
		plot.clf()
	else:
		plot.figure(figsize=(10, 6.39))
		ax = plot.subplot(111)

	plot.scatter(Energy, Sensitivity_tracked, color='darkgreen')
	plot.plot(Energy, Sensitivity_tracked, color='darkgreen', alpha=0.5, label='Compton (tracked)')
	plot.scatter(Energy, Sensitivity_untracked, color='blue')	
	plot.plot(Energy, Sensitivity_untracked, color='blue', alpha=0.5, label='Compton (untracked)')
	plot.scatter(Energy, Sensitivity_pair, color='darkred')
	plot.plot(Energy, Sensitivity_pair, color='darkred', alpha=0.5, label='Pair')	

	if xlog:
		plot.xscale('log')

	if ylog:
		plot.yscale('log')

	plot.legend(numpoints=1, scatterpoints=1, fontsize='small', frameon=True, loc='upper left')

	if save:
		plot.savefig('SourceSensitivity.png', bbox_inches='tight')

	if doplot:
		plot.show()

	return Energy, Sensitivity_tracked, Sensitivity_untracked, Sensitivity_pair 


##########################################################################################

def plotAllSourceSensitivities(data, angleSelection=0.8, plotIdeal=False, xlog=True, ylog=True, save=False):

	ComPairSensitivity=plotSourceSensitivity(data,angleSelection=angleSelection,doplot=False)
	ComPairIdealSensitivity=plotSourceSensitivity(data,angleSelection=angleSelection,ideal=True,doplot=False)
	ComPairGoodPSFSensitivity=plotSourceSensitivity(data,angleSelection=angleSelection,ideal=True,doPSF=2.0,doplot=False)

	a=ascii.read("digitized_alex_sensitivities.dat",names=['eng','sensit'])
	l=ascii.read("differential_flux_sensitivity_p8r2_source_v6_all_10yr_zmax100_n10.0_e1.50_ts25_000_090.txt",names=['emin','emax','e2diff','tmp'])

	#print a['energy']
	energy=a['eng']
	sens=a['sensit']

	erg2mev=624151.
	lateng=(l["emin"]+l["emax"])/2.

	plot.clf()
	#LAT
	plot.plot(lateng,l["e2diff"]*erg2mev,color='magenta',lw=2)
	plot.gca().set_xscale('log')
	plot.gca().set_yscale('log')
	plot.gca().set_xlim([1e-2,1e6])
	plot.gca().set_ylim([1e-8,1e-2])
	plot.gca().set_xlabel('Energy (MeV)')
	plot.gca().set_ylabel(r'$3\sigma$ Continuum Sensitivity $\times\ E^2$ [$\gamma$ MeV s$^{-1}$ cm$^{-2}$]')
	plot.annotate('Fermi-LAT', xy=(1e3,1e-6),xycoords='data',fontsize=16,color='magenta')

	#EGRET from https://heasarc.gsfc.nasa.gov/docs/cgro/egret/egret_tech.html#N_4_
	ind=numpy.arange(69,74,1)
	#plot.plot(energy[ind],sens[ind],color='blue',lw=2)
	egret_energy=numpy.array([35,100,200,500,3000,10000.])
	#egret_aeff=numpy.array([0.3,1.1,1.5,1.6,1.1,0.7])*1e3
	egret_aeff=numpy.array([0.07,0.9,1.4,1.5,1.1,0.7])*1e3
	egret_psf_eng_low=numpy.array([35,70,150,500,2000])
	egret_psf_eng_high=numpy.array([70,150,500,2000,30000])
	egret_psf_eng_mid=(egret_psf_eng_high+egret_psf_eng_low)/2.
	egret_psf_4eng=numpy.array([4.3,2.6,1.4,0.7,0.4])
	tck=interpolate.splrep(egret_psf_eng_mid,egret_psf_4eng,s=0)
	egret_psf=interpolate.splev(egret_energy,tck,der=0)
	egret_back=7.32e-9*(egret_energy/451.)**(-2.1) #from Sreekumar et al. 1998 in photons/cm2/s/MeV
	egret_exposure=86400.*7.*2.*0.4 #seconds in 2 weeks *0.4 efficiency
	egret_inclin=1/0.5 #from sensitivity vs inclination angle = 15 deg
	egret_sensitivity=egret_inclin*Isrc(egret_energy,egret_exposure,egret_aeff,3,omega(egret_psf),egret_back)
	plot.plot(egret_energy,egret_sensitivity,color='blue',lw=2)
	plot.annotate('EGRET', xy=(4e2,1e-4),xycoords='data',fontsize=16,color='blue')

	#SPI
	ind=numpy.arange(20,46)
	plot.plot(energy[ind],sens[ind],color='green',lw=2)
	plot.annotate('SPI', xy=(6e-2,1e-4),xycoords='data',fontsize=16,color='green')

	#COMPTEL
	comptel_energy=[0.73295844,0.8483429,1.617075,5.057877,16.895761,29.717747]
	comptel_sens=[6.566103E-4,3.6115389E-4,1.4393721E-4,1.6548172E-4,2.36875E-4,3.390693E-4]
	plot.plot(comptel_energy,comptel_sens,color='orange',lw=2)
	# from Schonfelder et al. (1993)
	#comptel_psf=1.247/(1-exp(-0.854*energy**0.9396))
	#comptel_sensitivity=Isrc(energy,comptel_exposure,comptel_aeff,3,omega(comptel_psf),comptel_back)
	comptel_energy_low=numpy.array([0.7,1,3,10])
	comptel_energy_high=numpy.array([1,3,10,30])
	comptel_energy=10**((numpy.log10(comptel_energy_low)+numpy.log10(comptel_energy_high))/2.)
	comptel_aeff=44.95*(numpy.log(comptel_energy)+0.853)**1.655*comptel_energy**(-0.8593)
#	comptel_back=numpy.array([9600,51000,16900,1150.])
#	comptel_exposure=numpy.array([2.3e6,3.9e6,7.1e6,10.8e6])
#	comptel_sensitivity=3*comptel_back**0.5/(comptel_aeff*comptel_exposure)#*comptel_energy**2
	comptel_sensitivity=numpy.array([1.2e-4,1.9e-4,7.9e-5,1.6e-5])*comptel_energy**2/(comptel_energy_high-comptel_energy_low)

	#plot.plot(comptel_energy,comptel_sensitivity,'r--',color='orange',lw=2)
	plot.annotate('COMPTEL', xy=(5,5e-4),xycoords='data',fontsize=16,color='orange')

	#NuSTAR
	ind=numpy.arange(84,147)
	plot.plot(energy[ind]*1e-3,sens[ind]*(energy[ind]/1e3)**2*1e3,color='purple',lw=2)
	plot.annotate('NuSTAR', xy=(0.1,3e-8),xycoords='data',fontsize=16,color='purple')

	#ComPair
	compair_eng=ComPairSensitivity[0]
	tracked=ComPairSensitivity[1]	
	untracked=ComPairSensitivity[2]
	pair=ComPairSensitivity[3]
	if plotIdeal:
		tracked_ideal=ComPairIdealSensitivity[1]
		pair_ideal=ComPairIdealSensitivity[3]
		pair_idealang=ComPairGoodPSFSensitivity[3]

	plot.plot(compair_eng,tracked,color='black',lw=3)	
	plot.plot(compair_eng,pair,'r--',color='black',lw=3)
	if plotIdeal:
		plot.plot(compair_eng,tracked_ideal,'r:',color='black',lw=3)
		plot.plot(compair_eng,pair_ideal,'r:',color='black',lw=3)
		plot.plot(compair_eng,pair_idealang,'r-.',color='black',lw=2)

	#Alex's numbers
	ind=numpy.arange(158,167)
	plot.plot(energy[ind],sens[ind],'r--',color='red',lw=2)
	plot.annotate('Previous ComPair', xy=(7e2,3e-6),xycoords='data',fontsize=14,color='red')

	#plot.plot(compair_eng,pair_idealang,'r-.',color='black',lw=3)
	plot.annotate('ComPair', xy=(1,1e-6),xycoords='data',fontsize=22)

	if save:
		plot.savefig('new_sensitivity.png')
	plot.show()

	return

##########################################################################################

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


	







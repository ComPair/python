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
import matplotlib
import matplotlib.pylab as plot
import math
from astropy.io import ascii
from astropy.table import Table
from scipy import interpolate, nansum

matplotlib.rcParams.update({'font.size': 14})
colors = ['red', 'blue', 'green', 'orange', 'brown', 'purple', 'darkred']


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
	print('parsing mimrec logs...')

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


	print('Done.')
	return data



##########################################################################################

def parseEventAnalysisLogs(directory, triggerEfficiencyFilename=None, silent=False):

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
            if silent==False:
                print("Parsing: %s" % analysisLog)
            with open(analysisLog) as filehandle:
                analysisLogLines = filehandle.readlines()

        except:
            print("*** Could not find log file: %s" % analysisLog)

        # Loop through the analysis log file
        numberOfPairEventsIdeal = None
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

            if "Events Not Reconstructed Flagged as Bad" in analysisLogLine:
                numberOfNotReconstructedEvents = analysisLogLine.split()[-1]

            if "Pair Events Ideal: " in analysisLogLine:
                numberOfPairEventsIdeal = analysisLogLine.split()[-1]
                print("pair events ideal:  ", numberOfPairEventsIdeal)

        # Add all the values to the results dictionary
        
        if numberOfPairEventsIdeal is None:
            numberOfPairEventsIdeal = numpy.nan

        holder=-999. # This is required for historic reasons

        if float(energy) <=30.0:
            data[simulationName].append([holder, holder, holder, holder, FWHM_energyComptonEvents, holder, holder, FWHM_angleComptonEvents, numberOfComptonEvents])
            data[simulationName].append([holder, holder, holder, holder, FWHM_energyTrackedComptonEvents, holder, holder, FWHM_angleTrackedComptonEvents, numberOfTrackedComptonEvents])
            data[simulationName].append([holder, holder, holder, holder, FWHM_energyUntrackedComptonEvents, holder, holder, FWHM_angleUntrackedComptonEvents, numberOfUntrackedComptonEvents])
        else:
            data[simulationName].append([numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan])
            data[simulationName].append([numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan])
            data[simulationName].append([numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan])

    #		if float(energy) >=1.0:
    #			data[simulationName].append([holder, holder, holder, holder, FWHM_pairEvents, holder, contaimentData_68, numberOfPairEvents, numberOfNotReconstructedEvents])

        if float(energy) >=3.0:
                data[simulationName].append([holder, holder, holder, holder, FWHM_pairEvents, numberOfPairEventsIdeal, contaimentData_68, numberOfPairEvents, numberOfNotReconstructedEvents])
        else:
            data[simulationName].append([numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan])

    return data


##########################################################################################

def plotAngularResolution(data, angleSelections=[1,0.9,0.8,0.7,0.6,0.5], xlog=True, ylog=False, 
	 						save=False, collapse=False, doplot=True, txtOutfileLabel='xxx'):
	
    if hasattr(angleSelections, '__iter__') == False:
        angleSelections = [angleSelections]

    plotNumber = 1
    if doplot:
        if len(angleSelections)==1:
            plot.figure(figsize=(8,6))
        elif collapse == False:
            # Create the new subplot
            plot.figure(figsize=(10,12))
        else:
            #print "collapse!"
            plot.figure(figsize=(10, 6.39))
            ax = plot.subplot(111)


    for angleSelection in angleSelections:
        results_txt_TC = open(f"{txtOutfileLabel}_AngRes_Cos{angleSelection}_TC.txt", 'w')
        results_txt_UC = open(f"{txtOutfileLabel}_AngRes_Cos{angleSelection}_UC.txt", 'w')
        results_txt_P =  open(f"{txtOutfileLabel}_AngRes_Cos{angleSelection}_P.txt", 'w')
    	
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
        t = [numpy.argsort(Energy)]
        Energy = Energy[t]
        FWHM_tracked = FWHM_tracked[t]
        FWHM_untracked = FWHM_untracked[t]
        Containment68 = Containment68[t]

        # remove nan's
        i=FWHM_tracked != 'nan'
        st=numpy.double(FWHM_tracked)

        j=FWHM_untracked != 'nan'
        sut=numpy.double(FWHM_untracked)

        k=Containment68 != 'nan'
        sp=numpy.double(Containment68)
        
        
    	# writing txt files	
        results_txt_TC.write("# Energy[MeV] AngRes_TkrCompton[deg]\n")
        for ii, en in enumerate(Energy[i]):
        	results_txt_TC.write("%.1f\t%.1f\n"%(en, st[i][ii]))
        results_txt_TC.close()
        print('Created %s_AngRes_Cos%s_TC.txt ...!'%(txtOutfileLabel, angleSelection))
        	
        results_txt_UC.write("# Energy[MeV] AngRes_UntkrCompton[deg]\n")
        for ii, en in enumerate(Energy[j]):
        	results_txt_UC.write("%.1f\t%.1f\n"%(en, sut[j][ii]))
        results_txt_UC.close()
        print('Created %s_AngRes_Cos%s_UC.txt ...!'%(txtOutfileLabel, angleSelection))
        
        results_txt_P.write("# Energy[MeV] AngRes_Pair[deg]\n")
        for ii, en in enumerate(Energy[k]):
        	results_txt_P.write("%.1f\t%.1f\n"%(en, sp[k][ii]))
        results_txt_P.close()
        print('Created %s_AngRes_Cos%s_P.txt ...!'%(txtOutfileLabel, angleSelection))


        # plot the data
        if doplot:
            if collapse == False:
                ax = plot.subplot( str(len(angleSelections)) + str(10 + plotNumber) )

                plot.scatter(Energy[i],st[i],color='darkgreen')
                plot.plot(Energy[i], st[i], color='darkgreen', alpha=0.5, label='Tracked Compton',lw=2)
                #plot.plot(Energy[i], st[i], color='darkgreen', alpha=0.5, label='Compton', lw=2)

                #plot.scatter(Energy[j],sut,color='blue')
                #plot.plot(Energy[j], sut, color='blue', alpha=0.5, label='Compton (untracked)')
                plot.scatter(Energy[j],sut[j],color='blue')
                plot.plot(Energy[j], sut[j], color='blue', alpha=0.5, label='Untracked Compton', lw=2)

                plot.scatter(Energy[k],sp[k],color='darkred')
                plot.plot(Energy[k],sp[k], color='darkred', alpha=0.5, label='Pair', lw=2)		

                #plot.text(0.015, 0.8, '%i$^\circ$' % round(numpy.degrees(numpy.arccos(angleSelection))),
                 #   	verticalalignment='bottom', horizontalalignment='left',
                 #  		transform=ax.transAxes, color='black', fontsize=16)

            else:
                angle = round(numpy.degrees(numpy.arccos(angleSelection)))
                plot.scatter(Energy[i],st[i], color=colors[plotNumber-1])
                plot.plot(Energy[i], st[i], color=colors[plotNumber-1], alpha=0.5, lw=2, label='Compton at %i$^\circ$' % angle)

                plot.scatter(Energy[j],sut, color=colors[plotNumber-1])
                plot.plot(Energy[j], sut, color=colors[plotNumber-1], alpha=0.5, linestyle='-.', lw=2)

                plot.scatter(Energy[k],sp[k], color=colors[plotNumber-1])
                plot.plot(Energy[k],sp[k], color=colors[plotNumber-1], alpha=0.5, linestyle='--', lw=2, label='Pair at %i$^\circ$' % angle)

            if plotNumber == len(angleSelections):
                #plot.title('Angular Resolution')			
                plot.legend(numpoints=1, scatterpoints=1, fontsize=16, frameon=True, loc='upper right')

            if xlog == True:
                plot.xscale('log')

            if ylog == True:
                plot.yscale('log')

            plot.ylabel('Angular Resolution ($^{\circ}$)')


            if plotNumber != len(angleSelections) and collapse== False:
                ax.set_xticklabels([])


            if plotNumber == len(angleSelections):
                plot.xlabel('Energy (MeV)', fontsize=16)


                plot.gca().set_ylim([0.,20.])

                plot.subplots_adjust(wspace=0, hspace=.2)

                if save:
                    plot.savefig('AngularResolution_Cos%s.pdf' % angleSelections[0])
                    plot.savefig('AngularResolution_Cos%s.png' % angleSelections[0])

                plot.show()

                plot.close()

            plotNumber = plotNumber + 1


    return Energy,st,sp

##########################################################################################

def plotAngularResolutionVsAngle(data, energySelections=None, xlog=False, ylog=False, save=False, collapse=False):

	plotNumber = 1

	if energySelections is None:
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


	if collapse == True:
		plot.figure(figsize=(10, 6.39))
		ax = plot.subplot(111)

	if len(energySelections)>6:
		energySelections=energySelections[[1,3,5,7,9,11]]
		print("plotting only every other energy: ", energySelections)


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
		if collapse==False:
			ax = plot.subplot( str(len(energySelections)) + str(10 + plotNumber) )

			plot.scatter(Angle,FWHM_tracked,color='darkgreen')
			plot.plot(Angle, FWHM_tracked, color='darkgreen', alpha=0.5, label='Compton', lw=2)

			#plot.scatter(Angle,FWHM_untracked,color='blue')
			#plot.plot(Angle, FWHM_untracked, color='blue', alpha=0.5, label='Compton (untracked)', lw=2)

			plot.scatter(Angle,sp,color='darkred')
			plot.plot(Angle,sp, color='darkred', alpha=0.5, label='Pair', lw=2)		

			plot.text(0.015, 0.8, '%s MeV' % energySelection,
		        	verticalalignment='bottom', horizontalalignment='left',
		        	transform=ax.transAxes,
		        	color='black', fontsize=16)

		else:
			if energySelection<3.:
				plot.scatter(Angle,FWHM_tracked, color=colors[plotNumber-1])
				plot.plot(Angle, FWHM_tracked, color=colors[plotNumber-1], alpha=0.5, 
						lw=2, label='Compton at %s  MeV' % energySelection)

			#plot.scatter(Angle,FWHM_untracked, color=colors[plotNumber-1])
			#plot.plot(Angle, FWHM_untracked, color=colors[plotNumber-1], alpha=0.5, lw=2, linestyle='-.')

			#print Containment68
			if energySelection>3.:
				#print Containment68[0], Angle, energySelection
				plot.scatter(Angle,Containment68, color=colors[plotNumber-1])
				plot.plot(Angle,Containment68, color=colors[plotNumber-1], alpha=0.5, label='Pair at %s MeV' % energySelection, 
						lw=2, linestyle='--')


		if plotNumber == len(energySelections):
			#plot.title('Angular Resolution')
			plot.legend(numpoints=1, scatterpoints=1, fontsize=16, frameon=True, loc='upper right')

		if xlog:
			plot.xscale('log')

		if ylog:
			plot.yscale('log')

		plot.ylabel('Angular Resolution ($^{\circ}$)', fontsize=16)

		if plotNumber == len(energySelections):
			plot.xlabel(r'$\theta$', fontsize=16)


		plotNumber = plotNumber + 1

	plot.ylim([1.0,10])

	plot.subplots_adjust(wspace=0, hspace=.2)

	if save:
		plot.savefig('AngularResolutionVsAngle_%sMeV.pdf' % energySelections[0])
		plot.savefig('AngularResolutionVsAngle_%sMeV.png' % energySelections[0])

	plot.show()

	plot.close()

##########################################################################################

def plotEnergyResolution(data, angleSelections=[1,0.9,0.8,0.7,0.6,0.5], xlog=True, ylog=False, 
                                         save=False, collapse=False, txtOutfileLabel='xxx'):

    if hasattr(angleSelections, '__iter__') == False:
        angleSelections = [angleSelections]

    plotNumber = 1
    if len(angleSelections)==1:
        plot.figure(figsize=(8,6))
    elif collapse == False:
        # Create the new subplot
        plot.figure(figsize=(10,12))
    else:
        plot.figure(figsize=(10, 6.39))
        ax = plot.subplot(111)

    for angleSelection in angleSelections:

        results_txt_TC = open( '%s_EnRes_Cos%s_TC.txt' % (txtOutfileLabel, angleSelection), 'w')
        results_txt_UC = open( '%s_EnRes_Cos%s_UC.txt' % (txtOutfileLabel, angleSelection), 'w')
        results_txt_P = open( '%s_EnRes_Cos%s_P.txt' % (txtOutfileLabel, angleSelection), 'w')

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

        # Remove the nan's
        i=Sigma_tracked != 'nan'
        #st=numpy.double(Sigma_tracked[i])/numpy.double(Energy[i])*1e-3
        st=numpy.double(Sigma_tracked[i])/numpy.double(Energy[i])*1e-3*2.355

        j=Sigma_untracked != 'nan'
        #sut=numpy.double(Sigma_untracked[j])/numpy.double(Energy[j])*1e-3
        sut=numpy.double(Sigma_untracked[j])/numpy.double(Energy[j])*1e-3*2.355

        k=Sigma_pair != 'nan'
        #sp=numpy.double(Sigma_pair[k])/numpy.double(Energy[k])*1e-3
        sp=numpy.double(Sigma_pair[k])/numpy.double(Energy[k])*1e-3*2.355


        # writing txt files	
        results_txt_TC.write("# Energy[MeV] EnRes_TkrCompton[FWHM/Energy]\n")
        for ii, en in enumerate(Energy[i]):
        	results_txt_TC.write("%.1f\t%.3f\n"%(en, st[ii]))
        results_txt_TC.close()
        print('Created %s_EnRes_Cos%s_TC.txt ...!'%(txtOutfileLabel, angleSelection))
        	
        results_txt_UC.write("# Energy[MeV] EnRes_UntkrCompton[FWHM/Energy]\n")
        for ii, en in enumerate(Energy[j]):
        	results_txt_UC.write("%.1f\t%.3f\n"%(en, sut[ii]))
        results_txt_UC.close()
        print('Created %s_EnRes_Cos%s_UC.txt ...!'%(txtOutfileLabel, angleSelection))
        
        results_txt_P.write("# Energy[MeV] EnRes_Pair[FWHM/Energy]\n")
        for ii, en in enumerate(Energy[k]):
         	results_txt_P.write("%.1f\t%.3f\n"%(en, sp[ii]))
        results_txt_P.close()
        print('Created %s_EnRes_Cos%s_P.txt ...!'%(txtOutfileLabel, angleSelection))
        

        # Print the data
        if collapse==False:
            ax = plot.subplot( str(len(angleSelections)) + str(10 + plotNumber) )

            plot.scatter(Energy[i],st,color='darkgreen')
            #plot.plot(Energy[i], st, color='darkgreen', alpha=0.5, lw=2, label='Compton')
            plot.plot(Energy[i], st, color='darkgreen', alpha=0.5, lw=2, label='Tracked Compton')

            plot.scatter(Energy[j],sut,color='blue')
            plot.plot(Energy[j], sut, color='blue', alpha=0.5, lw=2, label='Untracked Compton')

            #plot.scatter(Energy[k],sp,color='darkred')
            #plot.plot(Energy[k],sp, color='darkred', alpha=0.5, lw=2, label='Pair')

            #plot.text(1-0.015, 0.8, u'%i\N{DEGREE SIGN}' % round(numpy.degrees(numpy.arccos(angleSelection))),
             #   	verticalalignment='bottom', horizontalalignment='right',
             #   	transform=ax.transAxes,
             #   	color='black', fontsize=16)

        else:

            angle = round(numpy.degrees(numpy.arccos(angleSelection)))

            plot.scatter(Energy[i],st, color=colors[plotNumber-1])
            plot.plot(Energy[i], st, color=colors[plotNumber-1], alpha=0.5, lw=2, label='Compton at %i$^\circ$' % angle)

            plot.scatter(Energy[j][1:-1], sut[1:-1], color=colors[plotNumber-1])
            plot.plot(Energy[j][1:-1], sut[1:-1], color=colors[plotNumber-1], alpha=0.5, linestyle='-.', label='Compton Untracked at %i$^\circ$' % angle, lw=2)

            #plot.scatter(Energy[j][1:-1], sut[1:-1], color=colors[plotNumber-1])
            #plot.plot(Energy[j][1:-1], sut[1:-1], color=colors[plotNumber-1], alpha=0.5, linestyle='-.', label='Compton Untracked at %i$^\circ$' % angle, lw=2)

            plot.scatter(Energy[k],sp, color=colors[plotNumber-1])
            plot.plot(Energy[k],sp, color=colors[plotNumber-1], alpha=0.5, linestyle='--', lw=2, label='Pair at %i$^\circ$' % angle)


        if plotNumber == len(angleSelections):
            #plot.title('Energy Resolution')			
            plot.legend(numpoints=1, scatterpoints=1, fontsize=16, frameon=True, loc='lower left')
            plot.gca().set_ylim([0.,0.12])
            plot.gca().set_xlim([0.14, 11])

        if xlog:
            plot.xscale('log')
            plot.gca().set_ylim([0.0,0.12])

        if ylog:
            plot.yscale('log')

        #plot.ylabel(r'$\sigma$ / Energy')
        plot.ylabel(r'Energy Resolution (FWHM / Energy)')

        if plotNumber != len(angleSelections):
            ax.set_xticklabels([])

        if plotNumber == len(angleSelections):
            plot.xlabel(r'Energy (MeV)')

        plotNumber = plotNumber + 1


    plot.subplots_adjust(wspace=0, hspace=.2)

    if save:
        plot.savefig('EnergyResolution_Cos%s.pdf' % angleSelections[0])
        plot.savefig('EnergyResolution_Cos%s.png' % angleSelections[0])

    plot.show()

    plot.close()

    return Energy,st,sp

##########################################################################################

def plotEnergyResolutionVsAngle(data, energySelections=None, xlog=False, ylog=False, save=False, collapse=False):

	plotNumber = 1

	if energySelections is None:
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

	if collapse == True:
		plot.figure(figsize=(10, 6.39))
		ax = plot.subplot(111)

	if len(energySelections)>6:
		print("plotting only every other energy")
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

		# Removing nan's
		i=Sigma_tracked == 'nan'
		st=numpy.double(numpy.ma.array(Sigma_tracked,mask=i).compressed())/numpy.double(energySelection)*1.0e-3

		j=Sigma_untracked == 'nan'
		sut=numpy.double(numpy.ma.array(Sigma_untracked,mask=j).compressed())/numpy.double(energySelection)*1e-3		
		
		k=Sigma_pair == 'nan'
		sp=numpy.double(numpy.ma.array(Sigma_pair,mask=k).compressed())/numpy.double(energySelection)*1e-3		

		# Plot the data
		if collapse==False:
			ax = plot.subplot( str(len(energySelections)) + str(10 + plotNumber) )

			plot.plot(numpy.ma.array(Angle,mask=i).compressed(), st, color='darkgreen', alpha=0.75, lw=2, label='Compton', marker='o')
			
			#plot.plot(numpy.ma.array(Angle,mask=j).compressed(), sut, color='blue', alpha=0.75, lw=2, label='Compton (untracked)', marker='o')
			
			plot.plot(numpy.ma.array(Angle,mask=k).compressed(), sp, color='darkred', alpha=0.75, lw=2, label='Pair', marker='o')
			
			plot.text(1-0.015, 0.8, '%s MeV' % energySelection,
		        	verticalalignment='top', horizontalalignment='right',
		        	transform=ax.transAxes,
		        	color='black', fontsize=16)

		else:

			if energySelection<3.:
				#plot.scatter(numpy.ma.array(Angle,mask=i).compressed(), st, colors[plotNumber-1])
				plot.plot(numpy.ma.array(Angle,mask=i).compressed(), st, colors[plotNumber-1], 
					alpha=0.75, lw=2, label='Compton at %s  MeV' % energySelection)
			
			#plot.scatter(numpy.ma.array(Angle,mask=j).compressed(), sut, colors[plotNumber-1])
			#plot.plot(numpy.ma.array(Angle,mask=j).compressed(), sut, colors[plotNumber-1], 
			#	alpha=0.75, lw=2, label='Compton untracked at %s  MeV' % energySelection, linestyle='-.')
			
			if energySelection>3.:
				#plot.scatter(numpy.ma.array(Angle,mask=k).compressed(), sp, colors[plotNumber-1])
				plot.plot(numpy.ma.array(Angle,mask=k).compressed(), sp, colors[plotNumber-1], 
					alpha=0.75, lw=2, label='Pair at %s MeV' % energySelection, linestyle='--')

		if plotNumber == len(energySelections):
			#plot.title('Energy Resolution')						
			plot.legend(numpoints=1, scatterpoints=1, fontsize=16, frameon=True, loc='upper right')

		if xlog:
			plot.xscale('log')

		if ylog:
			plot.yscale('log')

		plot.ylabel(r'$\sigma$ / Energy')
		plot.xlim([0,60])

		if plotNumber == len(energySelections):
			plot.xlabel(r'$\theta$ (deg)')

		plotNumber = plotNumber + 1

	plot.subplots_adjust(wspace=0, hspace=.2)

	if save:
		plot.savefig('EnergyResolutionVsAngle_%sMeV.pdf' % energySelections[0])
		plot.savefig('EnergyResolutionVsAngle_%sMeV.png' % energySelections[0])

	plot.show()

	plot.close()

##########################################################################################

def plotEffectiveArea(data, angleSelections=[1,0.9,0.8,0.7,0.6,0.5], ideal=False, xlog=True, 
                      ylog=False, save=False, show=True, collapse=False, 
                    SurroundingSphere=150, txtOutfileLabel='xxx'):

    if hasattr(angleSelections, '__iter__') == False:
        angleSelections = [angleSelections]

    plotNumber = 1
    if len(angleSelections)==1:
        plot.figure(figsize=(8,6))
    elif collapse == False:
        # Create the new subplot
        plot.figure(figsize=(10,12))
    else:
        plot.figure(figsize=(10, 6.39))
        ax = plot.subplot(111)

    for angleSelection in angleSelections:
        
        results_txt_TC = open( '%s_Aeff_Cos%s_TC.txt' % (txtOutfileLabel, angleSelection), 'w')
        results_txt_UC = open( '%s_Aeff_Cos%s_UC.txt' % (txtOutfileLabel, angleSelection), 'w')
        results_txt_P = open( '%s_Aeff_Cos%s_P.txt' % (txtOutfileLabel, angleSelection), 'w')
        
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
                    # This removes the event selection on the final Aeff calculation
                    # It does not change anything from the FWHM or the 68% containment
                    # Compton events are multiplied by the ratio of tracked vs. untracked
                    total_compton_events = float(data[key][1][-1])
                    if numpy.isnan(total_compton_events):
                        pair_to_total_ratio = 1.0
                    else:
                        pair_to_total_ratio  = float(data[key][4][7])/(float(data[key][4][7])+total_compton_events)
                    #numberOfReconstructedEvents_tracked = 100000.*float(data[key][2][-1])/(total_compton_events)
                    #numberOfReconstructedEvents_untracked = 100000.*float(data[key][3][-1])/(total_compton_events)
                    numberOfReconstructedEvents_tracked = float(data[key][2][-1])
                    numberOfReconstructedEvents_untracked = float(data[key][3][-1])
                    numberOfReconstructedEvents_pair = float(data[key][4][-1]) #float(data[key][4][7])+(float(data[key][4][-1])*pair_to_total_ratio)
                else:
                    numberOfReconstructedEvents_tracked = float(data[key][2][-1])
                    numberOfReconstructedEvents_untracked = float(data[key][3][-1])
                    numberOfReconstructedEvents_pair = float(data[key][4][-2])

                #numberOfReconstructedEvents_tracked = float(data[key][1][-1])
                #numberOfReconstructedEvents_untracked = float(data[key][2][-1])
                #numberOfReconstructedEvents_pair = float(data[key][3][-1])
                print("NumberOfReconstructedEvents Pair", numberOfReconstructedEvents_pair)                
                # Calculate the effective area
                effectiveArea_tracked = (numberOfReconstructedEvents_tracked/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2
                effectiveArea_untracked = (numberOfReconstructedEvents_untracked/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2
                effectiveArea_pair = (numberOfReconstructedEvents_pair/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2

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
    
        #EffectiveArea_UntrackedSiStart=numpy.array([96.0/2884228*70685.8, 14936.0/2132319*70685.8, 21777.0/1744552*70685.8,\
 #                                                   17861.0/1568899*70685.8, 13722.0/1567836*70685.8, #numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, \
#                                                    numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan])
        
        print("EffectiveArea Pair", EffectiveArea_Pair)
        
        if ideal:
            EffectiveArea_Pair[0]=numpy.nan
            EffectiveArea_Pair[1]=numpy.nan

        
        # writing txt files	
        results_txt_TC.write("# Energy[MeV] Aeff_TkrCompton[cm2]\n")
        for ii, en in enumerate(Energy):
        	results_txt_TC.write("%.1f\t%.1f\n"%(en, EffectiveArea_Tracked[ii]))
        results_txt_TC.close()
        print('Created %s_Aeff_Cos%s_TC.txt ...!'%(txtOutfileLabel, angleSelections))
        	
        results_txt_UC.write("# Energy[MeV] Aeff_UntkrCompton[cm2]\n")
        for ii, en in enumerate(Energy):
        	results_txt_UC.write("%.1f\t%.1f\n"%(en, EffectiveArea_Untracked[ii]))
        results_txt_UC.close()
        print('Created %s_Aeff_Cos%s_UC.txt ...!'%(txtOutfileLabel, angleSelection))
        
        results_txt_P.write("# Energy[MeV] Aeff_Pair[cm2]\n")
        for ii, en in enumerate(Energy):
        	results_txt_P.write("%.1f\t%.1f\n"%(en, EffectiveArea_Pair[ii]))
        results_txt_P.close()
        print('Created %s_Aeff_Cos%s_P.txt ...!'%(txtOutfileLabel, angleSelection))
        
        
        # Plot the data
        if collapse == False:
            ax = plot.subplot( str(len(angleSelections)) + str(10 + plotNumber) )

            plot.scatter(Energy, EffectiveArea_Tracked, color='darkgreen')
            #plot.plot(Energy, EffectiveArea_Tracked, color='darkgreen', alpha=0.5, lw=2, label='Compton')

            #plot.scatter(Energy, EffectiveArea_Untracked, color='blue')
            #plot.plot(Energy, EffectiveArea_Untracked, color='blue', alpha=0.5, lw=2, label='Compton (untracked)')
            plot.plot(Energy, EffectiveArea_Tracked, color='darkgreen', alpha=0.5, lw=3, label='Tracked Compton')

            plot.scatter(Energy, EffectiveArea_Untracked, color='blue')
            plot.plot(Energy, EffectiveArea_Untracked, color='blue', alpha=0.5, lw=3, label='Untracked Compton')

            #plot.scatter(Energy, EffectiveArea_UntrackedSiStart, color='blue')
            #plot.plot(Energy, EffectiveArea_UntrackedSiStart, color='blue', linestyle="--", alpha=0.5, lw=3, label='Untracked Compton in Silicon')

            plot.scatter(Energy, EffectiveArea_Pair, color='darkred')
            #plot.plot(Energy, EffectiveArea_Pair, color='darkred', alpha=0.5, lw=2, label='Pair')
            plot.plot(Energy, EffectiveArea_Pair, color='darkred', alpha=0.5, lw=3, label='Pair')

            #plot.text(1-0.015, 0.8, u'%i\N{DEGREE SIGN}' % round(numpy.degrees(numpy.arccos(angleSelection))),
             #   verticalalignment='bottom', horizontalalignment='right',
             #   transform=ax.transAxes,
             #   color='black', fontsize=16)	

        else:
            angle = round(numpy.degrees(numpy.arccos(angleSelection)))
            plot.scatter(Energy[1:], EffectiveArea_Tracked[1:], color='darkgreen')
            plot.plot(Energy[1:], EffectiveArea_Tracked[1:], color='darkgreen', alpha=0.5, lw=3, label='Tracked Compton at %i$^\circ$' % angle)

            #plot.scatter(Energy, EffectiveArea_Untracked, color=colors[plotNumber-1])
            #plot.plot(Energy, EffectiveArea_Untracked, color=colors[plotNumber-1], lw=2, alpha=0.5, linestyle='-.')
            plot.scatter(Energy, EffectiveArea_Untracked, color=colors[plotNumber-1])
            plot.plot(Energy, EffectiveArea_Untracked, color=colors[plotNumber-1], lw=3, alpha=0.5, linestyle='-.', label='Untracked Compton at %i$^\circ$' % angle)


            plot.scatter(Energy[5:], EffectiveArea_Pair[5:], color=colors[plotNumber-1])
            #plot.plot(Energy[5:], EffectiveArea_Pair[5:], color=colors[plotNumber-1], alpha=0.5, lw=2, linestyle='--', label='Pair at %i$^\circ$' % angle)
            plot.plot(Energy[5:], EffectiveArea_Pair[5:], color=colors[plotNumber-1], alpha=0.5, lw=3, linestyle='--', label='Pair at %i$^\circ$' % angle)


        if plotNumber == len(angleSelections):
            #plot.title('Effective Area')			
            #plot.legend(numpoints=1, scatterpoints=1, fontsize=16, frameon=True, loc='upper left')
            plot.legend(numpoints=1, scatterpoints=1, fontsize=16, frameon=True, loc='lower right')
            
        #plot.ylabel(r'A$_{\mathrm{eff}}$ (cm$^2$)')
        plot.ylabel('Effective Area (cm$^2$)')


        if xlog:
            plot.xscale('log')
            plot.gca().set_xlim([0.1, 10000])

        if ylog:
            plot.gca().set_ylim([1, 4000.])
            plot.yscale('log')
        else:
            plot.gca().set_ylim([0.001,6000.])
            plot.gca().set_xlim([0.12, 10000])


        if plotNumber != len(angleSelections):
            ax.set_xticklabels([])

        if plotNumber == len(angleSelections):
            plot.xlabel(r'Energy (MeV)')

        plotNumber = plotNumber + 1


    plot.subplots_adjust(wspace=0, hspace=.2)

    if save:
        plot.savefig('EffectiveArea_Cos%s.pdf' % angleSelections[0], bbox_inches='tight')
        plot.savefig('EffectiveArea_Cos%s.png' % angleSelections[0], bbox_inches='tight')

    if show:
        plot.show()

    return EffectiveArea_Untracked, EffectiveArea_Tracked, EffectiveArea_Pair, Energy

##########################################################################################

def tabulateEffectiveArea(data, angleSelections=[1,0.9,0.8,0.7,0.6,0.5], ideal=False, doPrint= True, SurroundingSphere=150):

    if hasattr(angleSelections, '__iter__') == False:
        angleSelections = [angleSelections]

    angleSelections = set(angleSelections)
        
    energy_data = dict()
        
    for key in data.keys():
        energy = float(key.split('_')[1].replace('MeV', ''))
        half = key.split('_')[2].replace('Cos', '')
        angle = float(half.replace('.inc1.id1.sim', ''))

        #print energy, angle

        if angle not in angleSelections:
            continue

        if energy not in energy_data:
            energy_data[energy] = dict()
            
        numberOfSimulatedEvents = float(data[key][0])

        if ideal:
            # This removes the event selection on the final Aeff calculation
            # It does not change anything from the FWHM or the 68% containment
            # Compton events are multiplied by the ratio of tracked vs. untracked
            total_compton_events = float(data[key][1][-1])
            if numpy.isnan(total_compton_events):
                pair_to_total_ratio = 1.0
            else:
                pair_to_total_ratio  = float(data[key][4][7])/(float(data[key][4][7])+total_compton_events)
            #numberOfReconstructedEvents_tracked = 100000.*float(data[key][2][-1])/(total_compton_events)
            #numberOfReconstructedEvents_untracked = 100000.*float(data[key][3][-1])/(total_compton_events)
            numberOfReconstructedEvents_tracked = float(data[key][2][-1])
            numberOfReconstructedEvents_untracked = float(data[key][3][-1])
            numberOfReconstructedEvents_pair =  float(data[key][4][-2])#float(data[key][4][7])+(float(data[key][4][-1])*pair_to_total_ratio)
        else:
            numberOfReconstructedEvents_tracked = float(data[key][2][-1])
            numberOfReconstructedEvents_untracked = float(data[key][3][-1])
            numberOfReconstructedEvents_pair = float(data[key][4][-2])

		#numberOfReconstructedEvents_tracked = float(data[key][1][-1])
		#numberOfReconstructedEvents_untracked = float(data[key][2][-1])
		#numberOfReconstructedEvents_pair = float(data[key][3][-1])

		# Calculate the effective area
        effectiveArea_tracked = (numberOfReconstructedEvents_tracked/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2
        effectiveArea_untracked = (numberOfReconstructedEvents_untracked/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2
        effectiveArea_pair = (numberOfReconstructedEvents_pair/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2
            
		# Store the results
        energy_data[energy][angle] = [ effectiveArea_tracked, effectiveArea_untracked, effectiveArea_pair ]

    # create a table that is sorted first by energy, and second by angle
    angleSelections = sorted(list(angleSelections))
    table_data = sorted(
        [ ( energy, sorted([ (angle, d) for angle, d in ds.iteritems() ], key = lambda entry: entry[0]) )
          for energy, ds in energy_data.iteritems() ],
        key = lambda entry: entry[0])

    # print that table
    if doPrint:
    	print('Energy    | %s' % ('                      | '.join(['% 24.1f' % x for x in angleSelections])))
    	print('         ' + ' |  tracked  | untracked |  compton  |    pair  ' * len(angleSelections))
    	print('------------' + '-' * 48*len(angleSelections))
    	for energy, ds in table_data:
        	print('% 8.1f  | %s' % (energy, ' | '.join(['% 9.1f | % 9.1f | % 9.1f | % 9.1f' % (x[1][0], x[1][1], x[1][0] + x[1][1], x[1][2]) for x in ds])))

    return table_data


##########################################################################################

##########################################################################################

def plotEffectiveAreaVsAngle(data, energySelections=None, ideal=False, xlog=False, ylog=False, save=False, collapse=False, SurroundingSphere=150):

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
		print("plotting only every other energy")
		energySelections=energySelections[[1,3,5,7,9,11]]

	plotNumber = 1
	if len(energySelections)==1:
		plot.figure(figsize=(8,6))
	elif collapse == False:
		# Create the new subplot
		plot.figure(figsize=(10,12))
	else:
		plot.figure(figsize=(10, 6.39))
		ax = plot.subplot(111)

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
					if numpy.isnan(total_compton_events):
						pair_to_total_ratio = 1.0
					else:
						pair_to_total_ratio  = float(data[key][4][7])/(float(data[key][4][7])+total_compton_events)
					numberOfReconstructedEvents_tracked = 100000.*float(data[key][2][-1])/(total_compton_events)
					numberOfReconstructedEvents_untracked = 100000.*float(data[key][3][-1])/(total_compton_events)
					numberOfReconstructedEvents_pair = float(data[key][4][-2])# float(data[key][4][7])+(float(data[key][4][-1])*pair_to_total_ratio)
				else:
					numberOfReconstructedEvents_tracked = float(data[key][2][-1])
					numberOfReconstructedEvents_untracked = float(data[key][3][-1])
					numberOfReconstructedEvents_pair = float(data[key][4][-2])

				numberOfSimulatedEvents = float(data[key][0])

				# Calculate the effective area
				effectiveArea_tracked = (numberOfReconstructedEvents_tracked/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2
				effectiveArea_untracked = (numberOfReconstructedEvents_untracked/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2
				effectiveArea_pair = (numberOfReconstructedEvents_pair/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2

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

		# Plot the data
		if collapse == False:
			ax = plot.subplot( str(len(energySelections)) + str(10 + plotNumber) )

			plot.scatter(Angle, EffectiveArea_Tracked, color='darkgreen')
			plot.plot(Angle, EffectiveArea_Tracked, color='darkgreen', alpha=0.5, lw=2, label='Compton')

			#plot.scatter(Angle, EffectiveArea_Untracked, color='blue')
			#plot.plot(Angle, EffectiveArea_Untracked, color='blue', alpha=0.5, lw=2, label='Compton (untracked)')

			plot.scatter(Angle, EffectiveArea_Pair, color='darkred')
			plot.plot(Angle, EffectiveArea_Pair, color='darkred', alpha=0.5, lw=2, label='Pair')			

			plot.text(1-0.015, 0.8, '%s MeV' % energySelection,
			        verticalalignment='bottom', horizontalalignment='right',
		       		transform=ax.transAxes,
			        color='black', fontsize=16)

		else:
			if energySelection<3.:
				plot.scatter(Angle, EffectiveArea_Tracked, color=colors[plotNumber-1])
				plot.plot(Angle, EffectiveArea_Tracked, color=colors[plotNumber-1], alpha=0.5, lw=2, label='Compton at %s MeV' % energySelection)

			#plot.scatter(Angle, EffectiveArea_Untracked, color=colors[plotNumber-1])
			#plot.plot(Angle, EffectiveArea_Untracked, color=colors[plotNumber-1], alpha=0.5, lw=2, linestyle='-.')

			if energySelection>3.:
				plot.scatter(Angle, EffectiveArea_Pair, color=colors[plotNumber-1])
				plot.plot(Angle, EffectiveArea_Pair, color=colors[plotNumber-1], alpha=0.5, lw=2, linestyle='--', label='Pair at %s MeV' % energySelection)


		if plotNumber == len(energySelections):
			#plot.title('Effective Area')			
			plot.legend(numpoints=1, scatterpoints=1, fontsize=16, frameon=True, loc='upper left')

		#plot.ylabel(r'A$_{\mathrm{eff}}$ (cm$^2$)')
		plot.ylabel('Effective Area (cm$^2$)')

		if xlog:
			plot.xscale('log')

		if ylog:
			plot.yscale('log')
			#if len(energySelections)==1:
			plot.gca().set_ylim([100.,5000.])


		if plotNumber != len(energySelections) and collapse == False:
			ax.set_xticklabels([])

		if plotNumber == len(energySelections) or collapse == True:
			plot.xlabel(r'$\theta$')

		plotNumber = plotNumber + 1


	plot.subplots_adjust(wspace=0, hspace=.2)

	if save:
		plot.savefig('EffectiveAreaVsAngle_%sMeV.pdf' % energySelections[0])
		plot.savefig('EffectiveAreaVsAngle_%sMeV.png' % energySelections[0])

	plot.show()

##########################################################################################
def resultsToFits(data, outfile='output.fits'):

	Energy = []	
	Angle = []
	aeff_vs_energy = []
	for key in data.keys():
		energy = float(key.split('_')[1].replace('MeV',''))
		half = key.split('_')[2].replace('Cos','')
		angle = float(half.replace('.inc1.id1.sim',''))
		if energy not in Energy:
			Energy.append(energy)
		if angle not in Angle:
			Angle.append(angle)
	Energy=numpy.array(Energy)
	Angle=numpy.array(Angle)
	i = [numpy.argsort(Energy)]
	j = [numpy.argsort(Angle)]

	aeff_vs_energy[0] = plotEffectiveArea(data, angleSelections=Angle[j], ideal=True, show=False)
	aeff_vs_energy[0] = plotEffectiveArea(data, angleSelections=Angle[j], ideal=True, show=False)
	
	print(aeff_vs_energy[0][0])

	#print energySelections, angleSelections


	#t.write(outfile,format='fits',overwrite=True)

########################################################################################## 

def plotSourceSensitivity(data, angleSelection=0.8, exposure = 3.536*10**7, ideal=False, doPSF=None, \
    xlog=True, ylog=True, save=False, doplot=False, showbackground=False, uniterg=False,doRealBkg=True, \
    SurroundingSphere=150):

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
                # This removes the event selection on the final Aeff calculation
                # It does not change anything from the FWHM or the 68% containment
                # Compton events are multiplied by the ratio of tracked vs. untracked
                total_compton_events = float(data[key][1][-1])
                if numpy.isnan(total_compton_events):
                    pair_to_total_ratio = 1.0
                else:
                    pair_to_total_ratio  = float(data[key][4][7])/(float(data[key][4][7])+total_compton_events)
                if total_compton_events == 0.:
                    total_compton_events=1.

                numberOfReconstructedEvents_tracked = 100000.*float(data[key][2][-1])/(total_compton_events)
                numberOfReconstructedEvents_untracked = 100000.*float(data[key][3][-1])/(total_compton_events)
                numberOfReconstructedEvents_pair = float(data[key][4][-4])#+(float(data[key][4][-1])*pair_to_total_ratio)
            else:
                numberOfReconstructedEvents_tracked = float(data[key][2][-1])
                numberOfReconstructedEvents_untracked = float(data[key][3][-1])
                numberOfReconstructedEvents_pair = float(data[key][4][-2])

			# Get the angular resolution
            fwhm_tracked = data[key][2][7]
            fwhm_untracked = data[key][3][7]
            containment68 = data[key][4][6]

            # Calculate the effective area
            effectiveArea_tracked = (numberOfReconstructedEvents_tracked/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2
            effectiveArea_untracked = (numberOfReconstructedEvents_untracked/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2
            effectiveArea_pair = (numberOfReconstructedEvents_pair/numberOfSimulatedEvents) * math.pi * SurroundingSphere**2

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
    
    #Note: this only changes the PSF for Pair events, not Compton!
    if doPSF != None:
        Containment68 = doPSF

    # digitized from Stong, Moskalenko & Reimer 2000, Figure 8, top right
    # high latitude |b|>5 deg
    # multiply by 10 to vaguely account for the albedo background
    oldeng2=numpy.array([0.10355561,0.3534914,1.2920963,4.659387,8.969312,18.735151,38.081676,69.40132,\
        144.98259,227.4451,342.42523,462.24567,725.01324,939.413,1908.1061,28725.793])
    olde2int2=numpy.array([2.7943178E-4,3.57757E-4,4.8821748E-4,6.806025E-4,8.0072926E-4,9.1560354E-4,\
        0.0010469892,0.0011638523,0.0013691497,0.0015439879,0.0016334692,0.0017039803,0.0018284274,\
        0.0018672496,0.0017879958,0.0014717471])

    # Alex's new background numbers from Gruber et al. (1999) and Weidenspointer et al. (2000) and >100 MeV from Ackermann et al. (2015)
    alex_eng2=numpy.array([0.5,0.8,1.0,2.0,3.0,5.0,8.0,10.0,50.0,100.,200,500])
    alex_e2int2=numpy.array([2e-2,1e-2,7e-3,3e-3,2e-3,8e-4,4e-4,3e-4,3e-5,3.2e-6,2e-6,6e-7])*alex_eng2

    # From Acero et al. (2016) - arxiv:1602.07246 |b| > 10 deg Galactic Diffuse
    lateng=numpy.array([58.665302,83.7944,127.701385,212.20918,296.02475,493.0605,740.34045,1265.6293,\
                        2019.2109,3006.2268,4828.027,8546.594,18742.852,42185.098,152450.55,496614.97])
    late2int_galactic=numpy.array([8.653016E-4,0.0011343559,0.0015828605,0.0020333533,0.0022578337,\
                                   0.002416496,0.0023796277,0.002305653,0.0019558307,0.0016045898,0.0011626304,7.918069E-4,\
                                   4.5331568E-4,2.5003447E-4,1.3304557E-4,7.2556504E-5])
    # From Ackermann et al. (2015) - ApJ 799 86 isotropic EGB
    lateng_igrb=numpy.array([120,170,240,340,490,690,900,1350,1950,2750,3850,5450,7750,11050,15500,22200,\
                             31000,43500,61500,86000,120000,170000,245000,350000,495000,700000.])
    lat_igrb=numpy.array([3.7e-6,2.3e-6,1.5e-6,9.7e-7,6.7e-7,4.9e-7,3e-7,1.8e-7,1.1e-7,6.9e-8,4.2e-8,\
                          2.6e-8,1.7e-8,1.2e-8,6.8e-9,4.4e-9,2.7e-9,1.8e-9,1.1e-9,6.2e-10,3.1e-10,1.9e-10,8.9e-11,6.3e-11,\
                          2.1e-11,9.7e-12])
    late2int_igrb0=lat_igrb*lateng_igrb
    lateng_igrb.sort()
    late2int_igrb0.sort()
    tck=interpolate.splrep(numpy.log10(lateng_igrb),numpy.log10(late2int_igrb0),s=0)
    late2int_igrb=10**interpolate.splev(numpy.log10(lateng),tck,der=0)
    late2int=late2int_galactic+late2int_igrb

    # COMPTEL * EGRET Galactic Diffuse from Gruber et al. (1999)
    gruber_eng=numpy.array([2.978623,5.1983213,9.07216,13.32116,19.94295,32.241817,44.707794,72.33151,\
        136.47008,278.06522,545.20044,1132.5265,3079.0847,6774.5522,17384.865,41301.277,105963.19,\
        317014.44,1315024.2,6868901.5,2.2191038E7,8.6879376E7])*1e-3
    gruber_e2int=numpy.array([5.278219,4.1341214,3.2380166,2.592378,1.8563249,1.2433603,0.8325035,\
        0.36483333,0.13988705,0.049062684,0.01799201,0.00590178,0.0014491071,4.9711246E-4,1.3645743E-4,\
        4.6819663E-5,1.4694342E-5,4.219293E-6,8.482257E-7,1.4922263E-7,4.098341E-8,9.63159E-9])*gruber_eng

    # COMPTEL extragalactic background from Weidenspointner (2001)
    wp_eng=numpy.array([0.10258658,0.15835446,0.34505856,0.56054914,0.9957529,1.9590781,\
                        3.4359305,10.434804,57.75009,135.85852,377.31998])
    wp_e2int=numpy.array([0.024865912,0.016798664,0.010788648,0.0067556994,0.00433873,0.002821951,\
                          0.0022472343,0.0021910856,0.002363885,0.0015768929,0.0014071865])

    ## Combining things
    eng2=numpy.append(gruber_eng[0:16],lateng[2:])
    e2int2=numpy.append(gruber_e2int[0:16],late2int[2:])
    # interpolate background at our energies
    tck=interpolate.splrep(numpy.log10(eng2),numpy.log10(e2int2),s=0)
    logbackground=interpolate.splev(numpy.log10(Energy),tck,der=0)
    background=10.**logbackground


    if showbackground:
        import ScienceSims
        if doPSF is None:
            hist_CompARM,hist_PairARM=ScienceSims.plot_AMEGO_background_sim(dir='../Simulations/BackgroundFiles/Sims_100s/',angleSelection=angleSelection,energy=Energy,st=FWHM_tracked,sp=Containment68)

        plot.figure()
        plot.plot(oldeng2,olde2int2,color='red')
        plot.scatter(oldeng2,olde2int2,color='red')
        plot.annotate('Strong, Moskalenko, Reimer (2000)',xy=(1e-2,1e-4),xycoords='data',fontsize=12,color='red')
        #plot.plot(alex_eng2,alex_e2int2,color='blue')
        #plot.scatter(alex_eng2,alex_e2int2,color='blue')
        #plot.annotate('Alex (new)',xy=(1e2,2e-4),xycoords='data',fontsize=12,color='blue')
        plot.plot(lateng,late2int,color='magenta')
        plot.scatter(lateng,late2int,color='magenta')
        plot.plot(lateng,late2int_galactic,'r--',color='magenta')
        plot.scatter(lateng,late2int_galactic,color='magenta')
        plot.plot(lateng,late2int_igrb,'r:',color='magenta',lw=2)
        plot.scatter(lateng,late2int_igrb,color='magenta')
        plot.annotate('Ackermann et al. (2016) - LAT Extragalactic',xy=(5,1e-5),xycoords='data',fontsize=12,color='magenta')
        plot.annotate('Acero et al. (2016) - LAT Galactic',xy=(1e2,3e-3),xycoords='data',fontsize=12,color='magenta')
        plot.plot(gruber_eng,gruber_e2int,color='blue')
        plot.scatter(gruber_eng,gruber_e2int,color='blue')
        plot.annotate('Gruber et al. (1999) - HEAO, COMPTEL, EGRET',xy=(2e-3,5e-2),xycoords='data',fontsize=12,color='blue')
        plot.plot(wp_eng,wp_e2int,'r:',color='orange')
        plot.scatter(wp_eng,wp_e2int,color='orange')
        plot.annotate('Weidenspointner (2001)',xy=(2e-3,1e-3),xycoords='data',fontsize=12,color='orange')
        #plot.plot(eng2,e2int2,color='purple')
        #plot.scatter(eng2,e2int2,color='purple')	
        plot.plot(Energy,background,color='green')
        plot.scatter(Energy,background,color='green')

        if doPSF is None:
            l=len(hist_CompARM[0])
            wc=numpy.where(hist_CompARM[1] != 0.)
            wp=numpy.where(hist_PairARM[1] != 0.)
            plot.plot(hist_CompARM[0][wc],hist_CompARM[0][wc]*hist_CompARM[1][wc]/omega(FWHM_untracked[wc])/EffectiveArea_Untracked[wc],color='brown')
            plot.plot(hist_PairARM[0][wp],hist_PairARM[0][wp]*hist_PairARM[1][wp]/omega(Containment68[wp])/EffectiveArea_Pair[wp],color='brown')

        import ScienceSims
        if doRealBkg:
            real_back_eng,real_back_rate=ScienceSims.plot_AMEGO_background_sim(dir='../Simulations/BackgroundFiles/Sims_1000s/',doplot=False)
            plot.plot(real_back_eng,real_back_rate/omega(FWHM_untracked)/EffectiveArea_Untracked,color='brown')

        plot.annotate('Interpolated Bkg',xy=(1,1e-2),xycoords='data',fontsize=12,color='green')
        plot.xscale('log')
        plot.yscale('log')
        plot.gca().set_ylim([0.000005,0.1])
        plot.xlabel(r'Energy (MeV)')
        plot.ylabel(r'$E^2 \times$ Intensity (MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
        plot.title('Diffuse Background')
        plot.savefig('plots/Background.pdf', bbox_inches='tight')
        plot.savefig('plots/Background.png', bbox_inches='tight')
        plot.show()

    Sensitivity_tracked = Isrc(Energy, exposure, EffectiveArea_Tracked, 3., omega(FWHM_tracked), background*10.)
    Sensitivity_untracked = Isrc(Energy, exposure, EffectiveArea_Untracked, 3., omega(FWHM_untracked), background*50.)
    Sensitivity_pair = Isrc(Energy, exposure, EffectiveArea_Pair, 3., omega(Containment68), background)

    if doPSF:
        Sensitivity_pair[0]=numpy.nan
        Sensitivity_pair[1]=numpy.nan

    if doplot:

        if plot.fignum_exists(1):
            plot.clf()
        else:
            plot.figure()#figsize=(10, 6.39))
            ax = plot.subplot(111)

        if uniterg:
            mev2erg=1/624151.
            unit='erg'
            ylim=[1e-12,1e-9]
        else:
            mev2erg=1
            unit='MeV'
            ylim=[1e-7,5e-5]

        plot.scatter(Energy, Sensitivity_tracked*mev2erg, color='darkgreen')
        plot.plot(Energy, Sensitivity_tracked*mev2erg, color='darkgreen', alpha=0.5, label='Compton')
        plot.scatter(Energy, Sensitivity_untracked*mev2erg, color='blue')	
        plot.plot(Energy, Sensitivity_untracked*mev2erg, color='blue', alpha=0.5, label='Compton (untracked)')
        plot.scatter(Energy, Sensitivity_pair*mev2erg, color='darkred')
        plot.plot(Energy, Sensitivity_pair*mev2erg, color='darkred', alpha=0.5, label='Pair')	
        plot.xlabel(r'Energy (MeV)')
        plot.ylabel(r'$E^2 \times$ Intensity ('+unit+' cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
        plot.ylim(ylim)

        if xlog:
            plot.xscale('log')

        if ylog:
            plot.yscale('log')

        plot.legend(numpoints=1, scatterpoints=1, fontsize=16, frameon=True, loc='upper left')

        if save:
            plot.savefig('SourceSensitivity.pdf', bbox_inches='tight')
            plot.savefig('SourceSensitivity.png', bbox_inches='tight')

    if doplot:
        plot.show()

    return Energy, Sensitivity_tracked, Sensitivity_untracked, Sensitivity_pair 

##########################################################################################

def plotAllSourceSensitivities(data, angleSelection=0.8, plotIdeal=False, xlog=True, ylog=True, \
    save=False, showbackground=False, doplot=True, uniterg=False):

    ComPairSensitivity=plotSourceSensitivity(data,angleSelection=angleSelection,doplot=False)
    ComPairIdealSensitivity=plotSourceSensitivity(data,angleSelection=angleSelection,ideal=True,doplot=False, showbackground=showbackground)
    ComPairGoodPSFSensitivity=plotSourceSensitivity(data,angleSelection=angleSelection,ideal=True,doPSF=2.0,doplot=False)

    if uniterg:
        mev2erg=1/624151.
        unit='erg'
        ylim=[1e-14,1e-8]
    else:
        mev2erg=1
        unit='MeV'
        ylim=[1e-8,1e-2]

    plot.figure(figsize=(10, 6.39))

    a=ascii.read("data/digitized_alex_sensitivities.dat",names=['eng','sensit'])
    l=ascii.read("data/differential_flux_sensitivity_p8r2_source_v6_all_10yr_zmax100_n10.0_e1.50_ts25_000_090.txt",names=['emin','emax','e2diff','tmp'])

    #print a['energy']
    energy=a['eng']
    sens=a['sensit']

    erg2mev=624151.
    lateng=(l["emin"]+l["emax"])/2.

    plot.clf()
    #LAT
    lat_energy=numpy.array([31.636559,55.997124,99.194466,178.46051,315.91608,558.86786,1004.1673,\
        1801.5197,3179.561,5600.7397,10030.165,17674.44,31625.436,54096.105,98242.1,173146.97,305125.9,\
        546138.4,990550.4,1712016.6,3043571.2])
    lat_Aeff=numpy.array([0.07349323,0.19316217,0.3357588,0.48363942,0.6068357,0.71063167,0.8214752,\
        0.88822705,0.9250036,0.90534276,0.9209482,0.9118694,0.9027835,0.92017376,0.89697146,0.89318365,\
        0.8858685,0.8856009,0.8253616,0.6998877,0.50385296])*1e4
    lat_psf_energy=numpy.array([32.03726,56.471546,99.54146,178.16881,304.57852,545.1637,990.85016,\
        1720.002,3078.6238,5595.484,9863.066,17385.465,30179.168,53196.316,95215.85,173057.44,300407.8,\
        529523.9,947792.44,1670658.0,2944841.2])
    lat_psf0=numpy.array([11.73808,7.8302007,5.135973,3.368779,2.1726828,1.3098335,0.83064204,0.5179475,\
        0.32846078,0.23048219,0.16448157,0.1255744,0.106081836,0.100847006,0.1025626,0.1025626,0.095870495,\
        0.08961505,0.082366556,0.07570436,0.07196857])
    tck=interpolate.splrep(numpy.log10(lat_psf_energy),numpy.log10(lat_psf0),s=0)
    lat_psf=10**interpolate.splev(numpy.log10(lat_energy),tck,der=0)

    #LAT background
    #  Galactic from Acero et al. (2016)
    lateng0=numpy.array([58.665302,83.7944,127.701385,212.20918,296.02475,493.0605,740.34045,1265.6293,\
        2019.2109,3006.2268,4828.027,8546.594,18742.852,42185.098,152450.55,496614.97])
    late2int=numpy.array([8.653016E-4,0.0011343559,0.0015828605,0.0020333533,0.0022578337,0.002416496,\
        0.0023796277,0.002305653,0.0019558307,0.0016045898,0.0011626304,7.918069E-4,4.5331568E-4,\
        2.5003447E-4,1.3304557E-4,7.2556504E-5])
    tck=interpolate.splrep(numpy.log10(lateng0),numpy.log10(late2int),s=0)
    lat_background_galactic=10**interpolate.splev(numpy.log10(lat_energy),tck,der=0)

    #  Extragalactic from Ackermann et al. (2015) - ApJ 799 86
    lateng_igrb=numpy.array([120,170,240,340,490,690,900,1350,1950,2750,3850,5450,7750,11050,15500,22200,\
        31000,43500,61500,86000,120000,170000,245000,350000,495000,700000.])
    lat_igrb=numpy.array([3.7e-6,2.3e-6,1.5e-6,9.7e-7,6.7e-7,4.9e-7,3e-7,1.8e-7,1.1e-7,6.9e-8,4.2e-8,\
        2.6e-8,1.7e-8,1.2e-8,6.8e-9,4.4e-9,2.7e-9,1.8e-9,1.1e-9,6.2e-10,3.1e-10,1.9e-10,8.9e-11,6.3e-11,\
        2.1e-11,9.7e-12])
    late2int_igrb0=lat_igrb*lateng_igrb
    tck=interpolate.splrep(numpy.log10(lateng_igrb),numpy.log10(late2int_igrb0),s=0)
    late2int_igrb=10**interpolate.splev(numpy.log10(lat_energy),tck,der=0)
    lat_background=lat_background_galactic+late2int_igrb

    lat_exposure=5.*365.*86400.*0.2 # 5 years, 20% of sky in FoV, assuming minimal SAA
    lat_sensitivity=Isrc(lat_energy,lat_exposure,lat_Aeff,3,omega(lat_psf),lat_background)
    plot.plot(lat_energy,lat_sensitivity*mev2erg,color='magenta',lw=2)
    #plot.plot(lateng,l["e2diff"]*erg2mev,color='magenta',lw=2)
    plot.gca().set_xscale('log')
    plot.gca().set_yscale('log')
    plot.gca().set_xlim([1e-2,1e6])
    plot.gca().set_ylim(ylim)
    plot.gca().set_xlabel('Energy (MeV)', fontsize=16, fontname="Helvetica")
    plot.gca().set_ylabel(r'$3\sigma$ Continuum Sensitivity $\times\ E^2$ ($\gamma$ '+unit+' cm$^{-2}$ s$^{-1}$)', fontsize=16)
    plot.annotate('Fermi-LAT', xy=(1e3,(1e-6*mev2erg)),xycoords='data',fontsize=16,color='magenta')

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
    plot.plot(egret_energy,egret_sensitivity*mev2erg,color='blue',lw=2)
    plot.annotate('EGRET', xy=(4e2,(1e-4*mev2erg)),xycoords='data',fontsize=16,color='blue')

    #SPI
    ind=numpy.arange(20,46)
    plot.plot(energy[ind],sens[ind]*mev2erg,color='green',lw=2)
    plot.annotate('SPI', xy=(6e-2,(1e-4*mev2erg)),xycoords='data',fontsize=16,color='green')

    #COMPTEL
    comptel_energy=numpy.array([0.73295844,0.8483429,1.617075,5.057877,16.895761,29.717747])
    comptel_sens=numpy.array([6.566103E-4,3.6115389E-4,1.4393721E-4,1.6548172E-4,2.36875E-4,3.390693E-4])
    plot.plot(comptel_energy,comptel_sens*mev2erg,color='orange',lw=2)
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
    plot.annotate('COMPTEL', xy=(5,(5e-4*mev2erg)),xycoords='data',fontsize=16,color='orange')

    #NuSTAR
    ind=numpy.arange(84,147)
    plot.plot(energy[ind]*1e-3,sens[ind]*(energy[ind]/1e3)**2*1e3*mev2erg,color='purple',lw=2)
    plot.annotate('NuSTAR', xy=(0.1,(3e-8*mev2erg)),xycoords='data',fontsize=16,color='purple')

    #ComPair
    compair_eng=ComPairSensitivity[0]
    tracked=ComPairSensitivity[1]	
    untracked=ComPairSensitivity[2]
    pair=ComPairSensitivity[3]

    if plotIdeal:
        tracked_ideal=ComPairIdealSensitivity[1]
        untracked_ideal=ComPairIdealSensitivity[2]
        #pair_ideal=ComPairIdealSensitivity[3]
        pair=ComPairIdealSensitivity[3]
        #pair_idealang=ComPairGoodPSFSensitivity[3]

    combined=numpy.zeros(len(compair_eng))

    for e in range(len(compair_eng)):
        print("untracked and tracked and pair", untracked[e], tracked[e], pair[e])
        if numpy.isnan(pair[e]) and numpy.isnan(tracked[e]):
            combined[e]=untracked[e]
        if numpy.isnan(pair[e]) and numpy.isfinite(tracked[e]):
            if untracked[e] > tracked[e]:
                combined[e]=tracked[e]*(untracked[e]/(untracked[e]+tracked[e]))**0.5
            if untracked[e] < tracked[e]:
                combined[e]=untracked[e]*(tracked[e]/(tracked[e]+untracked[e]))**0.5
        if numpy.isfinite(pair[e]) and numpy.isfinite(tracked[e]):
            if tracked[e] > pair[e]:
                combined[e]=pair[e]*(tracked[e]/(tracked[e]+pair[e]))**0.5
                #print (tracked[e]/(tracked[e]+pair[e]))**0.5
            if tracked[e] < pair[e]:
                combined[e]=tracked[e]*(pair[e]/(tracked[e]+pair[e]))**0.5
                #print (pair[e]/(tracked[e]+pair[e]))**0.5
        if numpy.isnan(untracked[e]) and numpy.isnan(pair[e]):
            combined[e]=tracked[e]
        if numpy.isnan(untracked[e]) and numpy.isnan(tracked[e]):
            combined[e]=pair[e]

        print("combined", combined[e])



    #Combine with activation simulation numbers
    amego_activation_interp = numpy.interp(compair_eng[(compair_eng>0.5) & (compair_eng<=10)],amego_activation_energy, amego_activation_sens)
    #plot.plot(compair_eng[(compair_eng>0.3) & (compair_eng<=10)], amego_activation_interp,color='blue', lw=3, linestyle='--')
    for e in range(len(compair_eng[(compair_eng>0.5) & (compair_eng<=10)])):
        combined[e+2]= (amego_activation_interp[e] + combined[e+2])/2

    #plot.plot(compair_eng,tracked,'r--',color='grey',lw=2)	
    #plot.plot(compair_eng,pair,'r--',color='grey',lw=2)
    plot.plot(compair_eng[(compair_eng>0.1) & (compair_eng<=1e3)],combined[(compair_eng>0.1) & (compair_eng<=1e3)]*mev2erg,color='black',lw=4)

    print("Final AMEGO Numbers, energy:", compair_eng[(compair_eng>0.1) & (compair_eng<=1e3)])
    print("Final AMEGO Numbers, energy:", combined[(compair_eng>0.1) & (compair_eng<=1e3)]*mev2erg )

    #if plotIdeal:
        #plot.plot(compair_eng,tracked_ideal,'r:',color='black',lw=3)
        ##plot.plot(compair_eng,pair_ideal,'r:',color='black',lw=3)
        #plot.plot(compair_eng,pair_idealang,'r-.',color='black',lw=2)

    #Alex's numbers
    #ind=numpy.arange(158,167)
    #plot.plot(energy[ind],sens[ind],'r--',color='red',lw=2)
    #plot.annotate('Previous ComPair', xy=(7e2,3e-6),xycoords='data',fontsize=14,color='red')

    #plot.plot(compair_eng,pair_idealang,'r-.',color='black',lw=3)
    plot.annotate('AMEGO', xy=(1,(3e-7*mev2erg)),xycoords='data',fontsize=22)

    if save:
        plot.savefig('full_sensitivity_Cos%s.pdf' % angleSelection)
        plot.savefig('full_sensitivity_Cos%s.png' % angleSelection)
        fout=open('AMEGO_sensitivity_Cos%s.dat' %angleSelection,"w")
        fout.write('Energy (MeV)'+"\t"+ 'E2xSensitivity ('+unit+'cm-2 s-1)')
        for i in range(len(compair_eng[(compair_eng>0.1) & (compair_eng<=1e3)])):
            fout.write(str(compair_eng[(compair_eng>0.1) & (compair_eng<=1e3)][i])+"\t"+\
                str(combined[(compair_eng>0.1) & (compair_eng<=1e3)]*mev2erg))

    if doplot == True:
        plot.show()

    return compair_eng,combined

##########################################################################################

def applyARM(data,events,angleSelection=1.0,ARMcut=None,energy=None,st=None,sp=None):

    import EventAnalysis
    """
    A function to select events within ARM defined by simulatinos
    Example Usage: 
    eventsInARM=EventViewer.applyARM(data,events)
    """
    from scipy.interpolate import interp1d

    angle=angleSelection  #not sure if this should vary simply with angle given background Sims
    sourceTheta=numpy.arccos(angle)
    if (energy == None) & (st==None) & (sp==None):
        energy,st,sp=plotAngularResolution(data,angleSelections=angle,doplot=False)
    st=st.astype(numpy.float)
    sp=sp.astype(numpy.float)
    energy=energy.astype(numpy.float)


    #	energy=numpy.append(energy,numpy.array([1e5,1e6]))*1e3
    #	sp=numpy.append(sp,numpy.repeat(sp[len(sp)-1],2))
    #	source_angRes=numpy.append(st[source_energy<=10],sp[source_energy>10])
    mask=numpy.isfinite(st)
    energy_comp=energy[mask]
    st=st[mask]
    mask=numpy.isfinite(sp)
    energy_pair=energy[mask]
    sp=sp[mask]

    # Loop through the events and figure out if each event is inside ARM

    f = interp1d(energy_comp,st,fill_value="extrapolate",kind='linear')

    nCompEvents=events['numberOfComptonEvents']
    CompinARM=[]
    for i in range(nCompEvents): #loop through each Compton Event in file

        interpRes=f(events['energy_ComptonEvents'][i]*1e-3)
        if ARMcut == None:
            ARMcut=interpRes
    #		print numpy.degrees(events['phi_Tracker'][i]),interpRes,events['energy_ComptonEvents'][i],i
        if abs(numpy.degrees(events['phi_Tracker'][i])-sourceTheta)<ARMcut: # is event in ARM?
            CompinARM.append(i)

    f = interp1d(energy_pair,sp,fill_value="extrapolate",kind='linear')
    angles=EventAnalysis.getARMForPairEvents(events, sourceTheta=sourceTheta, numberOfBins=100, \
        angleFitRange=[0,10], anglePlotRange=[0,180], openingAngleMax=180., showPlots=False, \
        numberOfPlots=0, finishExtraction=True, qualityCut=1, energyCut=numpy.nan, weightByEnergy=True, \
        showDiagnosticPlots=False, filename=None, log=False, getScaledDeviation=False, onlyangles=True)
    nPairEvents=events['numberOfPairEvents']
    PairinARM=[]
    for i in range(nPairEvents):
        interpRes=f((events['energy_pairElectron'][i]+events['energy_pairPositron'][i])*1e-3)
        if ARMcut == None:
            ARMcut=interpRes
        if angles[i]<ARMcut:
            PairinARM.append(i)
    #		print angles[i],interpRes,events['energy_pairElectron'][i]+events['energy_pairPositron'][i]

    return CompinARM, PairinARM, energy


##########################################################################################

def convolveAresAeff(data, angleSelections=[1,0.9,0.8,0.7,0.6,0.5], xlog=True, ylog=False, save=False, doplot=True):

    AeffData=tabulateEffectiveArea(data, angleSelections=angleSelections, ideal=True, doPrint=False, SurroundingSphere=150)

    AeffSum=[]
    Acceptance=[]
    for i in range(len(AeffData)):
        #Sum up only results for tracked and pair events
        #If all events are desired use [0:3] instead
        AeffSum.append(nansum(AeffData[i][1][0][1][1:3]))
        Acceptance.append(AeffSum[i]*2.5)

    NormAeff=AeffSum/sum(AeffSum)*19.0

    print(NormAeff)
    print(Acceptance)

    AresData=plotAngularResolution(data, angleSelections=angleSelections, doplot=False)

    return NormAeff

##########################################################################################

if __name__=='__main__':


    print(dir_path)
    if len(sys.argv) == 1:
        print('Usage: FigureOfMeritParser.py filename')

    else:

        # Get the user supplied simulation output
        simulationFilename = sys.argv[1]

        results = run(simulationFilename)

        print("")
        print(results)










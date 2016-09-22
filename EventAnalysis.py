#!/usr/bin/env python
"""
------------------------------------------------------------------------

A script to replace the energy and angular resolution measurements performed by mimrec 

Author: Daniel Kocevski (dankocevski@gmail.com)
Date: September 3rd, 2016

Usage Examples:
import EventAnalysis

# Parse the .tra file obtained from revan
events = EventAnalysis.parse('EffectiveArea_2MeV.inc1.id1.tra')
 
# Calculate the angular resolution measurement (ARM) for Compton events
FWHM_ARM, dphi = EventAnalysis.getARMForComptonEvents(events, numberOfBins=100, phiRadius=5)
 
# Calculate the angular resolution measurement (ARM) for pair events
angles, openingAngles, contaimentData_68, contaimentBinned_68 = EventAnalysis.getARMForPairEvents(events, numberOfBins=100)

# Calculate the energy resolution for Compton events
mean, FWHM = EventAnalysis.getEnergyResolutionForComptonEvents(events, numberOfBins=100, energyPlotRange=[0,10000], energyFitRange=[1800,2100])
 
# Calculate the energy resolution for Pair events
fitMax, FWHM = EventAnalysis.getEnergyResolutionForPairEvents(events, numberOfBins=100, energyPlotRange=[0,10000], energyFitRange=[1800,2100])

# Display some diagnostic plots
EventAnalysis.plotDiagnostics(events)

# Visualize the pair events
EventAnalysis.visualizePairs(events, numberOfPlots=10)

------------------------------------------------------------------------
"""

import os
import time
import sys
import fileinput
import matplotlib.pyplot as plot
from mpl_toolkits.mplot3d import Axes3D
import numpy
from itertools import product, combinations
from collections import OrderedDict
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.optimize import leastsq
import scipy.optimize
import math
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import AutoMinorLocator
import glob


# Set the default title font dict
titleFormat = {'fontsize': 12, 'fontweight' : plot.rcParams['axes.titleweight'], 'verticalalignment': 'baseline', 'horizontalalignment': 'center'}


##########################################################################################

def DoubleLorentzAsymGausArm(x, par):
	"""
	DoubleLorentzAsymGausArm(x, par)
	Fucntion Parameters:
	par[0]: Offset
	par[1]: Mean
	par[2]: Lorentz Width 1
	par[3]: Lorentz Height 1
	par[4]: Lorentz Width 2
	par[5]: Lorentz Height 2
	par[6]: Gaus Height
	par[7]: Gaus Sigma 1
	par[8]: Gaus Sigma 2
	"""

	y = par[0];

	y += numpy.abs(par[3]) * (par[2]*par[2])/(par[2]*par[2]+(x-par[1])*(x-par[1]));
	y += numpy.abs(par[5]) * (par[4]*par[4])/(par[4]*par[4]+(x-par[1])*(x-par[1]));

	for datum in x:
		arg = 0;
		if (datum - par[1] >= 0):
			if (par[7] != 0):
				arg = (datum - par[1])/par[7];
				y += numpy.abs(par[6])*numpy.exp(-0.5*arg*arg);
			else:
				if (par[8] != 0):
					arg = (datum - par[1])/par[8];
					y += numpy.abs(par[6])*numpy.exp(-0.5*arg*arg);

	return y


##########################################################################################

def DoubleLorentzian(x, offset, mean, width1, height1, width2, height2):
	"""
	DoubleLorentzAsymGausArm(x, par)
	Fucntion Parameters:
	par[0]: Offset
	par[1]: Mean
	par[2]: Lorentz Width 1
	par[3]: Lorentz Height 1
	par[4]: Lorentz Width 2
	par[5]: Lorentz Height 2
	"""

	y = offset

	y += numpy.abs(height1) * (width1*width1)/(width1*width1+(x-mean)*(x-mean));
	y += numpy.abs(height2) * (width2*width2)/(width2*width2+(x-mean)*(x-mean));

	return y

##########################################################################################

def Lorentzian(x, par):

	x = numpy.array(range(-50,50,1))/10.
	par = [0,0.5]

	mean = par[0]
	Gamma = par[1]

	y = (1/(Gamma*math.pi)) * ( ( Gamma**2 ) / ( (x-mean)**2 + Gamma**2 ) )

	plot.plot(x,y)
	plot.show()

##########################################################################################

def skewedGaussian(x,h=1, e=0,w=1,a=0):
	t = (x-e) / w
	# return 2 / w * scipy.stats.norm.pdf(t) * scipy.stats.norm.cdf(a*t)
	return h * 2 * scipy.stats.norm.pdf(t) * scipy.stats.norm.cdf(a*t)

##########################################################################################

def plotCube(shape=[80,80,60], position=[0,0,0], ax=None, color='blue'):
	"""
	A helper function to plot a 3D cube. Note that if no ax object is provided, a new plot will be created.
	Example Usage: 
	EventViewer.plotCube(shape=[50*2,50*2,35.75*2], position=[0,0,35.0], color='red', ax=ax)
	 """

	# individual axes coordinates (clockwise around (0,0,0))
	x_bottom = numpy.array([1,-1,-1,1,1]) * (shape[0]/2.) + position[0]
	y_bottom = numpy.array([-1,-1,1,1,-1]) * (shape[1]/2.) + position[1]
	z_bottom = numpy.array([-1,-1,-1,-1,-1]) * (shape[2]/2.) + position[2]

	# individual axes coordinates (clockwise around (0,0,0))
	x_top = numpy.array([1,-1,-1,1,1]) * (shape[0]/2.) + position[0]
	y_top= numpy.array([-1,-1,1,1,-1]) * (shape[1]/2.) + position[1]
	z_top = numpy.array([1,1,1,1,1]) * (shape[2]/2.) + position[2]

	x_LeftStrut = numpy.array([1,1]) * (shape[0]/2.0) + position[0]
	x_RightStrut = numpy.array([-1,-1]) * (shape[0]/2.0) + position[0]
	y_FrontStrut = numpy.array([1,1]) * (shape[1]/2.0) + position[1]
	y_BackStrut = numpy.array([-1,-1]) * (shape[1]/2.0) + position[1]
	z_Strut = numpy.array([1,-1]) * (shape[2]/2.0)  + position[2]

	if ax == None:
		fig = plot.figure()
		ax = fig.add_subplot(111, projection='3d')

	ax.plot3D(x_bottom, y_bottom, z_bottom, c=color)
	ax.plot3D(x_top, y_top, z_top, c=color)
	ax.plot3D(x_LeftStrut, y_FrontStrut, z_Strut, c=color)
	ax.plot3D(x_LeftStrut, y_BackStrut, z_Strut, c=color)
	ax.plot3D(x_RightStrut, y_FrontStrut, z_Strut, c=color)
	ax.plot3D(x_RightStrut, y_BackStrut, z_Strut, c=color)


##########################################################################################

def plotSphere(radius=300, ax=None):
	"""
	A helper function to plot a 3D sphere. Note that if no ax object is provided, a new plot will be created.
	Example Usage: 
	EventViewer.plotSphere(radius=300, ax=ax)
	 """

	#draw sphere
	u, v = numpy.mgrid[0:2*numpy.pi:20j, 0:numpy.pi:10j]
	x=numpy.cos(u)*numpy.sin(v)
	y=numpy.sin(u)*numpy.sin(v)
	z=numpy.cos(v)

	if ax == None:
		fig = plot.figure()
		ax = fig.add_subplot(111, projection='3d')

	ax.plot_wireframe(x*radius, y*radius, z*radius, color="gray")


##########################################################################################

def angularSeperation(vector1, vector2):

	# Calculate the product of the vector magnitudes
	product = numpy.linalg.norm(vector2) * numpy.linalg.norm(vector1)

	# Make sure we don't devide by zero
	if product != 0:

		# Calculate the dot product
		dotProduct = numpy.dot(vector2, vector1)

		# Make sure we have sane results
		value = dotProduct/product 
		if (value >  1.0): value =  1.0;
		if (value < -1.0): value = -1.0;

		# Get the reconstructed angle (in degrees)
		angle = numpy.degrees(numpy.arccos(value))

	else:

		# Return zero in case the denominator is zero
		angle = 0.0


	return angle


##########################################################################################


def parse(filename):

	# Create the dictionary that will contain all of the results
	events = {}

	# Define the lists to store energy information for Compton events
	energy_ComptonEvents = []
	energy_ComptonEvents_error = []	
	energy_TrackedComtonEvents = []
	energy_UntrackedComtonEvents = []
	energy_firstScatteredPhoton = []
	energy_firstScatteredPhoton_error = []
	energy_recoiledElectron = []
	energy_recoiledElectron_error = []

	# Define the lists to store energy information for pair events
	energy_pairElectron = []
	energy_pairElectron_error = []
	energy_pairPositron = []
	energy_pairPositron_error = []
	energy_pairDepositedInFirstLayer = []
	energy_pairDepositedInFirstLayer_error = []

	# Define the lists to store position information
	position_originalPhoton = []
	position_firstInteraction = []
	position_firstInteraction_error = []
	position_secondInteraction = []
	position_secondInteraction_error = []
	position_pairConversion = []
	position_firstInteraction = []
	position_firstInteraction_error = []

	# Define the lists to store the direction vectors
	direction_pairElectron = []
	direction_pairPositron = []

	# Define other lists
	phi_Tracker = []
	qualityOfComptonReconstruction = []
	qualityOfPairReconstruction = []


	# Read the number of lines in the file
	command = 'wc %s' % filename
	output = os.popen(command).read()
	totalNumberOfLines = float(output.split()[0])

	# Start the various event and line counters
	numberOfUnknownEventTypes = 0
	numberOfComptonEvents = 0
	numberOfPairEvents = 0
	numberOfPhotoElectricEffectEvents = 0
	numberOfTrackedElectronEvents = 0
	numberOfUntrackedElectronEvents = 0
	lineNumber = 0
	eventNumber = 0

	# Create empty index lists to track event types
	index_tracked = []
	index_untracked = []

	# Start by collecting all events
	skipEvent = False

	# Loop through the .tra file
	print '\nParsing: %s' % filename	
	for line in fileinput.input([filename]):

		try:
			sys.stdout.write("Progress: %d%%   \r" % (lineNumber/totalNumberOfLines * 100) )
			sys.stdout.flush()
		except:
			pass

		if 'ET ' in line:

			# Split the line
			lineContents = line.split()	

			eventType = lineContents[1]

			# Skip the event if its of unknown type
			if 'UN' in eventType:

				numberOfUnknownEventTypes = numberOfUnknownEventTypes + 1
				skipEvent = True
				continue

			else:
				skipEvent = False

			if 'CO' in eventType:
				# Increment the Compton event counter
				numberOfComptonEvents = numberOfComptonEvents + 1

			if 'PA' in eventType:
				# Increment the pair event counter
				numberOfPairEvents = numberOfPairEvents + 1

			if 'PH' in eventType:
				# Increment the photo electron effect counter
				numberOfPhotoElectricEffectEvents = numberOfPhotoElectricEffectEvents + 1


		####### Compton Events #######

		# Extract the Compton event energy information
		if 'CE ' in line and skipEvent == False:

			# Split the line
			lineContents = line.split()	

			# Get the energy of the first scattered gamma-ray
			energy_firstScatteredPhoton.append(float(lineContents[1]))
			energy_firstScatteredPhoton_error.append(float(lineContents[2]))

			# Get the energy of the recoiled electron
			energy_recoiledElectron.append(float(lineContents[3]))
			energy_recoiledElectron_error.append(float(lineContents[4]))


		# Extract the Compton event hit information
		if 'CD ' in line and skipEvent == False:

			# Split the line
			lineContents = line.split()	

			# Get the position of the first scattered gamma-ray
			x1 = float(lineContents[1])
			y1 = float(lineContents[2])
			z1 = float(lineContents[3])

			# Get the position uncertainty of the first scattered gamma-ray		
			x1_error = float(lineContents[4])
			y1_error = float(lineContents[5])
			z1_error = float(lineContents[6])

			# Get the position of the second scattered gamma-ray
			x2 = float(lineContents[7])
			y2 = float(lineContents[8])
			z2 = float(lineContents[9])

			# Get the position uncertainty of the second scattered gamma-ray		
			x2_error = float(lineContents[10])
			y2_error = float(lineContents[11])
			z2_error = float(lineContents[12])

			# Get the origin position of the original gamma-ray
			x0 = x1
			y0 = y1
			z0 = 1000.0

			# Get the position of the second scattered gamma-ray
			x_electron = float(lineContents[13])
			y_electron = float(lineContents[14])
			z_electron = float(lineContents[15])

			# Get the position uncertainty of the second scattered gamma-ray		
			x_electron_error = float(lineContents[16])
			y_electron_error = float(lineContents[17])
			z_electron_error = float(lineContents[18])

			# Record the energy of the Compton event (regardless of whether it was has a tracked or untrack electron)
			energy_ComptonEvents.append(energy_firstScatteredPhoton[-1] + energy_recoiledElectron[-1])

			# Record the energy error of the Compton event (regardless of whether it was has a tracked or untrack electron)
			energy_ComptonEvents_error.append( math.sqrt( energy_firstScatteredPhoton_error[-1]**2 + energy_recoiledElectron_error[-1]**2 ) )


			# Determine if the recoil electron was tracked by the detector
			if x_electron != 0:

				# Record the energy of tracked events
				energy_TrackedComtonEvents.append(energy_firstScatteredPhoton[-1] + energy_recoiledElectron[-1])

				# Record the number of tracked events
				numberOfTrackedElectronEvents = numberOfTrackedElectronEvents + 1

				# Add this event to the index of tracked events
				index_tracked.append(eventNumber)

			else:

				# Record the energy of tracked events
				energy_UntrackedComtonEvents.append(energy_firstScatteredPhoton[-1] + energy_recoiledElectron[-1])

				# Record the number of untracked events
				numberOfUntrackedElectronEvents = numberOfUntrackedElectronEvents + 1

				# Add this event to the index of untracked events
				index_untracked.append(eventNumber)


			# Store the origin coordinates of the original gamma-ray 
			position0 = numpy.array([x0, y0, z0])
			# position0Error = numpy.array([x1_error,y1_error,z1_error])

			# Store the coordinates of the first interaction in an array
			position1 = numpy.array([x1, y1, z1])
			position1Error = numpy.array([x1_error, y1_error, z1_error])

			# Store the coordinates of the second interaction in an array		
			position2 = numpy.array([x2, y2, z2])
			position2Error = numpy.array([x2_error, y2_error, z2_error])

			# Calculate the vector between the second and first interaction
			directionVector2 = position2 - position1

			# Calculate the vector between the first interaction and the origin of the original gamma-ray
			directionVector1 = position1 - position0

			# Calculate the product of the vector magnitudes
			product = numpy.linalg.norm(directionVector1) * numpy.linalg.norm(directionVector2)

			# Make sure we don't devide by zero
			if product != 0:

				# Calculate the dot product
				dotProduct = numpy.dot(directionVector2, directionVector1)

				# Make sure we have sane results
				value = dotProduct/product 
				if (value >  1.0): value =  1.0;
				if (value < -1.0): value = -1.0;

				# Get the reconstructed phi angle (in radians)
				phi_Tracker.append(numpy.arccos(value))

			else:

				# Return zero in case the denominator is zero
				phi_Tracker.append(0.0)

			# Add the position information to their respective list of positions
			position_originalPhoton.append(position0)
			position_firstInteraction.append(position1)
			position_firstInteraction_error.append(position1Error)
			position_secondInteraction.append(position2)
			position_secondInteraction_error.append(position2Error)

			# Increment the event number
			eventNumber = eventNumber + 1


		####### Pair Events #######

		# Extract the pair conversion hit information
		if 'PC ' in line and skipEvent == False:

			# Split the line
			lineContents = line.split()	

			# Get the position of the pair conversion
			x1 = float(lineContents[1])
			y1 = float(lineContents[2])
			z1 = float(lineContents[3])

			# Save the position
			position_pairConversion.append([x1,y1,z1])



		# Extract the pair electron information
		if 'PE ' in line and skipEvent == False:


			# Split the line
			lineContents = line.split()

			# Get the electron information
			energy_pairElectron.append(float(lineContents[1]))
			energy_pairElectron_error.append(float(lineContents[2]))

			# Get the direction of the pair electron
			x = float(lineContents[3])
			y = float(lineContents[4])
			z = float(lineContents[5])

			# Store the direction of the pair electron
			direction_pairElectron.append([x,y,z])

		# Extract the pair positron information
		if 'PP ' in line and skipEvent == False:

			# Split the line
			lineContents = line.split()

			# Get the electron information
			energy_pairPositron.append(float(lineContents[1]))
			energy_pairPositron_error.append(float(lineContents[2]))

			# Get the direction of the pair electron
			x = float(lineContents[3])
			y = float(lineContents[4])
			z = float(lineContents[5])

			# Store the direction of the pair electron
			direction_pairPositron.append([x,y,z])

		# Extract the reconstruction quality
		if 'TQ ' in line and skipEvent == False:

			# Split the line
			lineContents = line.split()

			if 'CO' in eventType:

				# Get the reconstruction quality
				qualityOfComptonReconstruction.append(float(lineContents[1]))

			if 'PA' in eventType:

				# Get the reconstruction quality
				qualityOfPairReconstruction.append(float(lineContents[1]))


		# Increment the line number for the progress indicator
		lineNumber = lineNumber + 1


	# Add all of the lists to the results dictionary
	events['energy_ComptonEvents'] = numpy.array(energy_ComptonEvents).astype(float)
	events['energy_ComptonEvents_error'] = numpy.array(energy_ComptonEvents_error).astype(float)
	events['energy_TrackedComtonEvents'] = numpy.array(energy_TrackedComtonEvents).astype(float)
	events['energy_UntrackedComtonEvents'] = numpy.array(energy_UntrackedComtonEvents).astype(float)
	events['energy_firstScatteredPhoton'] = numpy.array(energy_firstScatteredPhoton).astype(float)
	events['energy_firstScatteredPhoton_error'] = numpy.array(energy_firstScatteredPhoton_error).astype(float)
	events['energy_recoiledElectron'] = numpy.array(energy_recoiledElectron).astype(float)
	events['energy_recoiledElectron_error'] = numpy.array(energy_recoiledElectron_error).astype(float)
	events['energy_pairElectron'] = numpy.array(energy_pairElectron).astype(float)
	events['energy_pairElectron_error'] = numpy.array(energy_pairElectron_error).astype(float)
	events['energy_pairPositron'] = numpy.array(energy_pairPositron).astype(float)
	events['energy_pairPositron_error'] = numpy.array(energy_pairPositron_error).astype(float)
	events['energy_pairDepositedInFirstLayer'] = numpy.array(energy_pairDepositedInFirstLayer).astype(float)
	events['energy_pairDepositedInFirstLayer_error'] = numpy.array(energy_pairDepositedInFirstLayer_error).astype(float)

	events['position_originalPhoton'] = numpy.array(position_originalPhoton).astype(float)
	events['position_firstInteraction'] = numpy.array(position_firstInteraction).astype(float)
	events['position_firstInteraction_error'] = numpy.array(position_firstInteraction_error).astype(float)
	events['position_secondInteraction'] = numpy.array(position_secondInteraction).astype(float)
	events['position_secondInteraction_error'] = numpy.array(position_secondInteraction_error).astype(float)
	events['position_pairConversion'] = numpy.array(position_pairConversion).astype(float)

	events['direction_pairElectron'] = numpy.array(direction_pairElectron).astype(float)
	events['direction_pairPositron'] = numpy.array(direction_pairPositron).astype(float)

	events['phi_Tracker'] = numpy.array(phi_Tracker).astype(float)
	events['numberOfComptonEvents'] =numberOfComptonEvents
	events['numberOfPairEvents'] = numberOfPairEvents
	events['index_tracked'] = numpy.array(index_tracked)
	events['index_untracked']= numpy.array(index_untracked)
	events['qualityOfComptonReconstruction'] = numpy.array(qualityOfComptonReconstruction).astype(float)
	events['qualityOfPairReconstruction'] = numpy.array(qualityOfPairReconstruction).astype(float)


	# Print some event statistics
	print "\n\nStatistics of Event Selection"
	print "***********************************"
	print "Total number of analyzed events: %s" % (numberOfComptonEvents + numberOfPairEvents)

	if numberOfComptonEvents + numberOfPairEvents == 0:
		print "No events pass selection"
		return

	print ""
	print "Number of unknown events: %s (%i%%)" % (numberOfUnknownEventTypes, 100*numberOfUnknownEventTypes/(numberOfComptonEvents + numberOfPairEvents + numberOfUnknownEventTypes))
	print "Number of pair events: %s (%i%%)" % (numberOfPairEvents, 100*numberOfPairEvents/(numberOfComptonEvents + numberOfPairEvents + numberOfUnknownEventTypes))	
	print "Number of Compton events: %s (%i%%)" % (numberOfComptonEvents, 100*numberOfComptonEvents/(numberOfComptonEvents + numberOfPairEvents + numberOfUnknownEventTypes))
	print " - Number of tracked electron events: %s (%i%%)" % (numberOfTrackedElectronEvents, 100.0*(float(numberOfTrackedElectronEvents)/numberOfComptonEvents))
	print " - Number of untracked electron events: %s (%i%%)" % (numberOfUntrackedElectronEvents, 100*(float(numberOfUntrackedElectronEvents)/numberOfComptonEvents))
	print ""
	print ""


	return events


##########################################################################################

def getARMForComptonEvents(events, numberOfBins=100, phiRadius=180, includeUntrackedElectrons=True, showPlots=True):

	# Set some constants
	electron_mc2 = 511.0		# KeV

	# Retrieve the event data
	energy_firstScatteredPhoton = events['energy_firstScatteredPhoton']
	energy_recoiledElectron = events['energy_recoiledElectron']
	phi_Tracker = events['phi_Tracker']
	index_tracked = events['index_tracked']

	# Check if the user wants to exclude the untracked events
	if includeUntrackedElectrons == False:
		energy_firstScatteredPhoton = energy_firstScatteredPhoton[index_tracked]
		energy_recoiledElectron = energy_recoiledElectron[index_tracked]
		phi_Tracker = phi_Tracker[index_tracked]

	# Calculate the Compton scattering angle
	value = 1 - electron_mc2 * (1/energy_firstScatteredPhoton - 1/(energy_recoiledElectron + energy_firstScatteredPhoton));

	# Keep only sane results
	# index = numpy.where( (value > -1) | (value < 1) )
	# value = value[index]

	# Calculate the theoretical phi angle (in radians)
	phi_Theoretical = numpy.arccos(value);

	# Calculate the difference between the tracker reconstructed scattering angle and the theoretical scattering angle
	dphi = numpy.degrees(phi_Tracker) - numpy.degrees(phi_Theoretical)

	# Fit only a subsample of the data
	selection = numpy.where( (dphi > (-1*phiRadius)) & (dphi < phiRadius) )

	# Set the plot size
	plot.figure(figsize=[10,9])
	# plot.rcParams['figure.figsize'] = 10, 9

	gs = gridspec.GridSpec(4,1)
	ax1 = plot.subplot(gs[:3, :])

	# Create the histogram
	histogram_angleResults = ax1.hist(dphi[selection], numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
	ax1.set_xlim([-1*phiRadius,phiRadius])

	# Extract the binned data and bin locations
	dphi_binned = histogram_angleResults[0]
	bins = histogram_angleResults[1]
	bincenters = 0.5*(bins[1:]+bins[:-1])

	# Set the initial parameters
	offset = 0
	mean = 0
	width1 = 1
	height1 = numpy.max(dphi_binned)
	width2 = 0
	height2 = 0

	# Fit the histogram data
	optimizedParameters, covariance = scipy.optimize.curve_fit(DoubleLorentzian, bincenters, dphi_binned, [offset, mean, width1, height1, width2, height2])

	# Calculate the optimized curve
	y_fit = DoubleLorentzian(bincenters, optimizedParameters[0], optimizedParameters[1], optimizedParameters[2], optimizedParameters[3], optimizedParameters[4], optimizedParameters[5])

	# Get the fwhm of the fit
	x1 = bincenters[numpy.where(y_fit >= numpy.max(y_fit)/2)[0][0]]
	x2 = bincenters[numpy.where(y_fit >= numpy.max(y_fit)/2)[0][-1]]
	FWHM = x2-x1

	# Print some statistics
	print "\n\nStatistics of ARM histogram and fit (Compton Events)"
	print "***********************************"
	print "Compton and pair events in ARM histogram: %s (%s%%)" % ( len(dphi[selection]), 100*len(dphi[selection])/(len(dphi)) )
	print ""
	print "Mean of fit: %s" % optimizedParameters[1]	
	print "FWHM of fit: %s" %  FWHM	
	print ""

	# Annotate the plot
	ax1.text(0.03, 0.9, "Mean = %.3f deg\nFWHM = %.3f deg" % (optimizedParameters[1], FWHM), verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes, color='black', fontsize=12)

	# Plot the data
	ax1.plot(bincenters, y_fit, color='darkred', linewidth=2)
	ax1.plot([x1,x2],[numpy.max(y_fit)/2.,numpy.max(y_fit)/2.], color='darkred', linestyle='--', linewidth=2)
	ax1.set_title('Angular Resolution (Compton Events)', fontdict=titleFormat)

	# ax1.axes.xaxis.set_ticklabels([])

	# Create a subplot for the residuals 
	ax2 = plot.subplot(gs[3, :])

	# Plot the residuals
	ax2.step(bincenters, dphi_binned-y_fit, color='#3e4d8b', alpha=0.9)	
	ax2.plot([bincenters[0],bincenters[-1]], [0,0], color='darkred', linewidth=2, linestyle='--')
	ax2.set_xlim([-1*phiRadius,phiRadius])
	ax2.set_xlabel('ARM - Compton Cone (deg)')

	# Set the minor ticks
	ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
	ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
	ax2.xaxis.set_minor_locator(AutoMinorLocator(4))


	# Show the plot
	if showPlots == True:
		plot.show()
	else:
		plot.close()


	return FWHM, dphi



##########################################################################################

def getARMForPairEvents(events, numberOfBins=100, angleFitRange=[0,180], anglePlotRange=[0,45], showPlots=True, numberOfPlots=0, finishExtraction=True, qualityCut=1, energyCut=numpy.nan, wieghtByEnergy=True, showDiagnosticPlots=True):

	# Define the list to contain the resulting angle measurements
	angles = []
	openingAngles = []

	# Start some counters
	plotNumber = 0
	numberOfRejectedEvents = 0

	# Create a list to contain an index of events passing the quality cut
	index_goodQuality = []

	# Loop through each event and calculate the reconstructed photon direction and it's offset to the true direction
	for index in range(events['numberOfPairEvents']):

		# Make a cut on the quality of the pair reconstruction
		if events['qualityOfPairReconstruction'][index] > qualityCut:
			numberOfRejectedEvents = numberOfRejectedEvents + 1
			continue

		# Get the electron and positron energies
		energy_pairElectron = events['energy_pairElectron'][index]
		energy_pairPositron = events['energy_pairPositron'][index]

		# Get total energy
		energy_PairSum = energy_pairElectron + energy_pairPositron

		# Make a cut on the total summed energy of the pair
		if energy_PairSum < energyCut:
			numberOfRejectedEvents = numberOfRejectedEvents + 1
			continue

		# Register this pair as 'good'
		index_goodQuality.append(index)

		# Get the position of the gamma conversion
		position_conversion = events['position_pairConversion'][index]

		# Get the origin position of the original gamma-ray
		position_source = [position_conversion[0], position_conversion[1], 1000.0]

		# Get the electron and positron direction vectors. These are unit vectors.
		direction_electron = events['direction_pairElectron'][index]
		direction_positron = events['direction_pairPositron'][index]

		# Weight the direction vectors by the electron and positron energies
		if wieghtByEnergy == True:		
			direction_electron = direction_electron * (energy_pairElectron/energy_PairSum)
			direction_positron = direction_positron * (energy_pairPositron/energy_PairSum)

		# Get the vector that bisects the electron and positron vectors
		direction_bisect = (direction_electron + direction_positron)/2.0

		# Invert the bisect vector to obtain the reconstructed source vector
		direction_source_reconstructed = -1*direction_bisect

		# Calculate the vector between the first interaction and the origin of the original gamma-ray
		direction_source = -1*(position_conversion - position_source)

		# Calculate the product of the vector magnitudes
		angle = angularSeperation(direction_source, direction_source_reconstructed)
		angles.append(angle)

		openingAngle = angularSeperation(direction_electron, direction_positron)
		openingAngles.append(openingAngle)

		# # Calculate the product of the vector magnitudes
		# product = numpy.linalg.norm(direction_source) * numpy.linalg.norm(direction_source_reconstructed)

		# # Make sure we don't devide by zero
		# if product != 0:

		# 	# Calculate the dot product
		# 	dotProduct = numpy.dot(direction_source, direction_source_reconstructed)

		# 	# Make sure we have sane results
		# 	value = dotProduct/product 
		# 	if (value >  1.0): value =  1.0;
		# 	if (value < -1.0): value = -1.0;

		# 	# Get the reconstructed angle (in degrees)
		# 	angle = numpy.degrees(numpy.arccos(value))

		# 	# Store the angle 
		# 	angles.append(angle)

		# else:

		# 	# Return zero in case the denominator is zero
		# 	angle = 0.0
		# 	angles.append(angle)


		if plotNumber != numberOfPlots:

			print 'Plot number: %s' % index
			print ""
			print "Photon conversion coordiantes: %s, %s, %s" % (position_conversion[0], position_conversion[1], position_conversion[2])
			print ""
			print 'Source vector (True): %s, %s, %s' % (direction_source[0], direction_source[1], direction_source[2])
			print 'Source vector (Reconstructed): %s, %s, %s' % (direction_source_reconstructed[0], direction_source_reconstructed[1], direction_source_reconstructed[2])
			print 'Angle = %s' % angle

			fig = plot.figure()
			ax = fig.add_subplot(111, projection='3d')


			# Plot the geometry
			plotCube(shape=[50*2,50*2,35.75*2], position=[0,0,35.0], color='red', ax=ax)
			plotCube(shape=[40*2,40*2,30*2], position=[0,0,29.25], color='blue', ax=ax)
			plotCube(shape=[40*2,40*2,5.0*2], position=[0,0,-8.0], color='green', ax=ax)

			# Set the plot limits
			ax.set_xlim3d(-60,60)
			ax.set_ylim3d(-60,60)
			ax.set_zlim3d(-50,100)

			# Set the plot labels
			ax.set_xlabel('x')
			ax.set_ylabel('y')
			ax.set_zlabel('z')

			# Plot the electron and positron vectors							
			ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], direction_electron[0], direction_electron[1], direction_electron[2], pivot='tail', arrow_length_ratio=0.05, color='darkblue', length=50)
			ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], direction_positron[0], direction_positron[1], direction_positron[2], pivot='tail', arrow_length_ratio=0.05, color='darkred', length=50)

			# Plot the reconstructed photon direction and the true photon direction
			ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], direction_source_reconstructed[0], direction_source_reconstructed[1], direction_source_reconstructed[2], pivot='tail', arrow_length_ratio=0, color='green', linestyle='--', length=100)
			ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], direction_source[0], direction_source[1], direction_source[2], pivot='tail', arrow_length_ratio=0, color='red', linestyle='--', length=100)

			# Plot the vector bisecting the electron and positron trajectories
			ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], direction_bisect[0], direction_bisect[1], direction_bisect[2], pivot='tail', arrow_length_ratio=0, color='green', linestyle='--', length=50)
			ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], -1*direction_bisect[0], -1*direction_bisect[1], -1*direction_bisect[2], pivot='tail', arrow_length_ratio=0, color='green', linestyle='--', length=50)


			plotNumber = plotNumber + 1
			plot.show()

		elif finishExtraction == False:

			return

	if showDiagnosticPlots == True:
		plot.hist(events['energy_pairElectron'][index_goodQuality]/(events['energy_pairElectron'][index_goodQuality]+events['energy_pairPositron'][index_goodQuality]), bins=100, alpha=0.8, color='#3e4d8b', label=r'$e^{-}$')
		plot.hist(events['energy_pairPositron'][index_goodQuality]/(events['energy_pairElectron'][index_goodQuality]+events['energy_pairPositron'][index_goodQuality]), bins=100, alpha=0.8, color='darkred', label=r'$e^{+}$')
		plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper center')	
		plot.xlabel('Fractional Energy')
		ax = plot.subplot(111)
		ax.xaxis.set_minor_locator(AutoMinorLocator(4))
		ax.yaxis.set_minor_locator(AutoMinorLocator(4))
		plot.show()

		plot.hist(events['energy_pairElectron'][index_goodQuality], bins=100, alpha=0.8, color='#3e4d8b', label=r'$e^{-}$')
		plot.hist(events['energy_pairPositron'][index_goodQuality], bins=100, alpha=0.75, color='darkred', label=r'$e^{+}$')
		plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper right')	
		plot.xlabel('Energy (KeV)')
		ax = plot.subplot(111)
		ax.xaxis.set_minor_locator(AutoMinorLocator(4))
		ax.yaxis.set_minor_locator(AutoMinorLocator(4))		
		plot.show()


	# Conver the list of angles to a numpy array
	angles = numpy.array(angles)

	# # Select the events within the desired quality range
	# selection_quality = numpy.where( events['qualityOfPairReconstruction'] <= qualityCut )

	# # Apply the selection filter
	# angles_unfiltered = angles
	# angles = angles[selection_quality]

	# Select the events within the desired energy range
	selection_fit = numpy.where( (angles >= angleFitRange[0]) & (angles <= angleFitRange[1]) )

	# Apply the fit selection
	angles_fit = angles[selection_fit]

	# Set the plot size
	plot.figure(figsize=[10,7])

	# Create the axis
	ax1 = plot.subplot(111)

	# Create the histogram
	histogramResults = ax1.hist(angles_fit, bins=numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
	ax1.set_xlabel('Angular resolution (deg)')
	ax1.set_xlim(anglePlotRange)
	ax1.set_title('Angular Resolution (Pair Events)', fontdict=titleFormat)

	# Setup the minor ticks
	ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
	ax1.yaxis.set_minor_locator(AutoMinorLocator(4))

	# Convert the list of angles into a numpy array
	angles_fit = numpy.array(angles_fit)

	# Sort the angle array
	angles_fit.sort()

	# Find the 68% containment of the cumulative sum of the angle distribution
	contaimentData_68 = angles_fit[int(len(angles_fit)*.68)]

	# Extract the binned data and bin locations
	angles_binned = histogramResults[0]
	bins = histogramResults[1]
	bincenters = 0.5*(bins[1:]+bins[:-1])
	
	# Find the peak of the historgram
	angles_max = numpy.max(angles_binned)

	# Find the cumulative sum of the angle distribution and its max
	angles_binned_cumulativeSum = numpy.cumsum(angles_binned)
	angles_binned_cumulativeMax = angles_binned_cumulativeSum[-1]

	# Find the index that corresponds to 68% of the max
	index_68 = numpy.where(angles_binned_cumulativeSum >= angles_binned_cumulativeMax*0.68)[0][0]

	# Get the 68% containment of the cumulative sum of the binned angle distribution
	contaimentBinned_68 = bincenters[index_68]

	# Add the containment values to the plot
	ax1.plot([contaimentData_68,contaimentData_68], [0,ax1.get_ylim()[1]], color='green', linewidth=1.5, linestyle='--', label="68%% (data): %.2f deg" % contaimentData_68)
	ax1.plot([contaimentBinned_68,contaimentBinned_68], [0,ax1.get_ylim()[1]], color='darkred', linewidth=1.5, linestyle='--', label="68%% (histogram): %.2f deg" %  contaimentBinned_68)
	plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper right')


	# Print some statistics
	print "\n\nStatistics of ARM histogram and fit (Pair Events)"
	print "*****************************************************"
	print ""
	print "Total number of pair events: %s" % events['numberOfPairEvents']
	print "Number of pair events passing quality cut: %s (%s%%)" % ( len(angles), 100*len(angles)/(events['numberOfPairEvents']) ) 
	print "Number of pair events included in ARM histogram: %s (%s%%)" % ( len(angles_fit), 100*len(angles_fit)/(len(angles_fit)) )
	print ""
	print "Maximum: %s" % angles_max	
	print "68%% containment (data): %.2f deg" % contaimentData_68
	print "68%% containment (histogram): %.2f deg" %  contaimentBinned_68		
	print ""


	# Show the plot
	if showPlots == True:
		plot.show()
	else:
		plot.close()

	return angles, contaimentData_68, contaimentBinned_68


##########################################################################################

def getEnergyResolutionForPairEvents(events, numberOfBins=100, energyPlotRange=None, energyFitRange=[0,1e5], showPlots=True, qualityCut=1.0):

	# Retrieve the event data
	energy_pairElectron = events['energy_pairElectron']
	energy_pairPositron = events['energy_pairPositron']
	energy_pairElectron_error = events['energy_pairElectron_error']
	energy_pairPositron_error = events['energy_pairPositron_error']
	qualityOfPairReconstruction = events['qualityOfPairReconstruction']

	# Estimate the energy of the incoming photon and its associated error
	energy_pairReconstructedPhoton = energy_pairElectron + energy_pairPositron
	energy_pairReconstructedPhoton_error = numpy.sqrt((energy_pairElectron_error/energy_pairElectron)**2 + (energy_pairPositron_error/energy_pairPositron)**2) * energy_pairReconstructedPhoton

	# Select the events within the desired energy range
	selection = numpy.where( (energy_pairReconstructedPhoton >= energyFitRange[0]) & (energy_pairReconstructedPhoton <= energyFitRange[1]) & (qualityOfPairReconstruction <= qualityCut))

	# Create the histogram
	histogramResults = plot.hist(energy_pairReconstructedPhoton[selection], bins=numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
	plot.xlabel('Energy (keV)')
	plot.xlim(energyPlotRange)
	plot.close()

	# Extract the binned data and bin locations
	energy_binned = histogramResults[0]
	bins = histogramResults[1]
	bincenters = 0.5*(bins[1:]+bins[:-1])

	# Get the bin center containing the maximum of the histogram
	bin_max = bincenters[numpy.argmax(energy_binned)]

	# Set the initial parameters
	height = numpy.max(energy_binned)	# h
	scale = 1e4							# w
	location = bin_max					# epsilon
	shape = -2							# alpha

	# Fit the histogram data
	optimizedParameters, covariance = scipy.optimize.curve_fit(skewedGaussian, bincenters, energy_binned, [height, scale, location, shape])

	# Calculate the optimized curve
	try:
		y_fit = skewedGaussian(bincenters, optimizedParameters[0], optimizedParameters[1], optimizedParameters[2], optimizedParameters[3])
	except Exception, message:
	   print message



	# Get the max of the fit
	fitMax = bincenters[numpy.argmax(y_fit)]

	# Get the fwhm of the fit
	x1 = bincenters[numpy.where(y_fit >= numpy.max(y_fit)/2)[0][0]]
	x2 = bincenters[numpy.where(y_fit >= numpy.max(y_fit)/2)[0][-1]]
	FWHM = x2-x1

	# Print some statistics
	print "\n\nStatistics of energy histogram and fit (pair events)"
	print "********************************************************"
	print "Number of Compton and pair events in histogram: %s (%s%%)" % ( len(energy_pairReconstructedPhoton[selection]), 100*len(energy_pairReconstructedPhoton[selection])/(len(energy_pairReconstructedPhoton)) )
	print ""
	print "Max of fit: %s keV" % fitMax	
	print "FWHM of fit: %s keV" % FWHM	
	print ""


	# Set the plot size
	plot.figure(figsize=[10,9])

	# Setup the plot grid
	gs = gridspec.GridSpec(4,1)

	# Create the first axis
	ax1 = plot.subplot(gs[:3, :])

	# Plot the histogram
	histogramResults = plot.hist(energy_pairReconstructedPhoton[selection], bins=numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
	ax1.set_title('Energy Resolution (Pair Events)', fontdict=titleFormat)
	if energyPlotRange != None:
		ax1.set_xlim(energyPlotRange)

	# Overplot the fit
	ax1.plot(bincenters, y_fit, color='darkred', linewidth=2)
	ax1.plot([x1,x2],[numpy.max(y_fit)/2.,numpy.max(y_fit)/2.], color='darkred', linestyle='--', linewidth=2)

	# Annotate the plot
	ax1.text(0.03, 0.9, "Max = %.3f keV\nFWHM = %.3f keV" % (fitMax, FWHM), verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes, color='black', fontsize=12)

	# Create a subplot for the residuals 
	ax2 = plot.subplot(gs[3, :])

	# Plot the residuals
	ax2.step(bincenters, energy_binned-y_fit, color='#3e4d8b', alpha=0.9)	
	ax2.plot([bincenters[0],bincenters[-1]], [0,0], color='darkred', linewidth=2, linestyle='--')
	ax2.set_xlabel('Energy (keV)')
	if energyPlotRange != None:		
		ax2.set_xlim(energyPlotRange)

	ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
	ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
	ax2.xaxis.set_minor_locator(AutoMinorLocator(4))

	# Show the plot
	if showPlots == True:
		plot.show()
	else:
		plot.close()

	return fitMax, FWHM

##########################################################################################

def getEnergyResolutionForComptonEvents(events, numberOfBins=100, energyPlotRange=None, energyFitRange=[1800,2100], includeUntrackedElectrons=True, showPlots=True):

	# Retrieve the event data
	energy_ComptonEvents = events['energy_ComptonEvents']
	energy_recoiledElectron = events['energy_recoiledElectron']
	phi_Tracker = events['phi_Tracker']
	index_tracked = events['index_tracked']
	numberOfComptonEvents = events['numberOfComptonEvents']
	numberOfPairEvents = events['numberOfPairEvents']

	# Determine whether to include untracked electrons
	if includeUntrackedElectrons == False:
		energy_ComptonEvents = energy_TrackedComtonEvents

	# Convert the list to a numpy array
	energy_ComptonEvents = numpy.array(energy_ComptonEvents)

	# Select the events within the desired energy range
	if energyPlotRange != None:
		energySelection_plot = numpy.where( (energy_ComptonEvents >= energyPlotRange[0]) & (energy_ComptonEvents <= energyPlotRange[1]) )
	else:
		energySelection_plot = numpy.arange(len(energy_ComptonEvents))

	# Create the binned data
	histogram_energyResults = plot.hist(energy_ComptonEvents[energySelection_plot], numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
	plot.close()

	# Extract the binned data and bin locations
	energy_binned = histogram_energyResults[0]
	bins_energy = histogram_energyResults[1]
	bincenters_energy = 0.5*(bins_energy[1:]+bins_energy[:-1])

	# Set the range of the energy fit by finding the inflection points in the histogram. This will not work with poorly sampled data
	if energyFitRange == None:

		bin_max = numpy.argmax(energy_binned)

		for i in range(bin_max):
			# print  bin_max-i, bin_max-i-1, energy_binned[bin_max-i], energy_binned[bin_max-i-1]
			if energy_binned[bin_max-i-1] > energy_binned[bin_max-i]:
				bin_start = bincenters_energy[bin_max-i]
				break

		for i in range(len(bins_energy-1)-bin_max-2):
			# print  bin_max+i, bin_max+i+1, energy_binned[bin_max+i], energy_binned[bin_max+i+1]
			if energy_binned[bin_max+i+1] > energy_binned[bin_max+i]:
				bin_stop = bincenters_energy[bin_max+i]
				break
			else:
				bin_stop = bincenters_energy[bin_max+i]

		energyFitRange = [bin_start, bin_stop]

	# Set the range of the energy fit to a hard coded percentage of the maximum value
	# if energyFitRange == None:
		# energyFitRange = [bincenters_energy[numpy.argmax(energy_binned)]*0.91, bincenters_energy[numpy.argmax(energy_binned)]*1.15]


	# Fit a gaussian to the energy data within the user specified energy range
	energySelection_fit = numpy.where( (energy_ComptonEvents >= energyFitRange[0]) & (energy_ComptonEvents <= energyFitRange[1]) )	
	mu_Guassian, sigma = norm.fit(energy_ComptonEvents[energySelection_fit])

	# Create the Guassian fit line
	x = numpy.array(range(int(energyFitRange[0]),int(energyFitRange[1]), 1))
	y_fit = mlab.normpdf( x, mu_Guassian, sigma)

	# Calculate the Gaussian fit statistics
	FWHM_Guassian = 2*math.sqrt(2*math.log(2))*sigma


	############  Begin asymmetric fit to the binned data ############  

	# Create the binned data
	histogram_energyResults = plot.hist(energy_ComptonEvents, numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
	plot.close()	

	# Extract the binned data and bin locations
	energy_binned = histogram_energyResults[0]
	bins_energy = histogram_energyResults[1]
	bincenters = 0.5*(bins_energy[1:]+bins_energy[:-1])

	# Select only the data within the desires energy fit range 
	energySelection_fit = numpy.where( (bincenters >= energyFitRange[0]) & (bincenters <= energyFitRange[1]) )	
	energy_binned = energy_binned[energySelection_fit]
	bincenters = bincenters[energySelection_fit]

	# Get the bin center containing the maximum of the histogram
	bin_max = bincenters[numpy.argmax(energy_binned)]

	# Set the initial parameters
	height = numpy.max(energy_binned)	# h
	scale = 1e4							# w
	location = bin_max					# epsilon
	shape = -2							# alpha

	# Fit the histogram data
	optimizedParameters, covariance = scipy.optimize.curve_fit(skewedGaussian, bincenters, energy_binned, [height, scale, location, shape])

	# Calculate the optimized curve to an asymmetric gaussian
	try:
		y_fit2 = skewedGaussian(bincenters, optimizedParameters[0], optimizedParameters[1], optimizedParameters[2], optimizedParameters[3])
	except Exception, message:
	   print message

	# Get the max of the fit
	fitMax_skewedGuassian = bincenters[numpy.argmax(y_fit2)]

	# Get the fwhm of the fit
	x1 = bincenters[numpy.where(y_fit2 >= numpy.max(y_fit2)/2)[0][0]]
	x2 = bincenters[numpy.where(y_fit2 >= numpy.max(y_fit2)/2)[0][-1]]
	FWHM_skewedGuassian = x2-x1

	############  End asymmetric fit to the binned data ############  


	# Plot the histogram and the fit line, normalized to match the histogram data
	histogram_energyResults = plot.hist(energy_ComptonEvents[energySelection_plot], numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
	plot.title('Energy Resolution (Compton Events)', fontdict=titleFormat)
	plot.xlabel('Energy (keV)')

	# Set the minor ticks
	ax = plot.gca()
	ax.xaxis.set_minor_locator(AutoMinorLocator(4))
	ax.yaxis.set_minor_locator(AutoMinorLocator(4))

	# Overplot the Gaussian fit
	y_fitNormalized = y_fit/numpy.max(y_fit)*numpy.max(energy_binned)
	plot.plot(x, y_fitNormalized, color='darkred', linewidth=2)
	plot.plot([mu_Guassian-(FWHM_Guassian/2.), mu_Guassian+(FWHM_Guassian/2.)], [numpy.max(y_fitNormalized)/2.,numpy.max(y_fitNormalized)/2.], color='darkred', linestyle='dotted', linewidth=2)

	# Overplot the asymetric Gaussian fit
	plot.plot(bincenters, y_fit2, color='green', linewidth=2)
	plot.plot([x1,x2],[numpy.max(y_fit2)/2.,numpy.max(y_fit2)/2.], color='green', linestyle='dotted', linewidth=2)

	# Plot the range of the fit
	plot.plot([energyFitRange[0],energyFitRange[0]], [0,ax.get_ylim()[1]], color='lightblue', linestyle='--', linewidth=1)
	plot.plot([energyFitRange[1],energyFitRange[1]], [0,ax.get_ylim()[1]], color='lightblue', linestyle='--', linewidth=1)

	# Annotate the plot
	ax0 = plot.subplot(111)	
	ax0.text(0.03, 0.85, "Gaussian\nMean = %.3f keV\nFWHM = %.3f keV" % (mu_Guassian, FWHM_Guassian), verticalalignment='bottom', horizontalalignment='left', transform=ax0.transAxes, color='black', fontsize=12)
	ax0.text(0.03, 0.70, "Skewed Gaussian\nMax = %.3f keV\nFWHM = %.3f keV" % (fitMax_skewedGuassian, FWHM_skewedGuassian), verticalalignment='bottom', horizontalalignment='left', transform=ax0.transAxes, color='black', fontsize=12)

	# Print some statistics
	print "\n\nStatistics of energy histogram and fit (Compton events)"
	print "***********************************************************"
	print "Compton and pair events in energy histogram: %s (%s%%)" % ( len(energy_ComptonEvents[energySelection_fit]), 100*len(energy_ComptonEvents[energySelection_fit])/(numberOfComptonEvents + numberOfPairEvents))
	print ""
	print "Mean of Guassian fit: %s keV" % mu_Guassian
	print "FWHM of Guassian fit: %s keV" % FWHM_Guassian	
	print ""
	print "Max of asymmetric Guassian fit: %s keV" % fitMax_skewedGuassian
	print "FWHM of asymmetric Guassian fit: %s keV" % FWHM_skewedGuassian	

	# Show the plot
	if showPlots == True:
		plot.show()
	else:
		plot.close()

	return mu_Guassian, FWHM_Guassian, fitMax_skewedGuassian, FWHM_skewedGuassian


##########################################################################################

def plotDiagnostics(events, showPlots=True):

	# Get the angular resolution measurements
	FWHM, dphi = getARMForComptonEvents(events, showPlots=False)

	# Retrieve the event information
	energy_ComptonEvents = events['energy_ComptonEvents']
	energy_ComptonEvents_error = events['energy_ComptonEvents_error']
	energy_firstScatteredPhoton = events['energy_firstScatteredPhoton']
	energy_firstScatteredPhoton_error = events['energy_firstScatteredPhoton_error']	
	energy_recoiledElectron = events['energy_recoiledElectron']
	energy_recoiledElectron_error = events['energy_recoiledElectron_error']

	index_untracked = events['index_untracked']
	index_tracked = events['index_tracked']

	# Photon Energy vs ARM (Compton events)
	plot.figure(figsize=[11,9])
	if len(index_tracked) != 0:
		plot.scatter(dphi[index_tracked], energy_ComptonEvents[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plot.scatter(dphi[index_untracked], energy_ComptonEvents[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plot.ylabel('Reconstructed Photon Energy (keV)')
	plot.xlabel(r'ARM ($\Delta\theta$)')

	# Set the minor ticks
	ax = plot.gca()
	ax.xaxis.set_minor_locator(AutoMinorLocator(4))
	ax.yaxis.set_minor_locator(AutoMinorLocator(4))


	# Incoming photon energy error vs ARM (Compton events)
	plot.figure(figsize=[11,9])
	percentError_ComptonEvents = 100. * energy_ComptonEvents_error/energy_ComptonEvents
	if len(index_tracked) != 0:
		plot.scatter(dphi[index_tracked], percentError_ComptonEvents[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plot.scatter(dphi[index_untracked], percentError_ComptonEvents[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plot.ylabel('Reconstructed Photon Energy Error (%)')
	plot.xlabel(r'ARM ($\Delta\theta$)')

	# Set the minor ticks
	ax = plot.gca()
	ax.xaxis.set_minor_locator(AutoMinorLocator(4))
	ax.yaxis.set_minor_locator(AutoMinorLocator(4))


	# First scattered photon energy vs ARM (Compton events)
	plot.figure(figsize=[11,9])
	if len(index_tracked) != 0:
		plot.scatter(dphi[index_tracked], energy_firstScatteredPhoton[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plot.scatter(dphi[index_untracked], energy_firstScatteredPhoton[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plot.ylabel('Scattered Photon Energy (keV)')
	plot.xlabel(r'ARM ($\Delta\theta$)')

	# Set the minor ticks
	ax = plot.gca()
	ax.xaxis.set_minor_locator(AutoMinorLocator(4))
	ax.yaxis.set_minor_locator(AutoMinorLocator(4))


	# Electron recoil energy vs ARM (Compton events)
	plot.figure(figsize=[11,9])
	if len(index_tracked) != 0:
		plot.scatter(dphi[index_tracked], energy_recoiledElectron[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plot.scatter(dphi[index_untracked], energy_recoiledElectron[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plot.ylabel('Electron Recoil Energy (keV)')
	plot.xlabel(r'ARM ($\Delta\theta$)')

	# Set the minor ticks
	ax = plot.gca()
	ax.xaxis.set_minor_locator(AutoMinorLocator(4))
	ax.yaxis.set_minor_locator(AutoMinorLocator(4))


	# First scattered photon energy error vs ARM (Compton events)
	plot.figure(figsize=[11,9])
	percentError_firstScatteredPhoton = 100. * energy_firstScatteredPhoton_error/energy_firstScatteredPhoton
	if len(index_tracked) != 0:
		plot.scatter(dphi[index_tracked], percentError_firstScatteredPhoton[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plot.scatter(dphi[index_untracked], percentError_firstScatteredPhoton[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plot.ylabel('Scattered Photon Energy Error (%)')
	plot.xlabel(r'ARM ($\Delta\theta$)')

	# Set the minor ticks
	ax = plot.gca()
	ax.xaxis.set_minor_locator(AutoMinorLocator(4))
	ax.yaxis.set_minor_locator(AutoMinorLocator(4))


	# Electron recoil energy error vs ARM (Compton events)
	plot.figure(figsize=[11,9])
	percentError_recoiledElectron = 100. * energy_recoiledElectron_error/energy_recoiledElectron
	if len(index_tracked) != 0:
		plot.scatter(dphi[index_tracked], percentError_recoiledElectron[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plot.scatter(dphi[index_untracked], percentError_recoiledElectron[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plot.ylabel('Electron Recoil Energy Error (%)')
	plot.xlabel(r'ARM ($\Delta\theta$)')

	# Set the minor ticks
	ax = plot.gca()
	ax.xaxis.set_minor_locator(AutoMinorLocator(4))
	ax.yaxis.set_minor_locator(AutoMinorLocator(4))


	# Show the plot
	if showPlots == True:
		plot.show()
	else:
		plot.close()


##########################################################################################

def visualizePairs(events, numberOfPlots=10):

	getARMForPairEvents(events, numberOfPlots=numberOfPlots, finishExtraction=False, showDiagnosticPlots=False)

	return


##########################################################################################


def performCompleteAnalysis(filename=None, directory=None, energies=None, angles=None, showPlots=True, energySearchUnit='MeV', maximumComptonEnergy=10, minimumPairEnergy=10, energyRangeCompton=None, phiRadiusCompton=5):
	"""
	A function to plot the cosima output simulation file.
	Example Usage: 
	EventViewer.performCompleteAnalysis(filename='FarFieldPointSource_100MeV_Cos1.inc1.id1.tra')
	"""

	if filename == None and directory == None:
		print "*** No filename or directory provide ***"
		print "Please provide a  filename, a list of filenames, or a directory name"
		return

	# Check to see if the user supplied a directory.  If so, include all .tra files in the directory
	if directory != None:
		filenames = glob.glob(directory + '/*.tra')


	# Check if the user supplied a single file vs a list of files
	if isinstance(filename, list) == False and filename != None:
		filenames = [filename]


	# Try to get the energy from the filename
	if energies == None:
		energies = []
		for filename in filenames:
			try:
				energy = float(filename.split('_')[1].replace(energySearchUnit,''))
				energies.append(energy)

			except:
				print "*** Unable to resolve energy from filename ***"
				print "Expected filename format: FarFieldPointSource_100MeV_Cos1.inc1.id1.tra"
				return

	# Try to get the angle from the filename
	if angles == None:
		angles = []
		for filename in filenames:
			try:
				angle = float(filename.split('_')[2].split('.inc1')[0].replace('Cos',''))
				angles.append(angle)

			except:
				print "*** Unable to resolve angle from filename ***"
				print "Expected filename format: FarFieldPointSource_100MeV_Cos1.inc1.id1.tra"
				return

	# Print the identified files
	print "\nFiles identified for analysis:"
	for filename, energy, angle in zip(filenames, energies, angles):
		print "%s %s Cos %s %s" % (energy, energySearchUnit, angle, filename)
	print ""

	# Loop through the user specified filename(s) and extract the energy and angular resolution measurements
	for filename, energy, angle in zip(filenames, energies, angles):

		print "Parsing: %s %s Cos %s %s" % (energy, energySearchUnit, angle, filename)

		# Parse the .tra file obtained from revan
		events = parse(filename)

		# Don't bother measuring the energy and angular resolutuon values for Compton events above the specified maximumComptonEnergy
		if energy <= maximumComptonEnergy:

			# Calculate the energy resolution for Compton events
			print "Calculating the energy resolution for Compton events..."
			print "EventAnalysis.getEnergyResolutionForComptonEvents(events, numberOfBins=100, energyPlotRange=None, energyFitRange=%s)" % (energyRangeCompton)
			mean, FWHM_energyComptonEvents, fitMax, FWHM_skewed_energyComptonEvents = getEnergyResolutionForComptonEvents(events, numberOfBins=100, energyPlotRange=None, energyFitRange=energyRangeCompton, showPlots=showPlots)
		 
			# Calculate the angular resolution measurement (ARM) for Compton events
			print "\n\nCalculating the angular resolution measurement for Compton events..."
			print "EventAnalysis.getARMForComptonEvents(events, numberOfBins=100, phiRadius=%s)" % (phiRadiusCompton)			
			FWHM_angleComptonEvents, dphi = getARMForComptonEvents(events, numberOfBins=100, phiRadius=phiRadiusCompton, showPlots=showPlots)
		 
		else:

			mean = numpy.nan
			FWHM_energyComptonEvents = numpy.nan
			FWHM_angleComptonEvents = numpy.nan

		# Don't bother measuring the energy and angular resolutuon values for pair events below the specified minimumPairEnergy
		if energy >= minimumPairEnergy:

			# Calculate the energy resolution for Pair events
			print "\n\nCalculating the energy resolution for pair events..."
			print "EventAnalysis.getEnergyResolutionForPairEvents(events, numberOfBins=100)"
			fitMax, FWHM_pairComptonEvents = getEnergyResolutionForPairEvents(events, numberOfBins=100, showPlots=showPlots)

			# Calculate the angular resolution measurement (ARM) for pair events
			print "\n\nCalculating the angular resolution measurement for pair events..."
			print "EventAnalysis.getARMForPairEvents(events, numberOfBins=100, showDiagnosticPlots=False)"						
			angles, contaimentData_68, contaimentBinned_68 = getARMForPairEvents(events, numberOfBins=100, showDiagnosticPlots=False, showPlots=showPlots)
		 
		else:

			contaimentData_68 = numpy.nan
			FWHM_pairComptonEvents = numpy.nan

		# Open the results filename for writing
		output_filename = filename.replace('.tra','.log')
		output = open(output_filename, 'w')
		
		# Write the results to disk
		output.write("Results for simulation: %s %s Cos %s %s\n" % (energy, energySearchUnit, angle, filename))
		output.write("Compton Energy Resolution (keV): %s\n" % FWHM_energyComptonEvents)
		output.write("Compton Angular Resolution (deg): %s\n" % FWHM_angleComptonEvents)
		output.write("Pair Energy Resolution (keV): %s\n" % FWHM_pairComptonEvents)
		output.write("Pair Angular Containment (68%%): %s\n" % contaimentData_68)

		# Close the file
		output.close()

		print "Results saved to:\n%s" % output_filename


##########################################################################################

def getTriggerEfficiency(filename=None, directory=None, save=True, savefile=None):
	"""
	A function to extract the number of simulated events from a cosima .sim file
	Usage Examples: 
	EventViewer.getNumberOfSimulatedEvents(filename='FarFieldPointSource_100MeV_Cos1.inc1.id1.sim')
	EventViewer.getNumberOfSimulatedEvents(directory='./Simulations/MySimulations/')
	"""

	if filename == None and directory == None:
		print "*** No filename or directory provide ***"
		print "Please provide a  filename, a list of filenames, or a directory name"
		return

	# Check to see if the user supplied a directory.  If so, include all .tra files in the directory
	if directory != None:
		print "\nSearching: %s\n" % directory
		filenames = glob.glob(directory + '/*.sim')

	# Check if the user supplied a single file vs a list of files
	if isinstance(filename, list) == False and filename != None:
		filenames = [filename]

	# Create a dictionary to return the results
	triggerEfficiency = {} 		


	# Loop through each file
	for filename in filenames:

		print "Parsing: %s" % filename

		# Use grep to find all lines with ID and pipe the results to the tail command to read only the last line
		command = "grep 'ID ' %s | tail -n 1" % filename
		output = os.popen(command).read()

		# Extract the number of triggers
		numberOfTriggers = float(output.split()[1])

		# Extract the number of simulated events
		numberOfSimulatedEvents = float(output.split()[2])

		# Add the values to the results dictionary
		triggerEfficiency[filename] = [numberOfTriggers, numberOfSimulatedEvents]


	print "\nTrigger Efficiencies:"
	for filename in filenames:

		# Extract the values
		numberOfTriggers = triggerEfficiency[filename][0]
		numberOfSimulatedEvents = triggerEfficiency[filename][1]
	
		# Print the results
		print filename, numberOfTriggers, numberOfSimulatedEvents, ' %.2f%%' % (100*numberOfTriggers/numberOfSimulatedEvents)

	# Check to see if the results should be written to disk
	if save == True:

		# Set a default filename if none was provided
		if savefile == None:
			savefile = 'TriggerEfficiency.txt'

		# Open the file for writing
		output = open(savefile, 'w')

		# Write the results to disk
		for filename in filenames:

			# Extract the values
			numberOfTriggers = triggerEfficiency[filename][0]
			numberOfSimulatedEvents = triggerEfficiency[filename][1]

			# Write out the values
			output.write("%s %s %s\n" % (filename, numberOfTriggers, numberOfSimulatedEvents))

		# Close the file
		output.close()

		print "\nResults saved to:\n./%s\n" % savefile


	return triggerEfficiency

##########################################################################################


# if __name__ == "__main__":





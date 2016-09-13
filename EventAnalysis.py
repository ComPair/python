#!/usr/bin/env python

import os
import time
import sys
import fileinput
import matplotlib.pyplot as plt
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


"""
------------------------------------------------------------------------

Script to reconstruct events created by cosima:

Title: 
Reference: 
Link:

Author: Daniel Kocevski (dankocevski@gmail.com)
Date: September 3rd, 2016

Usage Examples:
import EventReconstruction
EventReconstruction.getGammaARM('MyComPair_Tower.inc1.id1.tra')

------------------------------------------------------------------------
"""


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

	plt.plot(x,y)
	plt.show()


##########################################################################################

class Simulation(object):

	def __init__(self):	

		self.events = []


##########################################################################################

class Event(object):

	def __init__(self):	

		self.id_trigger = None
		self.id_simulatedEvent = None
		self.time = None
		self.initialEnergy = None
		self.depositedEnergy = None
		self.escapedEnergy = None
		self.depositedEnergy_NonSensitiveMaterial = None

		self.interactions = None
		self.hits = None
		self.particleInformation = {}


##########################################################################################

class Interactions(object):

	def __init__(self):	

		self.interactionType = []
		self.ID_interaction = []
		self.ID_parentInteraction = []
		self.ID_detector = []
		self.timeStart = []
		self.x = []
		self.y = []
		self.z = []
		self.ID_parentParticleType = []
		self.x_newDirection_OriginalParticle = []
		self.y_newDirection_OriginalParticle = []
		self.z_newDirection_OriginalParticle = []
		self.x_polarization_OriginalParticle = []
		self.y_polarization_OriginalParticle = []
		self.z_polarization_OriginalParticle = []
		self.newKineticEnergy_OriginalParticle = []
		self.ID_childParticleType = []
		self.x_direction_NewParticle = []
		self.y_direction_NewParticle = []
		self.z_direction_NewParticle = []
		self.x_polarization_NewParticle = []
		self.y_polarization_NewParticle = []
		self.z_polarization_NewParticle = []
		self.newKineticEnergy_NewParticle = []


##########################################################################################

class Hits(object):

	def __init__(self):	

		self.x = []
		self.y = []
		self.z = []
		self.energy = []
		self.detector = []

        def __str__(self):

                str = ""
                hits = zip(self.detector,self.x,self.y,self.z,self.energy)
                for hit in hits:
                        str += "HTsim {:d}; {:f}; {:f}; {:f}; {:f}\n".format(*hit)
                return str



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
		fig = plt.figure()
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
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

	ax.plot_wireframe(x*radius, y*radius, z*radius, color="gray")


##########################################################################################

def parse(filename):

	# Create the dictionary that will contain all of the results
	events = {}

	# Define the lists to store data
	energy_ComptonEvents = []
	energy_ComptonEvents_error = []	
	energy_TrackedComtonEvents = []
	energy_UntrackedComtonEvents = []
	energy_firstScatteredPhoton = []
	energy_firstScatteredPhoton_error = []
	energy_recoiledElectron = []
	energy_recoiledElectron_error = []
	positions_originalPhoton = []
	positions_firstInteraction = []
	positions_firstInteraction_error = []
	positions_secondInteraction = []
	positions_secondInteraction_error = []
	phi_Tracker = []

	# Read the number of lines in the file
	command = 'wc %s' % filename
	output = os.popen(command).read()
	totalNumberOfLines = float(output.split()[0])

	# Start the various event and line counters
	numberOfComptonEvents = 0
	numberOfPairEvents = 0
	numberOfTrackedElectronEvents = 0
	numberOfUntrackedElectronEvents = 0
	lineNumber = 0
	eventNumber = 0

	# Create empty index lists to track event types
	index_tracked = []
	index_untracked = []

	# Loop through the .tra file
	print '\nParsing: %s' % filename	
	for line in fileinput.input([filename]):

		sys.stdout.write("Progress: %d%%   \r" % (lineNumber/totalNumberOfLines * 100) )
		sys.stdout.flush()

		# Extract the Compton event energy information
		if 'CE ' in line:

			# Split the line
			LineContents = line.split()	

			# Get the energy of the first scattered gamma-ray
			energy_firstScatteredPhoton.append(float(LineContents[1]))
			energy_firstScatteredPhoton_error.append(float(LineContents[2]))

			# Get the energy of the recoiled electron
			energy_recoiledElectron.append(float(LineContents[3]))
			energy_recoiledElectron_error.append(float(LineContents[4]))

			numberOfComptonEvents = numberOfComptonEvents + 1

		# Extract the Compton event hit information
		if 'CD ' in line:

			# Split the line
			LineContents = line.split()	

			# Get the position of the first scattered gamma-ray
			x1 = float(LineContents[1])
			y1 = float(LineContents[2])
			z1 = float(LineContents[3])

			# Get the position uncertainty of the first scattered gamma-ray		
			x1_error = float(LineContents[4])
			y1_error = float(LineContents[5])
			z1_error = float(LineContents[6])

			# Get the position of the second scattered gamma-ray
			x2 = float(LineContents[7])
			y2 = float(LineContents[8])
			z2 = float(LineContents[9])

			# Get the position uncertainty of the second scattered gamma-ray		
			x2_error = float(LineContents[10])
			y2_error = float(LineContents[11])
			z2_error = float(LineContents[12])

			# Get the origin position of the original gamma-ray
			x0 = x1
			y0 = y1
			z0 = 1000.0

			# Get the position of the second scattered gamma-ray
			x_electron = float(LineContents[13])
			y_electron = float(LineContents[14])
			z_electron = float(LineContents[15])

			# Get the position uncertainty of the second scattered gamma-ray		
			x_electron_error = float(LineContents[16])
			y_electron_error = float(LineContents[17])
			z_electron_error = float(LineContents[18])

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


			# Store the coordinates of the first interaction in an array
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
				dotProcution = numpy.dot(directionVector2, directionVector1)

				# Make sure we have sane results
				value = dotProcution/product 
				if (value >  1.0): value =  1.0;
				if (value < -1.0): value = -1.0;

				# Get the reconstructed phi angle (in radians)
				phi_Tracker.append(numpy.arccos(value))

			else:

				# Return zero in case the denominator is zero
				phi_Tracker.append(0.0)

			# Add the position information to their respective list of positions
			positions_originalPhoton.append(position0)
			positions_firstInteraction.append(position1)
			positions_firstInteraction_error.append(position1Error)
			positions_secondInteraction.append(position2)
			positions_secondInteraction_error.append(position2Error)

			# Increment the event number
			eventNumber = eventNumber + 1

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
	events['positions_originalPhoton'] = numpy.array(positions_originalPhoton).astype(float)
	events['positions_firstInteraction'] = numpy.array(positions_firstInteraction).astype(float)
	events['positions_firstInteraction_error'] = numpy.array(positions_firstInteraction_error).astype(float)
	events['positions_secondInteraction'] = numpy.array(positions_secondInteraction).astype(float)
	events['positions_secondInteraction_error'] = numpy.array(positions_secondInteraction_error).astype(float)
	events['phi_Tracker'] = numpy.array(phi_Tracker).astype(float)
	events['numberOfComptonEvents'] =numberOfComptonEvents
	events['numberOfPairEvents'] = numberOfPairEvents
	events['index_tracked'] = numpy.array(index_tracked)
	events['index_untracked']= numpy.array(index_untracked)


	# Print some event statistics
	print "\n\nStatistics of Event Selection"
	print "***********************************"
	print "Total number of analyzed events: %s" % (numberOfComptonEvents + numberOfPairEvents)
	print ""
	print "Number of pair events: %s (%i%%)" % (numberOfPairEvents, 100*numberOfPairEvents/(numberOfComptonEvents + numberOfPairEvents))	
	print "Number of Compton events: %s (%i%%)" % (numberOfComptonEvents, 100*numberOfComptonEvents/(numberOfComptonEvents + numberOfPairEvents))
	print " - Number of tracked electron events: %s (%i%%)" % (numberOfTrackedElectronEvents, 100.0*(float(numberOfTrackedElectronEvents)/numberOfComptonEvents))
	print " - Number of untracked electron events: %s (%i%%)" % (numberOfUntrackedElectronEvents, 100*(float(numberOfUntrackedElectronEvents)/numberOfComptonEvents))
	print ""
	print ""


	return events


##########################################################################################

def getARM(events, numberOfBins=100, phiRadius=180, includeUntrackedElectrons=True, showPlots=True):

	# Set some constants
	electron_mc2 = 511.0

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
	plt.figure(figsize=[10,9])
	# plt.rcParams['figure.figsize'] = 10, 9

	gs = gridspec.GridSpec(4,1)
	ax1 = plt.subplot(gs[:3, :])

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
	print "\n\nStatistics of ARM histogram and fit"
	print "***********************************"
	print "Compton and pair events in ARM histogram: %s (%s%%)" % ( len(dphi[selection]), 100*len(dphi[selection])/(len(dphi)) )
	print ""
	print "Mean of fit: %s" % optimizedParameters[1]	
	print "FWHM of fit: %s" %  FWHM	
	print ""

	# Annotate the plot
	ax1.text(0.03, 0.9, "Mean = %.3f\nFWHM = %.3f" % (optimizedParameters[1], FWHM), verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes, color='black', fontsize=12)

	# Plot the data
	ax1.plot(bincenters, y_fit, color='darkred', linewidth=2)
	ax1.plot([x1,x2],[numpy.max(y_fit)/2.,numpy.max(y_fit)/2.], color='darkred', linestyle='--', linewidth=2)
	# ax1.axes.xaxis.set_ticklabels([])

	# Create a subplot for the residuals 
	ax2 = plt.subplot(gs[3, :])

	# Plot the residuals
	ax2.step(bincenters, dphi_binned-y_fit, color='#3e4d8b', alpha=0.9)	
	ax2.plot([bincenters[0],bincenters[-1]], [0,0], color='darkred', linewidth=2, linestyle='--')
	ax2.set_xlim([-1*phiRadius,phiRadius])
	ax2.set_xlabel('ARM - Compton Cone')

	# Show the plot
	if showPlots == True:
		plt.show()
	else:
		plt.close()


	return FWHM, dphi


##########################################################################################

def getEnergyResolution(events, numberOfBins=100, energyPlotRange=[0,10000], energyFitRange=[1800,2100], includeUntrackedElectrons=True, showPlots=True):

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
	energySelection_plot = numpy.where( (energy_ComptonEvents >= energyPlotRange[0]) & (energy_ComptonEvents <= energyPlotRange[1]) )

	# Create the binned data
	histogram_energyResults = plt.hist(energy_ComptonEvents[energySelection_plot], numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')

	# Extract the binned data and bin locations
	energy_binned = histogram_energyResults[0]
	bins_energy = histogram_energyResults[1]
	bincenters_energy = 0.5*(bins_energy[1:]+bins_energy[:-1])

	# Fit a gaussian to the energy data within the user specified energy range
	energySelection_fit = numpy.where( (energy_ComptonEvents >= energyFitRange[0]) & (energy_ComptonEvents <= energyFitRange[1]) )	
	mu, sigma = norm.fit(energy_ComptonEvents[energySelection_fit])

	# Create a fit line
	x = numpy.array(range(energyFitRange[0],energyFitRange[1], 1))
	y_fit = mlab.normpdf( x, mu, sigma)

	# Close the plots before remaking them
	plt.close()

	# Plot the histogram and the fit line, normalized to match the histogram data
	histogram_energyResults = plt.hist(energy_ComptonEvents[energySelection_plot], numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
	plt.plot(x,(y_fit/numpy.max(y_fit))*numpy.max(energy_binned), color='darkred', linewidth=2)
	plt.xlabel('Energy (keV)')

	# Calculate the fit statistics
	FWHM = 2*math.sqrt(2*math.log(2))*sigma

	# Annotate the plot
	ax0 = plt.subplot(111)	
	ax0.text(0.03, 0.9, "Mean = %.3f\nFWHM = %.3f" % (mu, FWHM), verticalalignment='bottom', horizontalalignment='left', transform=ax0.transAxes, color='black', fontsize=12)

	# Print some statistics
	print "\n\nStatistics of energy histogram and fit"
	print "***********************************"
	print "Compton and pair events in energy histogram: %s (%s%%)" % ( len(energy_ComptonEvents[energySelection_fit]), 100*len(energy_ComptonEvents[energySelection_fit])/(numberOfComptonEvents + numberOfPairEvents))
	print ""
	print "Mean of fit: %s" % mu
	print "FWHM of fit: %s" % FWHM	
	print ""
	print ""

	# Show the plot
	if showPlots == True:
		plt.show()
	else:
		plt.close()

	return mu, FWHM

##########################################################################################

def plotDiagnostics(events, showPlots=True):

	# Get the angular resolution measurements
	FWHM, dphi = getARM(events, showPlots=False)

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
	plt.figure(figsize=[11,9])
	if len(index_tracked) != 0:
		plt.scatter(dphi[index_tracked], energy_ComptonEvents[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plt.scatter(dphi[index_untracked], energy_ComptonEvents[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plt.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plt.ylabel('Reconstructed Photon Energy (keV)')
	plt.xlabel(r'ARM ($\Delta\theta$)')


	# Incoming photon energy error vs ARM (Compton events)
	plt.figure(figsize=[11,9])
	percentError_ComptonEvents = 100. * energy_ComptonEvents_error/energy_ComptonEvents
	if len(index_tracked) != 0:
		plt.scatter(dphi[index_tracked], percentError_ComptonEvents[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plt.scatter(dphi[index_untracked], percentError_ComptonEvents[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plt.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plt.ylabel('Reconstructed Photon Energy Error (%)')
	plt.xlabel(r'ARM ($\Delta\theta$)')


	# First scattered photon energy vs ARM (Compton events)
	plt.figure(figsize=[11,9])
	if len(index_tracked) != 0:
		plt.scatter(dphi[index_tracked], energy_firstScatteredPhoton[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plt.scatter(dphi[index_untracked], energy_firstScatteredPhoton[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plt.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plt.ylabel('Scattered Photon Energy (keV)')
	plt.xlabel(r'ARM ($\Delta\theta$)')


	# Electron recoil energy vs ARM (Compton events)
	plt.figure(figsize=[11,9])
	if len(index_tracked) != 0:
		plt.scatter(dphi[index_tracked], energy_recoiledElectron[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plt.scatter(dphi[index_untracked], energy_recoiledElectron[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plt.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plt.ylabel('Electron Recoil Energy (keV)')
	plt.xlabel(r'ARM ($\Delta\theta$)')


	# First scattered photon energy error vs ARM (Compton events)
	plt.figure(figsize=[11,9])
	percentError_firstScatteredPhoton = 100. * energy_firstScatteredPhoton_error/energy_firstScatteredPhoton
	if len(index_tracked) != 0:
		plt.scatter(dphi[index_tracked], percentError_firstScatteredPhoton[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plt.scatter(dphi[index_untracked], percentError_firstScatteredPhoton[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plt.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plt.ylabel('Scattered Photon Energy Error (%)')
	plt.xlabel(r'ARM ($\Delta\theta$)')


	# Electron recoil energy error vs ARM (Compton events)
	plt.figure(figsize=[11,9])
	percentError_recoiledElectron = 100. * energy_recoiledElectron_error/energy_recoiledElectron
	if len(index_tracked) != 0:
		plt.scatter(dphi[index_tracked], percentError_recoiledElectron[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
	if len(index_untracked) != 0:
		plt.scatter(dphi[index_untracked], percentError_recoiledElectron[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
	plt.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
	plt.ylabel('Electron Recoil Energy Error (%)')
	plt.xlabel(r'ARM ($\Delta$\theta$)')


	# Show the plot
	if showPlots == True:
		plt.show()
	else:
		plt.close()


##########################################################################################

# if __name__ == "__main__":





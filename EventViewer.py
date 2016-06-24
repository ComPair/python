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


"""
------------------------------------------------------------------------

Script to visualize events created by cosima:

Title: 
Reference: 
Link:

Author: Daniel Kocevski (dankocevski@gmail.com)
Date: June 21st, 2016

Usage Examples:
import EventViewer
EventViewer.parse('MyComPair_Tower.inc1.id1.sim')

------------------------------------------------------------------------
"""

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


##########################################################################################

def PlotCube(shape=[80,80,60], position=[0,0,0], ax=None, color='blue'):

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

def PlotSphere(radius=300, ax=None):

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

# def ShowStats(filename):

	# for line in fileinput.input([filename]):


	# if 'SE' in line and currentEventNumber == 0:

	# 	# Create a new event
	# 	event = Event()

	# 	currentEventNumber = currentEventNumber + 1


	# fileinput.close()



##########################################################################################

def PlotInteractions(filename, eventNumber=1, ax=None, hidePhoto=True):

	parseLine = False
	currentEventNumber = 0
	lineNumber = 0

	interactionTypes = []
	x_interactions = []
	y_interactions = []
	z_interactions = []
	ID_interactions = []
	ID_parentInteractions = []	
	timeStarts = []
	ID_originalParticles = []

	particleInformation = {}

	interactionMarkers = {'INIT':u'$\u2193$', 'PAIR':'^', 'COMP':'s', 'PHOT': 'D', 'BREM': 'o', 'RAYL':'8', 'IONI':'d', 'INEL':'<', 'CAPT':'.', 'DECA':'h', 'ESCP':'x', 'ENTR':'+', 'EXIT':'x', 'BLAK':'-', 'ANNI':'*'}
	particleColors = {'0': 'gray', '1':'gray', '2':'red', '3':'blue', '4': 'lightgreen', '5': 'darkgreen'}
	particleLabels = {'0': r'$\gamma$', '1':r'$\gamma$', '2':r'$e^{+}$', '3':r'$e^{-}$', '4': r'$p^{+}$', '5': r'$p^{-}$'}

	for line in fileinput.input([filename]):

		lineNumber = lineNumber + 1

		# if 'IA' in line and 'PHOT' not in line:
		if 'IA' in line:

			# Skip photoelectric interactions
			if hidePhoto == True:
				if 'PHOT' in line:
					continue 

			LineContents = line.split(';')	

			interactionType = LineContents[0].split()[1].split()[0]
			ID_interaction = LineContents[0].split()[2].split()[0]
			ID_parentInteraction = LineContents[1].split()[0]
			ID_detector = LineContents[2]
			timeStart = float(LineContents[3])
			x_interaction = float(LineContents[4])
			y_interaction = float(LineContents[5])
			z_interaction = float(LineContents[6])
			ID_parentParticleType = LineContents[7].split()[0]
			x_newDirection_OriginalParticle = float(LineContents[8])
			y_newDirection_OriginalParticle = float(LineContents[9])
			z_newDirection_OriginalParticle = float(LineContents[10])
			x_polarization_OriginalParticle = LineContents[11]
			y_polarization_OriginalParticle = LineContents[12]
			z_polarization_OriginalParticle = LineContents[13]
			newKineticEnergy_OriginalParticle = LineContents[14]
			ID_childParticleType = LineContents[15]
			x_direction_NewParticle = float(LineContents[16])
			y_direction_NewParticle = float(LineContents[17])
			z_direction_NewParticle = float(LineContents[18])
			x_polarization_NewParticle = LineContents[19]
			y_polarization_NewParticle = LineContents[20]
			z_polarization_NewParticle = LineContents[21]
			newKineticEnergy_NewParticle = LineContents[22].rstrip()

			if 'INIT' in line:
				event.initialEnergy = newKineticEnergy_NewParticle

			# Create a unique particle id to track parent and child particles
			ID_parentParticle = ID_parentInteraction + '_' + ID_parentParticleType
			ID_childParticle = ID_interaction + '_' + ID_childParticleType

			# if ID_childParticleType == '1':
			# 	ID_childParticle = ID_parentInteraction + '_' + ID_childParticleType
			# else:
			# 	ID_childParticle = ID_interaction + '_' + ID_childParticleType

			# Store the information for the individual particles associated with this interaction
			if ID_parentParticle in particleInformation:
				particleInformation[ID_parentParticle]['x'].append(x_interaction)
				particleInformation[ID_parentParticle]['y'].append(y_interaction)
				particleInformation[ID_parentParticle]['z'].append(z_interaction)
				particleInformation[ID_parentParticle]['time'].append(timeStart)
			else:
				particleInformation[ID_parentParticle] = {}
				particleInformation[ID_parentParticle]['x'] = [x_interaction]
				particleInformation[ID_parentParticle]['y'] = [y_interaction]
				particleInformation[ID_parentParticle]['z'] = [z_interaction]
				particleInformation[ID_parentParticle]['time'] = [timeStart]

			if ID_childParticle in particleInformation:
				particleInformation[ID_childParticle]['x'].append(x_interaction)
				particleInformation[ID_childParticle]['y'].append(y_interaction)
				particleInformation[ID_childParticle]['z'].append(z_interaction)
				particleInformation[ID_childParticle]['time'].append(timeStart)
			else:
				particleInformation[ID_childParticle] = {}
				particleInformation[ID_childParticle]['x'] = [x_interaction]
				particleInformation[ID_childParticle]['y'] = [y_interaction]
				particleInformation[ID_childParticle]['z'] = [z_interaction]
				particleInformation[ID_childParticle]['time'] = [timeStart]

			# Add the current interaction information to the relevant arrays
			interactionTypes.append(interactionType)
			x_interactions.append(x_interaction)
			y_interactions.append(y_interaction)
			z_interactions.append(z_interaction)
			timeStarts.append(timeStart)
			ID_interactions.append(ID_interaction)

		if 'ID' in line and currentEventNumber != 0:

			event.id_trigger = line.split()[1]
			event.id_simulatedEvent = line.split()[2]

		if 'TI' in line and currentEventNumber != 0:

			event.time = line.split()[1]

		if 'ED' in line and currentEventNumber != 0:

			event.depositedEnergy = line.split()[1]

		if 'EC' in line and currentEventNumber != 0:

			event.escapedEnergy = line.split()[1]

		if 'NS' in line and currentEventNumber != 0:

			event.depositedEnergy_NonSensitiveMaterial = line.split()[1]			

		if 'SE' in line and currentEventNumber != 0:

			if eventNumber == currentEventNumber or eventNumber == 'All':

				if ax == None:
					fig = plt.figure(figsize=(8,5))
					# ax = fig.add_subplot(111, projection='3d', aspect='equal')
					ax = fig.add_subplot(111, projection='3d')

				ax.set_xlim3d(-60,60)
				ax.set_ylim3d(-60,60)
				ax.set_zlim3d(-50,100)

				ax.set_xlabel('x')
				ax.set_ylabel('y')
				ax.set_zlabel('z')

				# ax.axis('off')

				# Generate the color map
				norm = plt.Normalize()
				timeStarts[0] = timeStarts[1]
				colors = plt.cm.jet(norm(numpy.log(timeStarts)))

				for x, y, z, interactionType, ID, color in zip(x_interactions, y_interactions, z_interactions, interactionTypes, ID_interactions, colors):
					
					if interactionType == 'INIT':
						size = 200
						color='gray'
					elif interactionType == 'ANNI':
						size = 75
					else:
						size = 25

					ax.scatter(x, y, z, marker=interactionMarkers[interactionType], s=size, c=color, label=interactionType)
					ax.text(x, y, z, '  %s' % (ID), size=10, zorder=1)


				for ID in particleInformation.keys():

					time = numpy.array(particleInformation[ID]['time'])
					x = numpy.array(particleInformation[ID]['x'])
					y = numpy.array(particleInformation[ID]['y'])
					z = numpy.array(particleInformation[ID]['z'])
					particleType = ID.split('_')[1]

					i_timeSorted = numpy.argsort(time)
					x_sorted = x[i_timeSorted]
					y_sorted = y[i_timeSorted]
					z_sorted = z[i_timeSorted]

					ax.plot(x_sorted, y_sorted, z_sorted, c=particleColors[particleType], label=particleLabels[particleType], linewidth=0.5)


				# Plot the legend
				handles, labels = ax.get_legend_handles_labels()
				by_label = OrderedDict(zip(labels, handles))
				ax.legend(by_label.values(), by_label.keys(), scatterpoints=1)

				# Plot the geometry
				PlotCube(shape=[50*2,50*2,35.75*2], position=[0,0,35.0], color='red', ax=ax)
				PlotCube(shape=[40*2,40*2,30*2], position=[0,0,29.25], color='blue', ax=ax)
				PlotCube(shape=[40*2,40*2,5.0*2], position=[0,0,-8.0], color='green', ax=ax)

				print '\nShowing Trigger: %s' % event.id_trigger
				print ' - Time: %s' % event.time				
				print ' - Simulated event: %s' % event.id_simulatedEvent
				print ' - Initial energy: %s keV' % event.initialEnergy
				print ' - Total deposited energy in sensitive material: %s keV (%.2f%%)' % (event.depositedEnergy, float(event.depositedEnergy)/float(event.initialEnergy) * 100)
				print ' - Total deposited energy in non-sensitive material: %s keV (%.2f%%)' % (event.depositedEnergy_NonSensitiveMaterial, float(event.depositedEnergy_NonSensitiveMaterial)/float(event.initialEnergy) * 100)
				print ' - Total escaped energy: %s keV (%.2f%%)' % (event.escapedEnergy, float(event.escapedEnergy)/float(event.initialEnergy) * 100)

				# Show the plot
				plt.show()


			interactionTypes = []
			x_interactions = []
			y_interactions = []
			z_interactions = []
			ID_interactions = []
			ID_parentInteractions = []	
			timeStarts = []
			ID_originalParticles = []

			particleInformation = {}

			# Increment the event number
			currentEventNumber = currentEventNumber + 1

			# Create a new event
			event = Event()


		if 'SE' in line and currentEventNumber == 0:

			# Create a new event
			event = Event()

			currentEventNumber = currentEventNumber + 1


	fileinput.close()


##########################################################################################

def PlotHits(filename, interaction=1, returnAxis=False, marker='o'):

	# Create the neccessary lists
	x_hit = []
	y_hit = []
	z_hit = []
	energy = []

	eventNumber = 0

	for line in fileinput.input([filename]):

		if 'HTsim' in line:

			# Split the line
			LineContents = line.split(';')	

			# Extract the hit information
			x_hit.append(float(LineContents[1]))
			y_hit.append(float(LineContents[2]))
			z_hit.append(float(LineContents[3]))
			energy.append(float(LineContents[4]))

		if 'SE' in line and eventNumber != 0:

			fig = plt.figure()
			ax = fig.add_subplot(111, projection='3d', aspect='equal')

			ax.set_xlim3d(-60,60)
			ax.set_ylim3d(-60,60)
			ax.set_zlim3d(-50,100)

			ax.set_xlabel('x')
			ax.set_ylabel('y')
			ax.set_zlabel('z')

			# ax.axis('off')

			# Plot the hits
			ax.scatter(x_hit, y_hit, z_hit, c=numpy.log(energy), marker=marker)

			# Plot the geometry
			PlotCube(shape=[50*2,50*2,35.75*2], position=[0,0,35.0], color='red', ax=ax)
			PlotCube(shape=[40*2,40*2,30*2], position=[0,0,29.25], color='blue', ax=ax)
			PlotCube(shape=[40*2,40*2,5.0*2], position=[0,0,-8.0], color='green', ax=ax)


			eventNumber = eventNumber + 1

			x_hit = []
			y_hit = []
			z_hit = []
			energy = []

			if returnAxis == True:
				fileinput.close()
				return ax
			else:
				plt.show()

		if 'SE' in line and eventNumber == 0:
			eventNumber = eventNumber + 1


##########################################################################################

def Plot(filename):

	ax = PlotHits(filename, returnAxis=True, marker='x')
	PlotInteractions(filename, ax=ax)



##########################################################################################


if __name__ == "__main__":

    # Read in the file name
    if (len(sys.argv) == 2):
        filename = sys.argv[1]

        parseEventFile(filename)



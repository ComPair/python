#!/usr/bin/env python
"""
------------------------------------------------------------------------

Script to visualize events created by cosima:

Author: Daniel Kocevski (dankocevski@gmail.com)
Date: June 21st, 2016

Usage Examples:
import EventViewer
EventViewer.plot('MyComPair_Tower.inc1.id1.sim')

------------------------------------------------------------------------
"""

import os
import time
import sys
import fileinput
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
from itertools import product, combinations
from collections import OrderedDict




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

def parse(filename, skipPhoto=True):
	"""
	A function to parse the cosima output simulation file.  The function returns a simulation object.
	Example Usage: 
	simulation = EventViewer.parse(filename)
	 """

	# Start an event counter
	currentEventNumber = 0

	# Loop through each line of the file
	for line in fileinput.input([filename]):

		# Create the first event
		if 'SE' in line and currentEventNumber == 0:

			# Create a new simulation object to store all of the events in this run
			simulation = Simulation()

			# Create a new event
			event = Event()

			# Create a new object to store the interactions for this event
			interactions = Interactions()

			# Create a new object to store the hits for this event
			hits = Hits()

			# Increment the event number
			currentEventNumber = currentEventNumber + 1

		# Store the existing event and create a new event
		elif 'SE' in line or 'EN' in line:
			
			# Store the interaction and hit objects in their respective event
			event.interactions = interactions
			event.hits = hits

			# Store the current event in the simulation object			
			simulation.events.append(event)

			# Create a new event
			event = Event()

			# Create a new object to store the interactions for the new event
			interactions = Interactions()

			# Create a new object to store the hits for the new event
			hits = Hits()

			# Increment the event number
			currentEventNumber = currentEventNumber + 1

		# Get the event ID
		if 'ID' in line and currentEventNumber != 0:

			event.id_trigger = line.split()[1]
			event.id_simulatedEvent = line.split()[2]

		# Get the event time
		if 'TI' in line and currentEventNumber != 0:

			event.time = line.split()[1]

		# Get the total deposited energy 
		if 'ED' in line and currentEventNumber != 0:

			event.depositedEnergy = line.split()[1]

		# Get the total escaped energy 
		if 'EC' in line and currentEventNumber != 0:

			event.escapedEnergy = line.split()[1]

		# Get the total deposited energy in non-sensative material
		if 'NS' in line and currentEventNumber != 0:

			event.depositedEnergy_NonSensitiveMaterial = line.split()[1]		

		# if 'IA' in line and 'PHOT' not in line:
		if 'IA' in line:

			# Skip photoelectric interactions
			if skipPhoto == True:
				if 'PHOT' in line:
					continue 

			# Split the line
			LineContents = line.split(';')	

			# Parse each line and place the extracted information into their respective arrays
			interactions.interactionType.append(LineContents[0].split()[1].split()[0])
			interactions.ID_interaction.append(LineContents[0].split()[2].split()[0])
			interactions.ID_parentInteraction.append(LineContents[1].split()[0])
			interactions.ID_detector.append(LineContents[2])
			interactions.timeStart.append(float(LineContents[3]))
			interactions.x.append(float(LineContents[4]))
			interactions.y.append(float(LineContents[5]))
			interactions.z.append(float(LineContents[6]))
			interactions.ID_parentParticleType.append(LineContents[7].split()[0])
			interactions.x_newDirection_OriginalParticle.append(float(LineContents[8]))
			interactions.y_newDirection_OriginalParticle.append(float(LineContents[9]))
			interactions.z_newDirection_OriginalParticle.append(float(LineContents[10]))
			interactions.x_polarization_OriginalParticle.append(LineContents[11])
			interactions.y_polarization_OriginalParticle.append(LineContents[12])
			interactions.z_polarization_OriginalParticle.append(LineContents[13])
			interactions.newKineticEnergy_OriginalParticle.append(LineContents[14])
			interactions.ID_childParticleType.append(LineContents[15])
			interactions.x_direction_NewParticle.append(float(LineContents[16]))
			interactions.y_direction_NewParticle.append(float(LineContents[17]))
			interactions.z_direction_NewParticle.append(float(LineContents[18]))
			interactions.x_polarization_NewParticle.append(LineContents[19])
			interactions.y_polarization_NewParticle.append(LineContents[20])
			interactions.z_polarization_NewParticle.append(LineContents[21])
			interactions.newKineticEnergy_NewParticle.append(LineContents[22].rstrip())

			if 'INIT' in line:
				event.initialEnergy = interactions.newKineticEnergy_NewParticle[-1]

			# Create a unique particle id to track parent and child particles
			ID_parentParticle = interactions.ID_parentInteraction[-1] + '_' + interactions.ID_parentParticleType[-1]
			ID_childParticle = interactions.ID_interaction[-1] + '_' + interactions.ID_childParticleType[-1]

			# if ID_childParticleType == '1':
			# 	ID_childParticle = ID_parentInteraction + '_' + ID_childParticleType
			# else:
			# 	ID_childParticle = ID_interaction + '_' + ID_childParticleType

			# Store the information for the individual particles associated with this interaction

			# Record the particle trajectory
			if ID_parentParticle in event.particleInformation:
				event.particleInformation[ID_parentParticle]['x'].append(interactions.x[-1])
				event.particleInformation[ID_parentParticle]['y'].append(interactions.y[-1])
				event.particleInformation[ID_parentParticle]['z'].append(interactions.z[-1])
				event.particleInformation[ID_parentParticle]['time'].append(interactions.timeStart[-1])
			else:
				event.particleInformation[ID_parentParticle] = {}
				event.particleInformation[ID_parentParticle]['x'] = [interactions.x[-1]]
				event.particleInformation[ID_parentParticle]['y'] = [interactions.y[-1]]
				event.particleInformation[ID_parentParticle]['z'] = [interactions.z[-1]]
				event.particleInformation[ID_parentParticle]['time'] = [interactions.timeStart[-1]]

			if ID_childParticle in event.particleInformation:
				event.particleInformation[ID_childParticle]['x'].append(interactions.x[-1])
				event.particleInformation[ID_childParticle]['y'].append(interactions.y[-1])
				event.particleInformation[ID_childParticle]['z'].append(interactions.z[-1])
				event.particleInformation[ID_childParticle]['time'].append(timeStart[-1])
			else:
				event.particleInformation[ID_childParticle] = {}
				event.particleInformation[ID_childParticle]['x'] = [interactions.x[-1]]
				event.particleInformation[ID_childParticle]['y'] = [interactions.y[-1]]
				event.particleInformation[ID_childParticle]['z'] = [interactions.z[-1]]
				event.particleInformation[ID_childParticle]['time'] = [interactions.timeStart[-1]]
	
		# Record the hit information
		if 'HTsim' in line:

			# Split the line
			LineContents = line.split(';')	
                        
			# Extract the hit information
                        hits.detector.append(int(LineContents[0].split(' ')[1]))
			hits.x.append(float(LineContents[1]))
			hits.y.append(float(LineContents[2]))
			hits.z.append(float(LineContents[3]))
			hits.energy.append(float(LineContents[4]))


	# Close the input file
	fileinput.close()

	return simulation


##########################################################################################

def plot(filename, showEvent=1, ax=None, hidePhoto=True, showInteractions=True, showHits=True, hitMarker='o', geometry=None):
	"""
	A function to plot the cosima output simulation file.
	Example Usage: 
	EventViewer.plot(filename)
	 """

	simulation = parse(filename, skipPhoto=True)

	# Create a dictionary of symbols and colors for the various interaction and particle types
	interactionMarkers = {'INIT':u'$\u2193$', 'PAIR':'^', 'COMP':'s', 'PHOT': 'D', 'BREM': 'o', 'RAYL':'8', 'IONI':'d', 'INEL':'<', 'CAPT':'.', 'DECA':'h', 'ESCP':'x', 'ENTR':'+', 'EXIT':'x', 'BLAK':'-', 'ANNI':'*'}
	particleColors = {'0': 'gray', '1':'gray', '2':'red', '3':'blue', '4': 'lightgreen', '5': 'darkgreen'}
	particleLabels = {'0': r'$\gamma$', '1':r'$\gamma$', '2':r'$e^{+}$', '3':r'$e^{-}$', '4': r'$p^{+}$', '5': r'$p^{-}$'}

	# Loop through all of the events
	for event in simulation.events:

		# Check and see if the user wants to display this event
		if showEvent == int(event.id_trigger) or showEvent == 'All':

			# Check to see if the interactions plot should be displayed
			if showInteractions == True:

				# Create a new 3D axis object
				fig = plt.figure(figsize=(15,12))
				fig.suptitle('Event Interactions')
				ax = fig.add_subplot(111, projection='3d')

				# Let's make sure a zero doesn't make it into the time array so that we can display it on a log scale
				event.interactions.timeStart[0] = event.interactions.timeStart[1]	

				# Generate the color map
				norm = plt.Normalize()
				colors = plt.cm.jet(norm(numpy.log(event.interactions.timeStart)))

				# Plot each individual interaction 
				for x, y, z, interactionType, ID, color in zip(event.interactions.x, event.interactions.y, event.interactions.z, event.interactions.interactionType, event.interactions.ID_interaction, colors):
					
					# Change the size and color of the plotting symbol based on the interaction
					if interactionType == 'INIT':
						size = 200
						color='gray'
					elif interactionType == 'ANNI':
						size = 75
					else:
						size = 25

					# Plot the interaction
					ax.scatter(x, y, z, marker=interactionMarkers[interactionType], s=size, c=color, label=interactionType)

					# Label the interaction
					ax.text(x, y, z, '  %s' % (ID), size=10, zorder=1)

				# Add a colorbar
				interactionPlot = ax.scatter(event.interactions.x, event.interactions.y, event.interactions.z, marker='x', c=numpy.log(event.interactions.timeStart), s=1)				
				colorbar = fig.colorbar(interactionPlot, shrink=1, aspect=20)
				colorbar.set_label('log Time (sec)')

				# Plot the particle trajectories
				for ID in event.particleInformation.keys():

					time = numpy.array(event.particleInformation[ID]['time'])
					x = numpy.array(event.particleInformation[ID]['x'])
					y = numpy.array(event.particleInformation[ID]['y'])
					z = numpy.array(event.particleInformation[ID]['z'])
					particleType = ID.split('_')[1]

					# Make sure the trajectories are sorted intime
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
				if geometry == None:

					print "\nWarning: No geometry file specified.\n"

				# Define a named geometry for Compair
				elif 'Compair' in geometry and '/' not in geometry:

					# Plot the default compair geometry
					plotCube(shape=[50*2,50*2,35.75*2], position=[0,0,35.0], color='red', ax=ax)
					plotCube(shape=[40*2,40*2,30*2], position=[0,0,29.25], color='blue', ax=ax)
					plotCube(shape=[40*2,40*2,5.0*2], position=[0,0,-8.0], color='green', ax=ax)

					# Set the plot limits
					ax.set_xlim3d(-60,60)
					ax.set_ylim3d(-60,60)
					ax.set_zlim3d(-50,100)

				elif 'Amego' in geometry and '/' not in geometry:

					# Plot the default compair geometry
					plotCube(shape=numpy.array([52.5, 52.5, 0.75])*2, position=[0.0, 0.0, 70.25], color='red', ax=ax)		# ACD
					plotCube(shape=numpy.array([40.0, 40.0, 30.5])*2, position=[0.0, 0.0, 30.5], color='blue', ax=ax)		# Tracker
					plotCube(shape=numpy.array([40.0, 40.0, 2.0])*2, position=[0.0, 0.0, -3.0], color='green', ax=ax)		# CalCZT
					plotCube(shape=numpy.array([40.0, 40.0, 6.5])*2, position=[0.0, 0.0, -12.0], color='green', ax=ax)		# CalCSI

					# Set the plot limits
					ax.set_xlim3d(-60,60)
					ax.set_ylim3d(-60,60)
					ax.set_zlim3d(-50,100)

				# Read the geometry from a file
				else:

					# Create a dictionary to contain the found volumes
					volumes = {}
					name = 'none'

					# Loop through each line of the geometry file and extract the properties of each volume
					for line in fileinput.input([geometry]):

						# Get the volume name
						if 'Volume ' in line:
							lineContents = line.split()							
							name = lineContents[1]
							# print "\n%s" % name
							volumes[name] = []

						# Get the volume dimensions
						if (name + '.Shape') in line:

							if 'World' in name:
								position = [0,0,0]
								# print position
								volumes[name].append(position)

							lineContents = line.split()
							dimensions = [float(lineContents[2]), float(lineContents[3]), float(lineContents[4])]
							# print dimensions

							# Add the dimensions of the volume to the dictionary
							volumes[name].append(dimensions)


						# Get the volume position
						if (name + '.Position') in line:
							lineContents = line.split()
							position = [float(lineContents[1]), float(lineContents[2]), float(lineContents[3])]
							# print position
							volumes[name].append(position)

					# Plot each of the volumes
					for key in volumes:

						# Skip the world volume
						if 'World' in key:
							continue
			
						if len(volumes[key]) == 2:
							print "\nPlotting: %s" % key
							print " - Dimensions: %s" % volumes[key][0]
							print " - Position: %s" % volumes[key][1]
							plotCube(shape=numpy.array(volumes[key][0])*2, position=volumes[key][1], ax=ax)
						else:
							print "\nSkipping: %s" % key


				# Set the plot labels
				ax.set_xlabel('x')
				ax.set_ylabel('y')
				ax.set_zlabel('z')

				plt.show()


			# Check to see if the hits plot should be displayed
			if showHits == False:

				plt.show()

			else:

				# Create a new 3D axis object
				fig = plt.figure(figsize=(15,12))
				fig.suptitle('Detector Hits')				
				ax = fig.add_subplot(111, projection='3d')				

				# Plot the hits
				hitPlot = ax.scatter(event.hits.x, event.hits.y, event.hits.z, c=numpy.log(event.hits.energy), marker=hitMarker, label='Detector Hit')

				# Add a colorbar
				colorbar = fig.colorbar(hitPlot, shrink=1, aspect=20)
				colorbar.set_label('log Energy (keV)')

				# Plot the legend
				handles, labels = ax.get_legend_handles_labels()
				by_label = OrderedDict(zip(labels, handles))
				ax.legend(by_label.values(), by_label.keys(), scatterpoints=1)


				# Plot the geometry
				if geometry == None:

					print "\nWarning: No geometry file specified.\n"

				# Define a named geometry for Compair
				elif 'Compair' in geometry and '/' not in geometry:

					# Plot the default compair geometry
					plotCube(shape=[50*2,50*2,35.75*2], position=[0,0,35.0], color='red', ax=ax)
					plotCube(shape=[40*2,40*2,30*2], position=[0,0,29.25], color='blue', ax=ax)
					plotCube(shape=[40*2,40*2,5.0*2], position=[0,0,-8.0], color='green', ax=ax)

					# Set the plot limits
					ax.set_xlim3d(-60,60)
					ax.set_ylim3d(-60,60)
					ax.set_zlim3d(-50,100)

				elif 'Amego' in geometry and '/' not in geometry:

					# Plot the default compair geometry
					plotCube(shape=numpy.array([52.5, 52.5, 0.75])*2, position=[0.0, 0.0, 70.25], color='red', ax=ax)		# ACD
					plotCube(shape=numpy.array([40.0, 40.0, 30.5])*2, position=[0.0, 0.0, 30.5], color='blue', ax=ax)		# Tracker
					plotCube(shape=numpy.array([40.0, 40.0, 2.0])*2, position=[0.0, 0.0, -3.0], color='green', ax=ax)		# CalCZT
					plotCube(shape=numpy.array([40.0, 40.0, 6.5])*2, position=[0.0, 0.0, -12.0], color='green', ax=ax)		# CalCSI

					# Set the plot limits
					ax.set_xlim3d(-60,60)
					ax.set_ylim3d(-60,60)
					ax.set_zlim3d(-50,100)


				# Read the geometry from a file
				else:

					# Create a dictionary to contain the found volumes
					volumes = {}
					name = 'none'

					# Loop through each line of the geometry file and extract the properties of each volume
					for line in fileinput.input([geometry]):

						# Get the volume name
						if 'Volume ' in line:
							lineContents = line.split()							
							name = lineContents[1]
							print "\n%s" % name
							volumes[name] = []

						# Get the volume dimensions
						if (name + '.Shape') in line:

							if 'World' in name:
								position = [0,0,0]
								print position
								volumes[name].append(position)

							lineContents = line.split()
							dimensions = [float(lineContents[2]), float(lineContents[3]), float(lineContents[4])]
							print dimensions

							# Add the dimensions of the volume to the dictionary
							volumes[name].append(dimensions)


						# Get the volume position
						if (name + '.Position') in line:
							lineContents = line.split()
							position = [float(lineContents[1]), float(lineContents[2]), float(lineContents[3])]
							print position
							volumes[name].append(position)

					# Plot each of the volumes
					for key in volumes:

						# Skip the world volume
						if 'World' in key:
							continue
			
						if len(volumes[key]) == 2:

							print "\nPlotting: %s" % key
							print " - Dimensions: %s" % volumes[key][0]
							print " - Position: %s" % volumes[key][1]
							plotCube(shape=volumes[key][0], position=volumes[key][1], ax=ax)

						else:
							print "\nSkipping: %s" % key

				# Display the event demographics
				print '\nShowing Trigger: %s' % event.id_trigger
				print ' - Time: %s sec' % event.time				
				print ' - Simulated event: %s' % event.id_simulatedEvent
				print ' - Initial energy: %s keV' % event.initialEnergy
				print ' - Total deposited energy in sensitive material: %s keV (%.2f%%)' % (event.depositedEnergy, float(event.depositedEnergy)/float(event.initialEnergy) * 100)
				print ' - Total deposited energy in non-sensitive material: %s keV (%.2f%%)' % (event.depositedEnergy_NonSensitiveMaterial, float(event.depositedEnergy_NonSensitiveMaterial)/float(event.initialEnergy) * 100)
				print ' - Total escaped energy: %s keV (%.2f%%)' % (event.escapedEnergy, float(event.escapedEnergy)/float(event.initialEnergy) * 100)

				# Set the plot limits
				ax.set_xlim3d(-60,60)
				ax.set_ylim3d(-60,60)
				ax.set_zlim3d(-50,100)

				# Set the plot labels
				ax.set_xlabel('x')
				ax.set_ylabel('y')
				ax.set_zlabel('z')

				# Show the plot
				plt.show()


##########################################################################################

# if __name__ == "__main__":





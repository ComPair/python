#!/usr/bin/env python
"""
------------------------------------------------------------------------

A Wrapper script to run the entire MEGAlib simulation and Compair analysis pipeline 

Author: Daniel Kocevski (dankocevski@gmail.com)
Date: September 23rd, 2016

Usage Examples:
import EventViewer
SimulatonWrapper.Run('file.txt')

------------------------------------------------------------------------
"""

import glob
import numpy
import os
import EventAnalysis


def run(filename=None, directory=None, revanConfigFile=None, seed=None):

	if filename == None and directory == None:
		print "*** No filename or directory provide ***"
		print "Please provide a  filename, a list of filenames, or a directory name"
		return


	# Check to see if the user supplied a directory.  If so, get the list of source files
	if directory != None:
		sourcefiles = glob.glob(directory + '/*.source')

	# Check if the user supplied a single file vs a list of files
	if isinstance(filename, list) == False and filename != None:
		sourcefiles = [filename]


	# Get the revan config file if one was not provided
	revanConfigFile = glob.glob(directory + '/*.cfg')[0]

	# Generate a list of seeds
	seeds = int(numpy.random(len(sourcefiles))*100)

	# Loop through each of the source files and run them through cosima and revan, and analyze their output with EventAnalysis
	for sourcefile, seed in zip(sourcefiles, seeds):

		# Generate the cosima command
		command_cosima = "cosima -s %s %s" % (seed, sourcefile)

		# Issued the cosima command
		print command_cosima
		output = os.system(command_cosima)


		# Generate the sim filename
		simfile = sourcefile.replace('.source','.sim')

		# Generate the revan command
		command_revan = "revan -f %s -c %s" % (simfile, revanConfigFile)

		# Issued the revan command
		print command_revan
		output = os.system(command_revan)


		# Extract the number of triggered and simulated events
		EventAnalysis.getTriggerEfficiency(filename=simfile)

		# Generate the .tra filename
		trafile = simfile.replace('.sim', '.tra')

		# Analyze the results of the .tra file
		# EventAnalysis.performCompleteAnalysis(filename=trafile)


	return
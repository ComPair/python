#!/usr/bin/env python
import fileinput
import glob
import sys

# Energy resolution
# Sigma vs Energy for cos Angle

# Angular resolution
# FWHM vs Energy for cos Angle

# Effective area
# (NumberOfReconstructedEvents/NumberOfSimulatedEvents) * math.pi * 300**2

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

	files = glob.glob(filename)

	for file in files:

		RMS = []
		FWHM = []
		Containment68 = None

		# Read the file line by line
		for line in fileinput.input([file]):

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


		fileinput.close()


		if 'Pair' in interactionType:

			return Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, Containment68, NumberOfReconstructedEvents

		else:

			return Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS[-2], FWHM[-2], NumberOfReconstructedEvents


##########################################################################################

def run(simulationFilename):

	# Create an dictionary to store the data
	results = {}

	# if 'all' in filename:
	# 	filename = 
	# 	parseMimrecLog('*.log')
	# else:
	# 	parseMimrecLog(filename)

	# Generate the mimrec log names
	mimrecFilename_tracked = simulationFilename.replace('.sim', '.mimrec_tracked.log')
	mimrecFilename_untracked = simulationFilename.replace('.sim', '.mimrec_untracked.log')
	mimrecFilename_pair = simulationFilename.replace('.sim', '.mimrec_pair.log')

	# Generate a key value
	key = simulationFilename.split('.inc')[0]

	# Enter an empty list for the initial value to which the results will be appended
	results[key] = []

	# Get the number of simulated events and add the value to the results dictionary
	NumberOfSimulatedEvents = parseSimLog(simulationFilename)
	results[key].append(NumberOfSimulatedEvents)

	# Get the mimrec figures of merit and add the dictionary
	print "Parsing: %s" % mimrecFilename_tracked
	Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents = parseMimrecLog(mimrecFilename_tracked, interactionType='Compton')
	results[key].append([Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents])

	# Get the mimrec figures of merit and add the dictionary
	print "Parsing: %s" % mimrecFilename_untracked	
	Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents = parseMimrecLog(mimrecFilename_untracked, interactionType='Compton')
	results[key].append([Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, RMS, FWHM, NumberOfReconstructedEvents])

	# Get the mimrec figures of merit and add the dictionary
	print "Parsing: %s" % mimrecFilename_pair		
	Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, Containment68, NumberOfReconstructedPairEvents = parseMimrecLog(mimrecFilename_pair, interactionType='Pair')
	results[key].append([Constant, ConstantError, Mean, MeanError, Sigma, SigmaError, Containment68, NumberOfReconstructedPairEvents])

	return results

##########################################################################################

if __name__=='__main__':


	if len(sys.argv) == 1:
		print 'Usage: FigureOfMeritParser.py filename'

	else:

		# Get the user supplied simulation output
		simulationFilename = sys.argv[1]

		results = run(simulationFilename)

		print ""
		print results


	







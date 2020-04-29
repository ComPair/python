#!/usr/bin/env python

'''
Tiny helper script for running performCompleteAnalysis on Revan output files. 
'''

import EventAnalysis      
import sys

if __name__ == "__main__":
	if len(sys.argv) != 2 and len(sys.argv) != 4:
		print(f"Usage: {sys.argv[0]} <path to directory | .tra file> [low bound for fit in keV] [upper bound for fit in keV]")
		exit(-1)

	fname = sys.argv[1]
	if len(sys.argv) == 4:
		elow = int(sys.argv[2])
		ehigh = int(sys.argv[3])	
		bounds = [elow, ehigh]
	else:
		bounds = None

	print(fname)
	print(fname[-4::])

	if fname[-4::] == '.tra':
		EventAnalysis.performCompleteAnalysis(filename=fname, 
			showPlots=True, energyRangeCompton=bounds)
	else:
		EventAnalysis.performCompleteAnalysis(directory=fname, 
			showPlots=False, energyRangeCompton=bounds)
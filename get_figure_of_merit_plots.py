
# Author: Grant Sommer (grant.sommer@gwu.edu)
# Date: January 15, 2026

# This code uses outputs of runCosima.sh and runRevan.sh to plot and save Figures of Merit

import os
import argparse
import EventAnalysis
import FigureOfMeritPlotter

def complete_analysis(directory):
    EventAnalysis.getTriggerEfficiency(directory=directory)
    EventAnalysis.performCompleteAnalysis(directory=directory, showPlots=False)

def figures_of_merit(directory):
    triggerEfficiencyFilename =f'{directory}/TriggerEfficiency.txt'
    data = FigureOfMeritPlotter.parseEventAnalysisLogs(directory, triggerEfficiencyFilename=triggerEfficiencyFilename)
    # angle selections aren't literally the angles selected, but a list like [1.0,0.9,0.8,0.7,....] that label the selected angles
    FigureOfMeritPlotter.plotAngularResolution(data, angleSelections=[1.0],save=True,working_directory=directory)
    FigureOfMeritPlotter.plotEnergyResolution(data, angleSelections=[1.0],save=True, working_directory=directory)
    FigureOfMeritPlotter.plotEffectiveArea(data, angleSelections=[1.0],save=True, show=False, working_directory=directory)

    # FigureOfMeritPlotter.plotAngularResolutionVsAngle(data, energySelections=[1.0],working_directory=directory)
    # FigureOfMeritPlotter.plotEnergyResolutionVsAngle(data, energySelections=[1.0],working_directory=directory)
    # FigureOfMeritPlotter.plotEffectiveAreaVsAngle(data, energySelections=[1.0],working_directory=directory)

def main(working_directory):
    # working_directory=os.getcwd()
    # complete_analysis(working_directory)
    figures_of_merit(working_directory)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("d", help="Working directory")
    args = parser.parse_args()
    working_directory=args.d

    main(working_directory)
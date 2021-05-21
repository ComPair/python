import EventAnalysis, FigureOfMeritPlotter, os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", default=os.getcwd(), help="Directory with .sim and .tra files. (Default: current directory)")
parser.add_argument("-t", "--tag", default="xxx", help="Tag for output txt files. (Default: xxx)")
args = parser.parse_args()


#this is needed so that the plots will be created in the directory with the .sim and .tra files, not the directory with the scripts.
os.chdir(args.dir )

if not os.path.exists( "./TriggerEfficiency.txt" ):

    EventAnalysis.getTriggerEfficiency( directory = "./" )

EventAnalysis.performCompleteAnalysis( directory="./", showPlots = False )

data = FigureOfMeritPlotter.parseEventAnalysisLogs("./", triggerEfficiencyFilename="./TriggerEfficiency.txt")

sel= [1.0, 0.8]

FigureOfMeritPlotter.plotAngularResolution(data, angleSelections=sel, save=True, txtOutfileLabel=args.tag)
FigureOfMeritPlotter.plotEnergyResolution(data, angleSelections=sel, save=True, txtOutfileLabel=args.tag)
FigureOfMeritPlotter.plotEffectiveArea(data, angleSelections=sel, ylog=True, save=True, txtOutfileLabel=args.tag)

for angle in sel:
    FigureOfMeritPlotter.plotSourceSensitivity(data, angleSelection=angle, save=True, uniterg=False, exposure = 0.2 * 3 * 365.25 * 24 * 3600, txtOutfileLabel=args.tag)


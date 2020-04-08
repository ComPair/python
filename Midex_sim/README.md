Performance plots.

Step 1: create_source_file.py\
This python code will generate the source files.

Step 2: create_cosima_run.py\
This python code will generate runCosima.sh.

Step 3: runCosima.sh\
e.g., cosima -s 120 FarFieldPointSource_0.158MeV_Cos1.0.source

Step 4: create_revan_run.py\
This python code will generate runRevan.sh.

Step 5: runRevan.sh\
e.g., revan -a -n -f FarFieldPointSource_0.158MeV_Cos1.0.inc1.id1.sim -g /data/slag2/dtak/Geometry/AMEGO_Midex/AmegoBase.geo.setup -c FarFieldPointSource_0.158MeV_Cos1.0.cfg

Step 6: EventAnalysis.py\
import EventAnalysis\
EventAnalysis.getTriggerEfficiency(directory='/data/slag2/dtak/PerformancePlotSimFiles/')\
events = EventAnalysis.performCompleteAnalysis('FarFieldPointSource_5.011MeV_Cos1.0.inc1.id1.tra', showPlots=True)\
Recommend running 'performCompleteAnalysis' for an individual tra file.

Step 7: FigureOfMeritPlotter.py\
triggerEfficiencyFilename ='/data/slag2/dtak/PerformancePlotSimFiles/TriggerEfficiency.txt'
directory='/data/slag2/dtak/PerformancePlotSimFiles/'\
data = FigureOfMeritPlotter.parseEventAnalysisLogs(directory, triggerEfficiencyFilename=triggerEfficiencyFilename)

Step 8: Check the log file and modify them if necessary.

Step 9: Generate plots\
FigureOfMeritPlotter.plotEffectiveArea(data,angleSelections=1.0, ylog=True)\
FigureOfMeritPlotter.plotEnergyResolution(data, angleSelections=1.0)\
FigureOfMeritPlotter.plotAngularResolution(data, angleSelections=1.0)

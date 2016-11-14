# python   
python scripts for manipulation of MEGAlib output  

#Files of note
EventViewer.py: Makes event display using output sim files. Currently only works for ComPair, not AMEGO

**EventViewer.py Usage Example**:  
import EventViewer  
filename = 'MyComPair_Tower.inc1.id1.sim'  
EventViewer.plot(filename)  

EventAnalysis.py: Analyzes the output from revan and creates performance and log files

**EventAnalysis.py Usage Example**:
EventAnalysis.performCompleteAnalysis(directory="../Simulations/AMEGO4x4PerformancePlotTraFiles/")

FigureOfMeritPlotter.py: Takes output log files from EventAnalysis and creates instrument performance files

please see the ipython notebook for example usage



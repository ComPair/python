"""### # ! /usr/bin/env python
"""
"""
------------------------------------------------------------------------

A script to replace the energy and angular resolution measurements performed by mimrec 

Author: Daniel Kocevski (dankocevski@gmail.com)
Date: September 3rd, 2016

Modified by Donggeun Tak (takdg123@gmail.com)
Date: April 3rd, 2020

Additional modifications by Sean Griffin (sean.c.griffin@gmail.com)
Data: April 10, 2020

Additional modifications by Donggeun Tak (takdg123@gmail.com)
Data: April 15, 2020

Usage Examples:
import EventAnalysis

# Parse the .tra file obtained from revan
events = EventAnalysis.parse('EffectiveArea_2MeV.inc1.id1.tra', sourceTheta=0.9)
 
# Calculate the angular resolution measurement (ARM) for Compton events
FWHM_ARM, dphi = EventAnalysis.getARMForComptonEvents(events, numberOfBins=100, phiRadius=5)
 
# Calculate the angular resolution measurement (ARM) for pair events
angles, openingAngles, contaimentData_68, contaimentBinned_68 = EventAnalysis.getARMForPairEvents(events, numberOfBins=100)

# Calculate the energy resolution for Compton events
mean, FWHM = EventAnalysis.getEnergyResolutionForComptonEvents(events, numberOfBins=100, energyPlotRange=[0,10000], energyFitRange=[1800,2100])
 
# Calculate the energy resolution for Pair events
fitMax, FWHM = EventAnalysis.getEnergyResolutionForPairEvents(events, numberOfBins=100, energyPlotRange=[0,10000], energyFitRange=[1800,2100])

# Display some diagnostic plots
EventAnalysis.plotDiagnostics(events)

# Visualize the pair events
EventAnalysis.visualizePairs(events, numberOfPlots=10)

------------------------------------------------------------------------
"""

import os
import time
import sys
import fileinput
import numpy
from itertools import product, combinations, repeat
from collections import OrderedDict
from scipy.stats import norm
from scipy.optimize import leastsq
import scipy.optimize
import math
import glob
import gzip

try:
    import matplotlib.pyplot as plot
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.mlab as mlab
    import matplotlib.gridspec as gridspec
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib.ticker import AutoMinorLocator
    from matplotlib.colors import LogNorm
    # Set the default title font dict
    titleFormat = {'fontsize': 12, 'fontweight' : plot.rcParams['axes.titleweight'], 'verticalalignment': 'baseline', 'horizontalalignment': 'center'}
except:
    print("\n**** Warning: matplotlib not found. Do not try to make plots or bad things will happen! ****")

try:
    from IPython.display import clear_output
except:
    pass
##########################################################################################

def getDetailsFromFilename(filename):

    '''Function to get the energy and angle from a filename.
    Really should be meta data.'''

    details = {}
    filepath, filename = os.path.split(filename)
    info = filename.split('_')

    details['MeV'] = info[1][:-3]

    angle = info[2].split('.')
    details['Cos'] = "{}.{}".format(angle[0][3:], angle[1])

    return details

def RadialGaussian(x, norm, sigma):
    """
    RadialGaussian(x, norm, sigma)
    Symmetric gaussian distribution in radial coordinates
    Funtcion Parameters:
    norm: Overall normalization
    sigma: sigma (width) of the underlying gaussian distribution
    """

    return x * norm / ((2*numpy.pi)*sigma**2) * numpy.exp(-numpy.power( x/sigma, 2) / 2.0 )

def RadialDoubleGaussian(x, norm1, sigma1, norm2, sigma2):
    """
    RadialGaussian(x, norm, sigma)
    Symmetric gaussian distribution in radial coordinates
    Funtcion Parameters:
    norm: Overall normalization
    sigma: sigma (width) of the underlying gaussian distribution
    """

    return RadialGaussian(x, norm1, sigma1 ) + RadialGaussian(x, norm2, sigma2)
    

def RadialKing(x, norm, sigma, gamma):
    
    return x * norm * 2 *numpy.pi * kingFunction(x, sigma, gamma)
    
def RadialDoubleKing(x, norm, f1, sigma1, gamma1, sigma2, gamma2):
    
    return x * norm * 2 *numpy.pi * (f1 * kingFunction(x, sigma1, gamma1) + (1-f1) * kingFunction(x, sigma2, gamma2)  )
    


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

    plot.plot(x,y)
    plot.show()

##########################################################################################

def skewedGaussian(x,h=1, e=0,w=1,a=0):
    t = (x-e) / w
    # return 2 / w * scipy.stats.norm.pdf(t) * scipy.stats.norm.cdf(a*t)
    return h * 2 * scipy.stats.norm.pdf(t) * scipy.stats.norm.cdf(a*t)

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
        fig = plot.figure()
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
        fig = plot.figure()
        ax = fig.add_subplot(111, projection='3d')

    ax.plot_wireframe(x*radius, y*radius, z*radius, color="gray")


##########################################################################################

def angularSeparation(vector1, vector2):

    # Calculate the product of the vector magnitudes
    product = numpy.linalg.norm(vector2) * numpy.linalg.norm(vector1)

    # Make sure we don't devide by zero
    if product != 0:

        # Calculate the dot product
        dotProduct = numpy.dot(vector2, vector1)

        # Make sure we have sane results
        value = dotProduct/product 
        if (value >  1.0): value =  1.0;
        if (value < -1.0): value = -1.0;

        # Get the reconstructed angle (in degrees)
        angle = numpy.degrees(numpy.arccos(value))

    else:

        # Return zero in case the denominator is zero
        angle = 0.0


    return angle


##########################################################################################

def extractCoordinates(coordinateData):

    x = []
    y = []
    z = []

    for coordinates in coordinateData:
        x.append(coordinates[0])
        y.append(coordinates[1])
        z.append(coordinates[2])

    return x, y, z

##########################################################################################

def plotPairConversionCoordinates(events):

    # Get the individual pair conversion coordinates
    x, y, z = extractCoordinates(events['position_pairConversion'])

    plot.figure(figsize=[10,9], color='#3e4d8b', alpha=0.9, histtype='stepfilled')
    plot.hist(x, bins=50)
    plot.xlabel('x')

    plot.figure(figsize=[10,9], color='#3e4d8b', alpha=0.9, histtype='stepfilled')
    plot.hist(y, bins=50)
    plot.xlabel('y')

    plot.figure(figsize=[10,9], color='#3e4d8b', alpha=0.9, histtype='stepfilled')
    plot.hist(z, bins=50)
    plot.xlabel('z')


##########################################################################################

def parse(filename, sourceTheta=1.0, testnum=-1):

    # Create the dictionary that will contain all of the results
    events = {}

    # Define the lists to store energy information for Compton events
    energy_ComptonEvents = []
    energy_ComptonEvents_error = [] 
    energy_TrackedComptonEvents = []
    energy_UntrackedComptonEvents = []
    energy_firstScatteredPhoton = []
    energy_firstScatteredPhoton_error = []
    energy_recoiledElectron = []
    energy_recoiledElectron_error = []

    # Define the lists to store energy information for pair events
    energy_pairElectron = []
    energy_pairElectron_error = []
    energy_pairPositron = []
    energy_pairPositron_error = []
    energy_pairDepositedInFirstLayer = []
    energy_pairDepositedInFirstLayer_error = []

    # Define the lists to store position information
    position_originalPhoton = []
    position_firstInteraction = []
    position_firstInteraction_error = []
    position_secondInteraction = []
    position_secondInteraction_error = []
    position_pairConversion = []
    position_firstInteraction = []
    position_firstInteraction_error = []

    # Define the lists to store the direction vectors
    direction_recoilElectron = []
    direction_pairElectron = []
    direction_pairPositron = []

    # Define other lists
    phi_Tracker = []
    qualityOfComptonReconstruction = []
    qualityOfPairReconstruction = []


    # Read the number of lines in the file
    command = 'wc %s' % filename
    output = os.popen(command).read()
    totalNumberOfLines = float(output.split()[0])

    # Start the various event and line counters
    numberOfUnknownEventTypes = 0
    numberOfComptonEvents = 0
    numberOfPairEvents = 0
    numberOfPhotoElectricEffectEvents = 0
    numberOfTrackedElectronEvents = 0
    numberOfUntrackedElectronEvents = 0
    numberOfBadEvents = 0
    lineNumber = 0
    eventNumber = 0

    # Create empty index lists to track event types
    index_tracked = []
    index_untracked = []

    # track time of event and times between events
    time = []
    dt = []
    tn=0

    # Start by collecting all events
    skipEvent = False

    # Loop through the .tra file
    print('\nParsing: %s' % filename)

    if filename[-2:] == "gz": #compressed file
        lines = gzip.open(filename, mode="rt")
    else:
        lines = fileinput.input([filename])
    
    for line in lines:

        try:
            sys.stdout.write("Progress: %d%%   \r" % (lineNumber/totalNumberOfLines * 100) )
            sys.stdout.flush()
        except:
            pass

        if 'TI' in line:
            lineContents = line.split() 
            time.append(float(lineContents[1]))

        if testnum > 0:

            if tn == testnum:
                break
            tn=tn+1

        if 'ET ' in line:

            # Split the line
            lineContents = line.split() 

            eventType = lineContents[1]

            # Skip the event if its of unknown type
            if 'UN' in eventType:

                numberOfUnknownEventTypes = numberOfUnknownEventTypes + 1
                skipEvent = True
                continue

            else:
                skipEvent = False

            if 'CO' in eventType:
                # Increment the Compton event counter
                numberOfComptonEvents = numberOfComptonEvents + 1

            if 'PA' in eventType:
                # Increment the pair event counter
                numberOfPairEvents = numberOfPairEvents + 1

            if 'PH' in eventType:
                # Increment the photo electron effect counter
                numberOfPhotoElectricEffectEvents = numberOfPhotoElectricEffectEvents + 1

        if 'BD' in line:

            # Split the line
            lineContents = line.split()

            whybad = lineContents[1]

            if whybad != 'None':
                # Events don't get reconstructed for a number of reasons
                # In pair, it's mostly 'TooManyHistInCSR'
                numberOfBadEvents = numberOfBadEvents + 1

        ####### Compton Events #######

        # Extract the Compton event energy information
        if 'CE ' in line and skipEvent == False:

            # Split the line
            lineContents = line.split() 

            # Get the energy of the first scattered gamma-ray
            energy_firstScatteredPhoton.append(float(lineContents[1]))
            energy_firstScatteredPhoton_error.append(float(lineContents[2]))

            # Get the energy of the recoiled electron
            energy_recoiledElectron.append(float(lineContents[3]))
            energy_recoiledElectron_error.append(float(lineContents[4]))


        # Extract the Compton event hit information
        if 'CD ' in line and skipEvent == False:

            # Split the line
            lineContents = line.split() 

            # Get the position of the first scattered gamma-ray
            x1 = float(lineContents[1])
            y1 = float(lineContents[2])
            z1 = float(lineContents[3])

            # Get the position uncertainty of the first scattered gamma-ray     
            x1_error = float(lineContents[4])
            y1_error = float(lineContents[5])
            z1_error = float(lineContents[6])

            # Get the position of the second scattered gamma-ray
            x2 = float(lineContents[7])
            y2 = float(lineContents[8])
            z2 = float(lineContents[9])

            # Get the position uncertainty of the second scattered gamma-ray        
            x2_error = float(lineContents[10])
            y2_error = float(lineContents[11])
            z2_error = float(lineContents[12])

            # Get the origin position of the original gamma-ray
            x0 = x1
            y0 = y1
            z0 = 1000.0

            # Get the position of the second scattered gamma-ray
            x_electron = float(lineContents[13])
            y_electron = float(lineContents[14])
            z_electron = float(lineContents[15])

            # Get the position uncertainty of the second scattered gamma-ray        
            x_electron_error = float(lineContents[16])
            y_electron_error = float(lineContents[17])
            z_electron_error = float(lineContents[18])

            # Record the energy of the Compton event (regardless of whether it was has a tracked or untrack electron)
            energy_ComptonEvents.append(energy_firstScatteredPhoton[-1] + energy_recoiledElectron[-1])

            # Record the energy error of the Compton event (regardless of whether it was has a tracked or untrack electron)
            energy_ComptonEvents_error.append( math.sqrt( energy_firstScatteredPhoton_error[-1]**2 + energy_recoiledElectron_error[-1]**2 ) )


            # Determine if the recoil electron was tracked by the detector
            if x_electron != 0:

                # Record the direction of the recoil electron
                direction_recoilElectron.append([x_electron, y_electron, z_electron])

                # Record the energy of tracked events
                energy_TrackedComptonEvents.append(energy_firstScatteredPhoton[-1] + energy_recoiledElectron[-1])

                # Record the number of tracked events
                numberOfTrackedElectronEvents = numberOfTrackedElectronEvents + 1

                # Add this event to the index of tracked events
                index_tracked.append(eventNumber)

            else:

                # Record the direction of the recoil electron
                direction_recoilElectron.append([numpy.nan,numpy.nan,numpy.nan])

                # Record the energy of tracked events
                energy_UntrackedComptonEvents.append(energy_firstScatteredPhoton[-1] + energy_recoiledElectron[-1])

                # Record the number of untracked events
                numberOfUntrackedElectronEvents = numberOfUntrackedElectronEvents + 1

                # Add this event to the index of untracked events
                index_untracked.append(eventNumber)


            # Store the origin coordinates of the original gamma-ray 
            position0 = numpy.array([x0, y0, z0])
            # position0Error = numpy.array([x1_error,y1_error,z1_error])

            # Get the x-axis offset based on the theta of the source.  This assumes phi=0
            # Note: For Compton events adjusting the theta of the source happens in the parser
            # for pair events, it happens in the ARM calculation. 
            dx = numpy.tan(numpy.arccos(sourceTheta)) * (z1 - z0)

            # Set the origin position of the original gamma-ray
            position0 = numpy.array([x1-dx, y1, z0])

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
                dotProduct = numpy.dot(directionVector2, directionVector1)

                # Make sure we have sane results
                value = dotProduct/product 
                if (value >  1.0): value =  1.0;
                if (value < -1.0): value = -1.0;

                # Get the reconstructed phi angle (in radians) and add the known angle to the source 
                phi_Tracker.append(numpy.arccos(value)) 

            else:

                # Return zero in case the denominator is zero
                phi_Tracker.append(0.0)

            # Add the position information to their respective list of positions
            position_originalPhoton.append(position0)
            position_firstInteraction.append(position1)
            position_firstInteraction_error.append(position1Error)
            position_secondInteraction.append(position2)
            position_secondInteraction_error.append(position2Error)

            # Increment the event number
            eventNumber = eventNumber + 1


        ####### Pair Events #######

        # Extract the pair conversion hit information
        if 'PC ' in line and skipEvent == False:

            # Split the line
            lineContents = line.split() 

            # Get the position of the pair conversion
            x1 = float(lineContents[1])
            y1 = float(lineContents[2])
            z1 = float(lineContents[3])

            # Save the position
            position_pairConversion.append([x1,y1,z1])



        # Extract the pair electron information
        if 'PE ' in line and skipEvent == False:


            # Split the line
            lineContents = line.split()

            # Get the electron information
            energy_pairElectron.append(float(lineContents[1]))
            energy_pairElectron_error.append(float(lineContents[2]))

            # Get the direction of the pair electron
            x = float(lineContents[3])
            y = float(lineContents[4])
            z = float(lineContents[5])

            # Store the direction of the pair electron
            direction_pairElectron.append([x,y,z])

        # Extract the pair positron information
        if 'PP ' in line and skipEvent == False:

            # Split the line
            lineContents = line.split()

            # Get the electron information
            energy_pairPositron.append(float(lineContents[1]))
            energy_pairPositron_error.append(float(lineContents[2]))

            # Get the direction of the pair electron
            x = float(lineContents[3])
            y = float(lineContents[4])
            z = float(lineContents[5])

            # Store the direction of the pair electron
            direction_pairPositron.append([x,y,z])


        if 'PI ' in line and skipEvent == False:

            # Split the line
            lineContents = line.split()

            # Get the electron information
            energy_pairDepositedInFirstLayer.append(float(lineContents[1]))




        # Extract the reconstruction quality
        if 'TQ ' in line and skipEvent == False:

            # Split the line
            lineContents = line.split()

            if 'CO' in eventType:

                # Get the reconstruction quality
                qualityOfComptonReconstruction.append(float(lineContents[1]))

            if 'PA' in eventType:

                # Get the reconstruction quality
                qualityOfPairReconstruction.append(float(lineContents[1]))

        # Increment the line number for the progress indicator
        lineNumber = lineNumber + 1


    # Add all of the lists to the results dictionary
    events['energy_ComptonEvents'] = numpy.array(energy_ComptonEvents).astype(float)
    events['energy_ComptonEvents_error'] = numpy.array(energy_ComptonEvents_error).astype(float)
    events['energy_TrackedComptonEvents'] = numpy.array(energy_TrackedComptonEvents).astype(float)
    events['energy_UntrackedComptonEvents'] = numpy.array(energy_UntrackedComptonEvents).astype(float)
    events['energy_firstScatteredPhoton'] = numpy.array(energy_firstScatteredPhoton).astype(float)
    events['energy_firstScatteredPhoton_error'] = numpy.array(energy_firstScatteredPhoton_error).astype(float)
    events['energy_recoiledElectron'] = numpy.array(energy_recoiledElectron).astype(float)
    events['energy_recoiledElectron_error'] = numpy.array(energy_recoiledElectron_error).astype(float)
    events['energy_pairElectron'] = numpy.array(energy_pairElectron).astype(float)
    events['energy_pairElectron_error'] = numpy.array(energy_pairElectron_error).astype(float)
    events['energy_pairPositron'] = numpy.array(energy_pairPositron).astype(float)
    events['energy_pairPositron_error'] = numpy.array(energy_pairPositron_error).astype(float)
    events['energy_pairDepositedInFirstLayer'] = numpy.array(energy_pairDepositedInFirstLayer).astype(float)
    events['energy_pairDepositedInFirstLayer_error'] = numpy.array(energy_pairDepositedInFirstLayer_error).astype(float)

    events['position_originalPhoton'] = numpy.array(position_originalPhoton).astype(float)
    events['position_firstInteraction'] = numpy.array(position_firstInteraction).astype(float)
    events['position_firstInteraction_error'] = numpy.array(position_firstInteraction_error).astype(float)
    events['position_secondInteraction'] = numpy.array(position_secondInteraction).astype(float)
    events['position_secondInteraction_error'] = numpy.array(position_secondInteraction_error).astype(float)
    events['position_pairConversion'] = numpy.array(position_pairConversion).astype(float)

    events['direction_recoilElectron'] = numpy.array(direction_recoilElectron).astype(float)
    events['direction_pairElectron'] = numpy.array(direction_pairElectron).astype(float)
    events['direction_pairPositron'] = numpy.array(direction_pairPositron).astype(float)

    events['phi_Tracker'] = numpy.array(phi_Tracker).astype(float)
    events['numberOfComptonEvents'] =numberOfComptonEvents
    events['numberOfUntrackedElectronEvents'] =numberOfUntrackedElectronEvents
    events['numberOfTrackedElectronEvents'] =numberOfTrackedElectronEvents
    events['numberOfPairEvents'] = numberOfPairEvents
    events['numberOfUnknownEventTypes'] = numberOfUnknownEventTypes
    events['numberOfBadEvents'] = numberOfBadEvents
    events['index_tracked'] = numpy.array(index_tracked)
    events['index_untracked']= numpy.array(index_untracked)
    events['qualityOfComptonReconstruction'] = numpy.array(qualityOfComptonReconstruction).astype(float)
    events['qualityOfPairReconstruction'] = numpy.array(qualityOfPairReconstruction).astype(float)
    events['time'] = numpy.array(time).astype(float)
    events['deltime'] = numpy.append(0.,events['time'][1:]-events['time'][0:len(events['time'])-1])

    # Print some event statistics
    print("\n\nStatistics of Event Selection")
    print("***********************************")
    print("Total number of analyzed events: %s" % (numberOfComptonEvents + numberOfPairEvents))

    if numberOfComptonEvents + numberOfPairEvents == 0:
        print("No events pass selection")
        events=False
        return events

    print("")
    print("Number of unknown events: %s (%i%%)" % (numberOfUnknownEventTypes, 100*numberOfUnknownEventTypes/(numberOfComptonEvents + numberOfPairEvents + numberOfUnknownEventTypes)))
    print("Number of pair events: %s (%i%%)" % (numberOfPairEvents, 100*numberOfPairEvents/(numberOfComptonEvents + numberOfPairEvents + numberOfUnknownEventTypes)))
    print("Number of Compton events: %s (%i%%)" % (numberOfComptonEvents, 100*numberOfComptonEvents/(numberOfComptonEvents + numberOfPairEvents + numberOfUnknownEventTypes)))
    if numberOfComptonEvents > 0:
        print(" - Number of tracked electron events: %s (%i%%)" % (numberOfTrackedElectronEvents, 100.0*(float(numberOfTrackedElectronEvents)/numberOfComptonEvents)))
        print(" - Number of untracked electron events: %s (%i%%)" % (numberOfUntrackedElectronEvents, 100*(float(numberOfUntrackedElectronEvents)/numberOfComptonEvents)))
    print("")
    print("")

    fileinput.close()

    return events


##########################################################################################

def getARMForComptonEvents(events, numberOfBins=100, phiRadius=10, onlyTrackedElectrons=False, onlyUntrackedElectrons=False, showPlots=True, filename=None, energyCutSelection=False, energyCut = [0, 0]):

    # Set some constants
    electron_mc2 = 511.0        # KeV

    # Retrieve the event data
#    energyCutSelection=False
    if energyCutSelection:
        energySelection = (events['energy_ComptonEvents']<(energyCut[0]+3*energyCut[1])) * (events['energy_ComptonEvents']>(energyCut[0]-3*energyCut[1]))
        index_tracked = [i for i in events['index_tracked'] if energySelection[i]==True]
        index_untracked = [i for i in events['index_untracked'] if energySelection[i]==True]
        
        energy_firstScatteredPhoton = events['energy_firstScatteredPhoton'][energySelection]
        energy_recoiledElectron = events['energy_recoiledElectron'][energySelection]
        phi_Tracker = events['phi_Tracker'][energySelection]
            
        if onlyTrackedElectrons == True:
            energy_firstScatteredPhoton = events['energy_firstScatteredPhoton'][index_tracked]
            energy_recoiledElectron = events['energy_recoiledElectron'][index_tracked]
            phi_Tracker = events['phi_Tracker'][index_tracked]
            if onlyUntrackedElectrons == True:
                print("select either tracked or untracked compton events!!")
                exit()
            
        if onlyUntrackedElectrons == True:
            energy_firstScatteredPhoton = events['energy_firstScatteredPhoton'][index_untracked]
            energy_recoiledElectron = events['energy_recoiledElectron'][index_untracked]
            phi_Tracker = events['phi_Tracker'][index_untracked]
            if onlyTrackedElectrons == True:
                print("select either tracked or untracked compton events!!")
                exit()
            
    else:
        energy_firstScatteredPhoton = events['energy_firstScatteredPhoton']
        energy_recoiledElectron = events['energy_recoiledElectron']
        phi_Tracker = events['phi_Tracker']
        index_tracked = events['index_tracked']
        index_untracked=events['index_untracked']

        # Determine whether to include only Tracked or Untracked electron events
        if onlyTrackedElectrons == True:
            energy_firstScatteredPhoton = energy_firstScatteredPhoton[index_tracked]
            energy_recoiledElectron = energy_recoiledElectron[index_tracked]
            phi_Tracker = phi_Tracker[index_tracked]
            if onlyUntrackedElectrons == True:
                print("select either tracked or untracked compton events!!")
                exit()

        if onlyUntrackedElectrons == True:
            energy_firstScatteredPhoton = energy_firstScatteredPhoton[index_untracked]
            energy_recoiledElectron = energy_recoiledElectron[index_untracked]
            phi_Tracker = phi_Tracker[index_untracked]
            if onlyTrackedElectrons == True:
                print("select either tracked or untracked compton events!!")
                exit()

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
    plot.figure(figsize=[10,9])
    # plot.rcParams['figure.figsize'] = 10, 9

    gs = gridspec.GridSpec(4,1)
    ax1 = plot.subplot(gs[:3, :])
    
    if len(dphi[selection]) < 10:
        print("The number of events is not sufficient (<10)")
        return numpy.nan, numpy.nan
    elif len(dphi[selection]) < 500:
        print("The number of events is not sufficient, so that 'numberOfBins' and 'phiRadius' changed.")
        numberOfBins = 25
        
    # Create the histogram
    histogram_angleResults = ax1.hist(dphi[selection], numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
    ax1.set_xlim([-1*phiRadius,phiRadius])

    # Extract the binned data and bin locations
    dphi_binned = histogram_angleResults[0]
    if max(dphi_binned) < 10:
        numberOfBins = 20
        histogram_angleResults = ax1.hist(dphi[selection], numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
        ax1.set_xlim([-1*phiRadius,phiRadius])
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
    try:
        optimizedParameters, covariance = scipy.optimize.curve_fit(DoubleLorentzian, bincenters, dphi_binned, [offset, mean, width1, height1, width2, height2])

        # Calculate the optimized curve
        y_fit = DoubleLorentzian(bincenters, optimizedParameters[0], optimizedParameters[1], optimizedParameters[2], optimizedParameters[3], optimizedParameters[4], optimizedParameters[5])

        # Get the fwhm of the fit
        x1 = bincenters[numpy.where(y_fit >= numpy.max(y_fit)/2)[0][0]]
        x2 = bincenters[numpy.where(y_fit >= numpy.max(y_fit)/2)[0][-1]]    

        mean = optimizedParameters[1]
        FWHM = x2-x1

        # Annotate the plot
        ax1.text(0.03, 0.9, "Mean = %.3f deg\nFWHM = %.3f deg" % (mean, FWHM), verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes, color='black', fontsize=12)

        # Plot the data
        ax1.plot(bincenters, y_fit, color='darkred', linewidth=2)
        ax1.plot([x1,x2],[numpy.max(y_fit)/2.,numpy.max(y_fit)/2.], color='darkred', linestyle='--', linewidth=2)
        titlestr = ''
        if onlyTrackedElectrons:
            titlestr = "Tracked"
        if onlyUntrackedElectrons:
            titlestr = "Untracked"
        if not onlyTrackedElectrons and not onlyUntrackedElectrons:
            titlestr = "Tracked + Untracked"
        ax1.set_title(f'Angular Resolution (Compton Events)\n{titlestr}', fontdict=titleFormat)

        # ax1.axes.xaxis.set_ticklabels([])

        # Create a subplot for the residuals 
        ax2 = plot.subplot(gs[3, :])

        # Plot the residuals
        ax2.step(bincenters, dphi_binned-y_fit, color='#3e4d8b', alpha=0.9) 
        ax2.plot([bincenters[0],bincenters[-1]], [0,0], color='darkred', linewidth=2, linestyle='--')
        ax2.set_xlim([-1*phiRadius,phiRadius])
        ax2.set_xlabel('ARM - Compton Cone (deg)')

        # Set the minor ticks
        ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
        ax2.xaxis.set_minor_locator(AutoMinorLocator(4))

    except:

        print("**** Warning: fit failed to converge! ****")
        mean = numpy.nan
        FWHM = numpy.nan


    # Print some statistics
    print("\n\nStatistics of ARM histogram and fit (Compton Events)")
    print("***********************************")
    print("Compton and pair events in ARM histogram: %s (%s%%)" % ( len(dphi[selection]), 100*len(dphi[selection])/(len(dphi)) ))
    print("")
    print("Mean of fit: %s" % mean) 
    print("FWHM of fit: %s" %  FWHM)
    print("")

    if filename is not None:
        f=getDetailsFromFilename(filename)
        f1,f2 = f['MeV'], f['Cos']

        ft=None
        if onlyTrackedElectrons == True:
            ft='tracked'
        else:
            ft='untracked'

        plot.savefig("%sMeV_Cos%s_angular_resolution_%s.png" % (f1,f2,ft))

    # Show the plot
    if showPlots == True:
        plot.show()
    else:
        plot.close()


    return FWHM, dphi


########################################################

def kingFunction(x, sigma, gamma):

    K = (1. / (2. * numpy.pi * sigma**2.)) * (1. - (1./gamma)) * ((1 + (1/(2*gamma)) * (x**2/sigma**2))**(-1*gamma))

    return K

########################################################

def fcore(N, sigma_tail, sigma_core):

    value = 1 / (1 + (N * sigma_tail**2 / sigma_core**2))

    return value

##########################################################################################

def latScaleFactor(energy, c0=0.153, c1=0.0057000001, beta=-0.8):

    scaleFactor = numpy.sqrt( (c0 * ((energy/100)**beta) )**2 + c1**2)

    return scaleFactor

##########################################################################################

def latPSF(scaledDeviation, NTAIL, STAIL, SCORE, GTAIL, GCORE):

    psf = (fcore(NTAIL, STAIL, SCORE) * kingFunction(scaledDeviation, SCORE, GCORE)) + (( 1 - fcore(NTAIL, STAIL, SCORE) ) * kingFunction(scaledDeviation, STAIL, GTAIL))

    return psf

##########################################################################################


def getScaledDeviation(events, sourceTheta=0):

    # Get the scaled angular distribution
    scaledDeviation, openingAngles, contaimentData_68, contaimentBinned_68 = getARMForPairEvents(events, sourceTheta=sourceTheta, numberOfBins=100, angleFitRange=[0,10**1.5], anglePlotRange=[0,10**1.5], openingAngleMax=180., showPlots=True, numberOfPlots=0, finishExtraction=True, qualityCut=1, energyCut=numpy.nan, weightByEnergy=True, showDiagnosticPlots=False, filename=None, log=True, getScaledDeviation=True)

    # angles, openingAngles, contaimentData_68, contaimentBinned_68 = EventAnalysis.getARMForPairEvents(events, sourceTheta=sourceTheta, numberOfBins=100, angleFitRange=[0,10**1.5], anglePlotRange=[0,10**1.5], openingAngleMax=180., showPlots=True, numberOfPlots=0, finishExtraction=True, qualityCut=1, energyCut=numpy.nan, weightByEnergy=True, showDiagnosticPlots=False, filename=None, log=True, getScaledDeviation=True)

    return scaledDeviation


##########################################################################################

def getARMForPairEvents(events, sourceTheta=0, numberOfBins=100, angleFitRange=[0,30], anglePlotRange=[0,30], openingAngleMax=180., showPlots=True, numberOfPlots=0, finishExtraction=True, qualityCut=1, energyCut=numpy.nan, weightByEnergy=True, showDiagnosticPlots=True, filename=None, log=False, getScaledDeviation=False, onlyangles=False):


    # Define the list to contain the resulting angle measurements
    angles = []
    openingAngles = []
    distances = []

    # Start some counters
    plotNumber = 0
    numberOfRejectedEvents = 0

    # Define a default position to the source
    dz = 1000.0
    dx = 0
    dy = 0

    # Create a list to contain an index of events passing the quality cut
    index_goodQuality = []

    # Loop through each event and calculate the reconstructed photon direction and it's offset to the true direction
    for index in range(events['numberOfPairEvents']):

        # Make a cut on the quality of the pair reconstruction
        if events['qualityOfPairReconstruction'][index] > qualityCut:
            numberOfRejectedEvents = numberOfRejectedEvents + 1
            continue

        # Get the electron and positron energies
        energy_pairElectron = events['energy_pairElectron'][index]
        energy_pairPositron = events['energy_pairPositron'][index]

        # Get total energy
        energy_PairSum = energy_pairElectron + energy_pairPositron

        # Make a cut on the total summed energy of the pair
        if energy_PairSum < energyCut:
            numberOfRejectedEvents = numberOfRejectedEvents + 1
            continue

        # Register this pair as 'good'
        index_goodQuality.append(index)

        # Get the position of the gamma conversion
        position_conversion = events['position_pairConversion'][index]

        # Get the x-axis offset based on the theta of the source.  This assumes phi=0
        # Note: For Compton events adjusting the theta of the source happens in the parser
        # for pair events, it happens here in the ARM calculation. 
        dx = numpy.tan(numpy.radians(sourceTheta)) * (position_conversion[2] - dz)

        # Set the origin position of the original gamma-ray
        position_source = [position_conversion[0]-dx, position_conversion[1], dz]

        # Get the electron and positron direction vectors. These are unit vectors.
        direction_electron = events['direction_pairElectron'][index]
        direction_positron = events['direction_pairPositron'][index]

        # Weight the direction vectors by the electron and positron energies
        if weightByEnergy == True:      
            direction_electron = direction_electron * (energy_pairElectron/energy_PairSum)
            direction_positron = direction_positron * (energy_pairPositron/energy_PairSum)

        # Get the vector that bisects the electron and positron vectors
        direction_bisect = (direction_electron + direction_positron)/2.0

        # Invert the bisect vector to obtain the reconstructed source vector
        direction_source_reconstructed = -1*direction_bisect

        # Calculate the vector between the first interaction and the origin of the original gamma-ray
        direction_source = -1*(position_conversion - position_source)

        # Calculate the distance of the conversion point to the top-center of the spacecraft
        position_topCenter = numpy.array([0,0,60])
        distance = numpy.linalg.norm(position_conversion-position_topCenter)
        distances.append(distance)

        # Calculate the product of the vector magnitudes
        #print direction_source, direction_source_reconstructed
        angle = angularSeparation(direction_source, direction_source_reconstructed)
                
        #angle_shift = angle-float(sourceTheta)
        #print angle, angle_shift

        # Normalized by the incoming photon's reconstructed energy (for PSF fitting)
        if getScaledDeviation == True:

            # Convert from keV to MeV 
            energy_PairSum_MeV = energy_PairSum * 0.001 

            # Calculate the scale factor
            scaleFactor = latScaleFactor(energy_PairSum_MeV)

            # Get the energy scaled angular seperation
            angle = angle / scaleFactor


        # Store the angle
        angles.append(angle)

        openingAngle = angularSeparation(direction_electron, direction_positron)
        openingAngles.append(openingAngle)

        # # Calculate the product of the vector magnitudes
        # product = numpy.linalg.norm(direction_source) * numpy.linalg.norm(direction_source_reconstructed)

        # # Make sure we don't devide by zero
        # if product != 0:

        #   # Calculate the dot product
        #   dotProduct = numpy.dot(direction_source, direction_source_reconstructed)

        #   # Make sure we have sane results
        #   value = dotProduct/product 
        #   if (value >  1.0): value =  1.0;
        #   if (value < -1.0): value = -1.0;

        #   # Get the reconstructed angle (in degrees)
        #   angle = numpy.degrees(numpy.arccos(value))

        #   # Store the angle 
        #   angles.append(angle)

        # else:

        #   # Return zero in case the denominator is zero
        #   angle = 0.0
        #   angles.append(angle)


        if plotNumber != numberOfPlots:

            print('Plot number: %s' % index)
            print("")
            print("Photon conversion coordiantes: %s, %s, %s" % (position_conversion[0], position_conversion[1], position_conversion[2]))
            print("")
            print('Source vector (True): %s, %s, %s' % (direction_source[0], direction_source[1], direction_source[2]))
            print('Source vector (Reconstructed): %s, %s, %s' % (direction_source_reconstructed[0], direction_source_reconstructed[1], direction_source_reconstructed[2]))
            print('Angle = %s' % angle)

            fig = plot.figure()
            ax = fig.add_subplot(111, projection='3d')


            # Plot the geometry
            plotCube(shape=[50*2,50*2,35.75*2], position=[0,0,35.0], color='red', ax=ax)
            plotCube(shape=[40*2,40*2,30*2], position=[0,0,29.25], color='blue', ax=ax)
            plotCube(shape=[40*2,40*2,5.0*2], position=[0,0,-8.0], color='green', ax=ax)

            # Set the plot limits
            ax.set_xlim3d(-60,60)
            ax.set_ylim3d(-60,60)
            ax.set_zlim3d(-50,100)

            # Set the plot labels
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')


            ax.scatter( position_conversion[0], position_conversion[1], position_conversion[2], marker='o', label='Pair')
            # Plot the reconstructed photon direction and the true photon direction
            ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], direction_source[0], direction_source[1], direction_source[2], pivot='tail', arrow_length_ratio=0, color='purple', linestyle='--', length=100, label=r'$\gamma$')
            # ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], direction_source_reconstructed[0], direction_source_reconstructed[1], direction_source_reconstructed[2], pivot='tail', arrow_length_ratio=0, color='green', linestyle='--', length=100, label=r"$\gamma$'")

            # Plot the electron and positron vectors                            
            ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], direction_electron[0], direction_electron[1], direction_electron[2], pivot='tail', arrow_length_ratio=0.05, color='darkblue', length=50, label=r'$e^{-}$')
            ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], direction_positron[0], direction_positron[1], direction_positron[2], pivot='tail', arrow_length_ratio=0.05, color='darkred', length=50, label=r'$e^{+}$')

            # Plot the vector bisecting the electron and positron trajectories
            # ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], direction_bisect[0], direction_bisect[1], direction_bisect[2], pivot='tail', arrow_length_ratio=0, color='green', linestyle='--', length=50)
            # ax.quiver( position_conversion[0], position_conversion[1], position_conversion[2], -1*direction_bisect[0], -1*direction_bisect[1], -1*direction_bisect[2], pivot='tail', arrow_length_ratio=0, color='green', linestyle='--', length=50)

            # Plot the legend
            handles, labels = ax.get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys(), scatterpoints=1)

            plotNumber = plotNumber + 1
            plot.show()

        elif finishExtraction == False:

            return

    if showDiagnosticPlots == True:

        # Get the sum energy
        energy_PairSum = events['energy_pairElectron'] + events['energy_pairPositron']

        # Generate the color map
        norm = plot.Normalize()
        colors = plot.cm.jet(norm(energy_PairSum))
        # colors = plot.cm.jet(norm(numpy.log(energy_PairSum)))

        # plot.hist(energy_PairSum, bins=100, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
        # plot.xlabel('Energy (keV)')
        # plot.show()

        # plot.hist(openingAngles, bins=150, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
        # plot.xlabel('Opening Angle (deg)')
        # plot.xlim([0,45])
        # plot.show()


        # plot.figure(figsize=[10,7])
        # plot.scatter(angles, distances, marker='.', s=0.5, color='#3e4d8b')
        # plot.xlabel('Angular Resolution (deg)')       
        # plot.ylabel('distance')       
        # plot.show()

        x, y, z = extractCoordinates(events['position_pairConversion'])

        plot.figure(figsize=[10,7])
        plot.scatter(angles, z, marker='.', s=0.5, color='#3e4d8b')
        plot.xlabel('Angular Resolution (deg)')     
        plot.ylabel('z')        
        plot.show()


        plot.figure(figsize=[10,7])
        cm = plot.cm.get_cmap('jet')
        plot.scatter(angles, events['qualityOfPairReconstruction'], c=energy_PairSum, s=10, cmap=cm, edgecolor='none')
        cbar = plot.colorbar(pad = 0.02)
        cbar.set_label('Energy (keV)')
        plot.xlabel('Angular Resolution (deg)')
        plot.ylabel('Quality Factor')
        plot.xlim([0,10])
        plot.ylim([0,1])

        plot.figure(figsize=[10,7])
        plot.hist2d(angles, events['qualityOfPairReconstruction'], bins=300, norm=LogNorm())
        plot.xlabel('Angular Resolution (deg)')     
        plot.ylabel('Quality Factor')       
        plot.xlim([0,10])
        plot.ylim([0,1])
        plot.show()

        plot.figure(figsize=[11,9])
        plot.scatter(angles, events['energy_pairElectron'] / events['energy_pairPositron'], marker='.', s=0.5, color='#3e4d8b')
        plot.xlabel('Angular Resolution (deg)')
        plot.ylabel(r'$E_{e-}$/$E_{e+}$')       
        plot.yscale('log')
        plot.xlim([-5,90])
        plot.ylim([1e-3,1e3])
        plot.show()

        plot.figure(figsize=[11,9])
        plot.scatter(angles, events['energy_pairElectron'] + events['energy_pairPositron'], marker='.', s=0.5, color='#3e4d8b')
        plot.xlabel('Angular Resolution (deg)')     
        plot.ylabel(r'$E_{e-}$ + $E_{e+}$')     
        # plot.yscale('log')
        plot.xlim([-5,90])  
        plot.ylim([0,1e5])      
        plot.show()

        plot.figure(figsize=[11,9])
        plot.scatter(angles, events['energy_pairDepositedInFirstLayer'], marker='.', s=0.5, color='#3e4d8b')
        plot.xlabel('Angular Resolution (deg)')     
        plot.ylabel(r'Deposited Energy (keV)')      
        # plot.yscale('log')
        plot.xlim([-5,90])  
        # plot.ylim([0,1e5])        
        plot.show()


        plot.hist(events['energy_pairElectron'][index_goodQuality]/(events['energy_pairElectron'][index_goodQuality]+events['energy_pairPositron'][index_goodQuality]), bins=100, alpha=0.8, color='#3e4d8b', label=r'$e^{-}$')
        plot.hist(events['energy_pairPositron'][index_goodQuality]/(events['energy_pairElectron'][index_goodQuality]+events['energy_pairPositron'][index_goodQuality]), bins=100, alpha=0.8, color='darkred', label=r'$e^{+}$')
        plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper center')  
        plot.xlabel('Fractional Energy')
        ax = plot.subplot(111)
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax.yaxis.set_minor_locator(AutoMinorLocator(4))
        plot.show()

        plot.hist(events['energy_pairElectron'][index_goodQuality], bins=100, alpha=0.8, color='#3e4d8b', label=r'$e^{-}$')
        plot.hist(events['energy_pairPositron'][index_goodQuality], bins=100, alpha=0.75, color='darkred', label=r'$e^{+}$')
        plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper right')   
        plot.xlabel('Energy (KeV)')
        ax = plot.subplot(111)
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax.yaxis.set_minor_locator(AutoMinorLocator(4))     
        plot.show()

        plot.figure(figsize=[10,7])
        ax2=plot.subplot(111)
        ARMvOpeningAngle=ax2.plot(angles,openingAngles,'ro')
        ax2.set_xlabel('Angular Resolution (deg)')
        ax2.set_ylabel('Opening Angle (deg)')
        plot.show()
        #plot.clf()



    # Convert the list of angles and opening angles to numpy arrays 
    angles = numpy.array(angles)
    openingAngles = numpy.array(openingAngles)

    if onlyangles:
        return numpy.array(angles)


    # # Select the events within the desired quality range
    # selection_quality = numpy.where( events['qualityOfPairReconstruction'] <= qualityCut )

    # # Apply the selection filter
    # angles_unfiltered = angles
    # angles = angles[selection_quality]

    # Select the events within the desired angle range and opening angle
    selection_fit = numpy.where( (angles >= angleFitRange[0]) & (angles <= angleFitRange[1]) & (openingAngles <= openingAngleMax))
    oa_len = len(angles[selection_fit])

    # Apply the fit selection
    angles_fit = angles[selection_fit]

    if len(angles_fit) < 10:
        print("The number of events is too small (<10)")
        return numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, []

    # Set the plot size
    plot.figure(figsize=[10,7])

    gs = gridspec.GridSpec(4,1)
    ax1 = plot.subplot(gs[:3, :])

    # Create the axis
    # ax1 = plot.subplot(111)

    containmentData_90 = numpy.quantile( angles_fit, .90 )
    anglePlotRange[1] = containmentData_90

    # Create the histogram
    histogramResults = ax1.hist(angles_fit, bins=numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
    ax1.set_xlabel('Angular resolution (deg)')
    ax1.set_xlim(anglePlotRange)
    ax1.set_title('Angular Resolution (Pair Events)', fontdict=titleFormat)

    # Setup the minor ticks
    ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(4))

    # Convert the list of angles into a numpy array
    angles_fit = numpy.array(angles_fit)

    # Sort the angle array
    angles_fit.sort()

    # Find the 68% containment of the cumulative sum of the angle distribution
    contaimentData_68 = angles_fit[int(len(angles_fit)*.68)]

    # Extract the binned data and bin locations
    angles_binned = histogramResults[0]
    bins = histogramResults[1]
    bincenters = 0.5*(bins[1:]+bins[:-1])

    # Find the peak of the historgram
    angles_max = numpy.max(angles_binned)

    # Find the cumulative sum of the angle distribution and its max
    angles_binned_cumulativeSum = numpy.cumsum(angles_binned)
    angles_binned_cumulativeMax = angles_binned_cumulativeSum[-1]

    # Find the index that corresponds to 68% of the max
    index_68 = numpy.where(angles_binned_cumulativeSum >= angles_binned_cumulativeMax*0.68)[0][0]

    # Get the 68% containment of the cumulative sum of the binned angle distribution
    contaimentBinned_68 = bincenters[index_68]

    # Add the containment values to the plot
    ax1.axvline(contaimentData_68, color='green', linewidth=1.5, linestyle='--', label="68%% (data): %.2f deg" % contaimentData_68)
    ax1.axvline(contaimentBinned_68, color='#3e4d8b', linewidth=1.5, linestyle='--', label="68%% (histogram): %.2f deg" %  contaimentBinned_68)

    # Change to log scaling
    if log == True:
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        
    
    #Radial 2D gaussian is the product of 2 normal distributions, so 1 sigma corresponds to the 68% * 68% = 47% containment radius ec.
    #68% containment radius is given by from scipy.stats import norm; norm.ppf( (numpy.sqrt(0.68) + 1 )/2.0 ) = 1.3551214210715499 sigma
    sigma2D_to_68 = 1.3551214210715499
    
    try:
    #if True:
    
        #startingValues = [numpy.sum(angles_binned), contaimentBinned_68 / sigma2D_to_68, numpy.sum(angles_binned), 5.0  ]
        
        #startingValues = [ numpy.sum(angles_binned), contaimentBinned_68 / sigma2D_to_68, 3 ]
        
        startingValues = [ numpy.sum(angles_binned), 0.5, 0.5*contaimentBinned_68, 3, 2*contaimentBinned_68, 2 ]
       
       
        #F = RadialDoubleGaussian
        #F = RadialKing
        F = RadialDoubleKing
        optimizedParameters, covariance = scipy.optimize.curve_fit(F, bincenters, angles_binned, startingValues)

        # Calculate the optimized curve
        y_fit = F(bincenters, *optimizedParameters)

        angles_fit_cumulativeSum = numpy.cumsum(y_fit)
        angles_fit_cumulativeMax = angles_fit_cumulativeSum[-1]

        # Find the index that corresponds to 68% of the max
        index_68 = numpy.where(angles_fit_cumulativeSum >= angles_fit_cumulativeMax*0.68)[0][0]

        # Get the 68% containment of the cumulative sum of the binned angle distribution
        containmentFit_68 = bincenters[index_68]

        # Annotate the plot
        ax1.axvline(containmentFit_68, color='darkred', linewidth=1.5, linestyle='--', label="68%% (%s fit): %.2f deg" % (F.__name__, containmentFit_68 ))

        # Plot the fit results
        ax1.plot(bincenters, y_fit, color='darkred', linewidth=2)

        # Create a subplot for the residuals
        
        
        ax2 = plot.subplot(gs[3, :])

        # Plot the residuals
        ax2.step(bincenters, angles_binned-y_fit, color='#3e4d8b', alpha=0.9)
        ax2.plot([bincenters[0],bincenters[-1]], [0,0], color='darkred', linewidth=2, linestyle='--')
        ax2.set_xlim(anglePlotRange)
        ax2.set_ylabel( "Data - Fit" )
        ax2.set_xlabel('Angular resolution (deg)')

        # Set the minor ticks
        #ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
        #ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
        #ax2.xaxis.set_minor_locator(AutoMinorLocator(4))

    except Exception as E:
    #if False:

        print("**** Warning: fit failed to converge! ****")
        print(E)
        optimizedParameters = []
        containmentFit_68 = numpy.nan

    ax1.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper right')


    # Print some statistics
    print("\n\nStatistics of ARM histogram and fit (Pair Events)")
    print("*****************************************************")
    print("")
    print("Total number of pair events: %s" % events['numberOfPairEvents'])
    print("Number of pair events passing quality cut  : %s (%s%%)" % ( len(angles), 100*len(angles)/(events['numberOfPairEvents']) ) )
    print("Number of pair events passing opening angle: %s (%s%%)" % ( oa_len, 100*oa_len/(events['numberOfPairEvents']) ) )
    print("Number of pair events included in ARM histogram: %s (%s%%)" % ( len(angles_fit), 100*len(angles_fit)/(len(angles_fit)) ) )
    print("")
    print("Maximum: %s" % angles_max)
    print("68%% containment (data): %.2f deg" % contaimentData_68)
    print("68%% containment (histogram): %.2f deg" %  contaimentBinned_68)
    print("")
    print("%s fit parameters:" % F.__name__, optimizedParameters)
    print("68%% containment (fit): %.2f deg" %  containmentFit_68)

    if filename is not None:
        f=getDetailsFromFilename(filename)
        f1,f2 = f['MeV'], f['Cos']

        plot.savefig("%sMeV_Cos%s_angular_resolution_Pair.png" % (f1,f2))

    # Show the plot
    if showPlots == True:
        plot.show()
    else:
        plot.close()

    return angles, openingAngles, contaimentData_68, contaimentBinned_68, containmentFit_68, optimizedParameters


##########################################################################################

def getEnergyResolutionForPairEvents(events, numberOfBins=100, energyPlotRange=None, energyFitRange=False, showPlots=True, qualityCut=1.0,fileBase="5.011MeV_Cos0.8"):

    # Retrieve the event data
    energy_pairElectron = events['energy_pairElectron']
    energy_pairPositron = events['energy_pairPositron']
    energy_pairElectron_error = events['energy_pairElectron_error']
    energy_pairPositron_error = events['energy_pairPositron_error']
    qualityOfPairReconstruction = events['qualityOfPairReconstruction']

    # Estimate the energy of the incoming photon and its associated error
    energy_pairReconstructedPhoton = energy_pairElectron + energy_pairPositron
    energy_pairReconstructedPhoton_error = numpy.sqrt((energy_pairElectron_error/energy_pairElectron)**2 + (energy_pairPositron_error/energy_pairPositron)**2) * energy_pairReconstructedPhoton

    # Find maximum bin value to fit
    n, b, patches = plot.hist(energy_pairReconstructedPhoton,bins=numberOfBins)
    bin_max=numpy.where(n == n.max())
    if energyFitRange:
        minbin=0.5*b[bin_max][0]
        energyFitRange=[0.5*b[bin_max][0],1.5*b[bin_max][0]]
    else:
        energyFitRange=[0,1e5]

    # Select the events within the desired energy range
    selection = numpy.where( (energy_pairReconstructedPhoton >= energyFitRange[0]) & (energy_pairReconstructedPhoton <= energyFitRange[1]) & (qualityOfPairReconstruction <= qualityCut))

    # Create histrograms for the electrons and positrons
    histElectrons = plot.hist(energy_pairElectron[selection], 
                  bins=numberOfBins, 
                  alpha=0.5, histtype='stepfilled',label="Electrons")
    histPositrons = plot.hist(energy_pairPositron[selection], 
                  bins=numberOfBins, 
                  alpha=0.5, histtype='stepfilled',label="Positrons")
    plot.xlabel('Energy (keV)')
    plot.xlim(energyPlotRange)
    plot.legend()
    plot.savefig(fileBase+"_pairEnergies.png")
    plot.close()

    if len(energy_pairReconstructedPhoton[selection]) < 10:
        print("The number of events is too small (<10)")
        return numpy.nan, numpy.nan, numpy.nan
    elif len(energy_pairReconstructedPhoton[selection]) < 100:
        numberOfBins = 30
    
    # Create the histogram
    histogramResults = plot.hist(energy_pairReconstructedPhoton[selection], bins=numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
    plot.xlabel('Energy (keV)')
    plot.xlim(energyPlotRange)
    plot.savefig(fileBase+"_photonEnergies.png")
    plot.close()

    # Extract the binned data and bin locations
    energy_binned = histogramResults[0]
    bins = histogramResults[1]
    bincenters = 0.5*(bins[1:]+bins[:-1])

    # Get the bin center containing the maximum of the histogram
    bin_max = bincenters[numpy.argmax(energy_binned)]

    # Set the initial parameters
    height = numpy.max(energy_binned)   # h
    scale = 1e4                         # w
    location = bin_max                  # epsilon
    shape = -2                          # alpha

    # Fit the histogram data
    try:
        optimizedParameters, covariance = scipy.optimize.curve_fit(skewedGaussian, bincenters, energy_binned, [height, scale, location, shape])
    except Exception as msg:
        print("{}".format(msg))
    # Calculate the optimized curve
    try:
        y_fit = skewedGaussian(bincenters, optimizedParameters[0], optimizedParameters[1], optimizedParameters[2], optimizedParameters[3])
    except Exception as msg:
        print("{}".format(msg))
        return numpy.nan, numpy.nan, numpy.nan

    # Get the max of the fit
    fitMax = bincenters[numpy.argmax(y_fit)]

    # Get the fwhm of the fit
    x1 = bincenters[numpy.where(y_fit >= numpy.max(y_fit)/2)[0][0]]
    x2 = bincenters[numpy.where(y_fit >= numpy.max(y_fit)/2)[0][-1]]
    FWHM = x2-x1

    # Approximate sigma as FWHM/2
    sigma=FWHM/2.

    # Print some statistics
    print("\n\nStatistics of energy histogram and fit (pair events)")
    print("********************************************************")
    print("Number of Compton and pair events in histogram: %s (%s%%)" % ( len(energy_pairReconstructedPhoton[selection]), 100*len(energy_pairReconstructedPhoton[selection])/(len(energy_pairReconstructedPhoton)) ))
    print("")
    print("Fitting in range: ", energyFitRange[0], energyFitRange[1])
    print("Max of fit: %s keV" % fitMax)    
    print("FWHM of fit: %s keV" % FWHM)
    print("sigma of fit: %s keV" % sigma)
    print("")


    # Set the plot size
    plot.figure(figsize=[10,9])

    # Setup the plot grid
    gs = gridspec.GridSpec(4,1)

    # Create the first axis
    ax1 = plot.subplot(gs[:3, :])

    # Plot the histogram
    histogramResults = plot.hist(energy_pairReconstructedPhoton[selection], bins=numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
    ax1.set_title('Energy Resolution (Pair Events)', fontdict=titleFormat)
    if energyPlotRange != None:
        ax1.set_xlim(energyPlotRange)

    # Overplot the fit
    ax1.plot(bincenters, y_fit, color='darkred', linewidth=2)
    ax1.plot([x1,x2],[numpy.max(y_fit)/2.,numpy.max(y_fit)/2.], color='darkred', linestyle='--', linewidth=2)

    # Annotate the plot
    ax1.text(0.03, 0.9, "Max = %.3f keV\nFWHM = %.3f keV" % (fitMax, FWHM), verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes, color='black', fontsize=12)

    # Create a subplot for the residuals 
    ax2 = plot.subplot(gs[3, :])

    # Plot the residuals
    ax2.step(bincenters, energy_binned-y_fit, color='#3e4d8b', alpha=0.9)   
    ax2.plot([bincenters[0],bincenters[-1]], [0,0], color='darkred', linewidth=2, linestyle='--')
    ax2.set_xlabel('Energy (keV)')
    if energyPlotRange != None:     
        ax2.set_xlim(energyPlotRange)

    ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))

    # Show the plot
    if showPlots == True:
        plot.show()
    else:
        plot.close()

    return fitMax, FWHM, sigma

##########################################################################################

def getEnergyResolutionForComptonEvents(events, numberOfBins=100, energyHardCut = 5, energyPlotRange=None, energyFitRange=None, onlyTrackedElectrons=False, onlyUntrackedElectrons=False, showPlots=False, filename=None, inputEnergy=None):

    # Retrieve the event data
    energy_ComptonEvents = events['energy_ComptonEvents']
    numberOfComptonEvents = len(events['index_tracked'])+len(events['index_untracked'])
    numberOfPairEvents = events['numberOfPairEvents']

    # Determine whether to include only Tracked or Untracked electron events
    if onlyTrackedElectrons == True:
        energy_ComptonEvents  = events['energy_TrackedComptonEvents']
        numberOfComptonEvents = len(events['index_tracked']) 
        if onlyUntrackedElectrons == True:
            print("select either tracked or untracked compton events!!")
            exit()

    if onlyUntrackedElectrons == True:
        energy_ComptonEvents  = events['energy_UntrackedComptonEvents']
        numberOfComptonEvents = len(events['index_untracked'])
        if onlyTrackedElectrons == True:
            print("select either tracked or untracked compton events!!")
            exit()


    # Convert the list to a numpy array
    energy_ComptonEvents = numpy.array(energy_ComptonEvents)

    if numberOfComptonEvents < 20:
        print("The number of events is too small (<20)")
        return numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan
    elif numberOfComptonEvents < 100:
        numberOfBins = 30


    # Select the events within the desired energy range
    if energyPlotRange != None:
        energySelection_plot = numpy.where((energy_ComptonEvents >= energyPlotRange[0]) & (energy_ComptonEvents <= energyPlotRange[1]))
    else:
        energySelection_plot = numpy.arange(len(energy_ComptonEvents))

    # Create the binned data
    histogram_energyResults = plot.hist(energy_ComptonEvents[energySelection_plot], numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
    plot.close()
    
    '''
    # Extract the binned data and bin locations
    energy_binned = histogram_energyResults[0]
    
    if numpy.argmax(energy_binned) < 10:
        numberOfBins = 8
        #print(numpy.argmax(energy_binned))
        histogram_energyResults = plot.hist(energy_ComptonEvents[energySelection_plot], numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
        plot.close()
        #print(numpy.argmax(energy_binned))
    '''
    
    energy_binned = histogram_energyResults[0]
    bins_energy = histogram_energyResults[1]
    bincenters_energy = 0.5*(bins_energy[1:]+bins_energy[:-1])
    bin_start = 0
    # Set the range of the energy fit by finding the inflection points in the histogram. This will not work with poorly sampled data
    if energyFitRange == None:
        if max(energy_binned)>10:
            if inputEnergy == None:
                bin_max = numpy.argmax(energy_binned)
            else:
                bin_max = numpy.argmin(numpy.abs(bins_energy - inputEnergy*1000))-1
                
            print(numpy.argmax(energy_binned), numpy.argmin(numpy.abs(bins_energy - inputEnergy*1000))-1)
            for i in range(bin_max):
                if energy_binned[bin_max] > 500:
                    if (energy_binned[bin_max-i-1] > energy_binned[bin_max-i]) or (energy_binned[bin_max-i-1]< (energy_binned[bin_max]*0.7)):
                        bin_start = bincenters_energy[bin_max-i]
                        break
                else:
                    if ((energy_binned[bin_max-i-1] > energy_binned[bin_max-i]) and (energy_binned[bin_max-i-1]< (energy_binned[bin_max]/2.))) or (energy_binned[bin_max-i-1]< (energy_binned[bin_max]/2.)):
                        bin_start = bincenters_energy[bin_max-i]
                        if (bin_start > bincenters_energy[numpy.argmax(energy_binned)]*(1-energyHardCut/100.)):
                            bin_start = bincenters_energy[numpy.argmax(energy_binned)]*(1-energyHardCut/100.)
                        break

            for i in range(len(bins_energy-1)-bin_max-2):
                if (energy_binned[bin_max+i+1] > energy_binned[bin_max+i]):
                    bin_stop = bincenters_energy[bin_max+i]
                    break
                else:
                    bin_stop = bincenters_energy[bin_max+i]
                    
            
            
            energyFitRange = [min(bin_start,bincenters_energy[numpy.argmax(energy_binned)]*0.97) , max(bin_stop, bincenters_energy[numpy.argmax(energy_binned)]*1.15)]
            
        else:
            energyFitRange =[bincenters_energy[numpy.argmax(energy_binned)]*0.91, bincenters_energy[numpy.argmax(energy_binned)]*1.15]
    
    # Set the range of the energy fit to a hard coded percentage of the maximum value
    # if energyFitRange == None:
        # energyFitRange = [bincenters_energy[numpy.argmax(energy_binned)]*0.91, bincenters_energy[numpy.argmax(energy_binned)]*1.15]

    # Fit a gaussian to the energy data within the user specified energy range
    energySelection_fit1 = numpy.where( (energy_ComptonEvents >= energyFitRange[0]) & (energy_ComptonEvents <= energyFitRange[1]) ) 
    mu_Guassian, sigma = norm.fit(energy_ComptonEvents[energySelection_fit1])

    # Create the Guassian fit line
    x = numpy.array(range(int(energyFitRange[0]),int(energyFitRange[1]), 1))
    y_fit = norm.pdf( x, mu_Guassian, sigma)

    # Calculate the Gaussian fit statistics
    FWHM_Guassian = 2*math.sqrt(2*math.log(2))*sigma


    ############  Begin asymmetric fit to the binned data ############  

    # Create the binned data
    
    #if len(energy_ComptonEvents)>10000 and False:
    #    energy_ComptonEvents_temp = energy_ComptonEvents[(energy_ComptonEvents > energyFitRange[0])*(energy_ComptonEvents <energyFitRange[1])]
    #else:
    energy_ComptonEvents_temp = energy_ComptonEvents
        
    histogram_energyResults = plot.hist(energy_ComptonEvents_temp, numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')
    plot.close()    

    # Extract the binned data and bin locations
    energy_binned = histogram_energyResults[0]
    bins_energy = histogram_energyResults[1]
    bincenters = 0.5*(bins_energy[1:]+bins_energy[:-1])
    
    #if len(energy_ComptonEvents)>10000 and False:
    #    energy_binned_skew = histogram_energyResults[0]
    #    bins_energy_skew = histogram_energyResults[1]
    #    bincenters_skew = 0.5*(bins_energy[1:]+bins_energy[:-1])

    # Select only the data within the desires energy fit range 
    energySelection_fit = numpy.where( (bincenters >= energyFitRange[0]) & (bincenters <= energyFitRange[1]) )  
    energy_binned = energy_binned[energySelection_fit]
    bincenters = bincenters[energySelection_fit]

    # Get the bin center containing the maximum of the histogram
    bin_max = bincenters[numpy.argmax(energy_binned)]

    # Set the initial parameters
    # SG: Basing values on numbers from the 
    height = numpy.max(energy_binned)   # h
    scale = sigma #1e4                          # w
    location = mu_Guassian #bin_max                 # epsilon
    shape = -2                          # alpha

    # Fit the histogram data and Calculate the optimized curve to an asymmetric gaussian
    try:
        optimizedParameters, covariance = scipy.optimize.curve_fit(skewedGaussian, bincenters, energy_binned, [height, scale, location, shape])
        y_fit2 = skewedGaussian(bincenters, optimizedParameters[0], optimizedParameters[1], optimizedParameters[2], optimizedParameters[3])
    except Exception as msg:
        print("{}".format(msg))
        return numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan

    # Get the max of the fit
    fitMax_skewedGuassian = bincenters[numpy.argmax(y_fit2)]

    # Get the fwhm of the fit
    x1 = bincenters[numpy.where(y_fit2 >= numpy.max(y_fit2)/2)[0][0]]
    x2 = bincenters[numpy.where(y_fit2 >= numpy.max(y_fit2)/2)[0][-1]]
    FWHM_skewedGuassian = x2-x1

    ############  End asymmetric fit to the binned data ############  


    # Plot the histogram and the fit line, normalized to match the histogram data
    histogram_energyResults = plot.hist(energy_ComptonEvents[energySelection_plot], numberOfBins, color='#3e4d8b', alpha=0.9, histtype='stepfilled')

    #energy_binned = histogram_energyResults[0]
    #bins_energy = histogram_energyResults[1]
    #bincenters = 0.5*(bins_energy[1:]+bins_energy[:-1])
  
    titlestr = ''
    if onlyTrackedElectrons:
        titlestr = "Tracked"
    if onlyUntrackedElectrons:
        titlestr = "Untracked"
    if not onlyTrackedElectrons and not onlyUntrackedElectrons:
        titlestr = "Tracked + Untracked"
    plot.title(f'Energy Resolution (Compton Events)\n{titlestr}', fontdict=titleFormat)

    plot.xlabel('Energy (keV)')

    # Set the minor ticks
    ax = plot.gca()
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))

    # Overplot the Gaussian fit
    y_fitNormalized = y_fit/numpy.max(y_fit)*numpy.max(histogram_energyResults[0])
    plot.plot(x, y_fitNormalized, color='darkred', linewidth=2)
    plot.plot([mu_Guassian-(FWHM_Guassian/2.), mu_Guassian+(FWHM_Guassian/2.)], [numpy.max(y_fitNormalized)/2.,numpy.max(y_fitNormalized)/2.], color='darkred', linestyle='dotted', linewidth=2)

    # Overplot the asymetric Gaussian fit
    plot.plot(bincenters, y_fit2, color='green', linewidth=2)
 
    plot.plot([x1,x2],[numpy.max(y_fit2)/2.,numpy.max(y_fit2)/2.], color='green', linestyle='dotted', linewidth=2)

    # Plot the range of the fit
    plot.plot([energyFitRange[0],energyFitRange[0]], [0,ax.get_ylim()[1]], color='lightblue', linestyle='--', linewidth=1)
    plot.plot([energyFitRange[1],energyFitRange[1]], [0,ax.get_ylim()[1]], color='lightblue', linestyle='--', linewidth=1)

    # Annotate the plot
    ax0 = plot.subplot(111) 
    ax0.text(0.03, 0.97, "Gaussian\nMean = %.3f keV\nFWHM = %.3f keV" % (mu_Guassian, FWHM_Guassian), verticalalignment='top', horizontalalignment='left', transform=ax0.transAxes, color='black', fontsize=12)
    ax0.text(0.97, 0.97, "Skewed Gaussian\nMax = %.3f keV\nFWHM = %.3f keV" % (fitMax_skewedGuassian, FWHM_skewedGuassian), verticalalignment='top', horizontalalignment='right', transform=ax0.transAxes, color='black', fontsize=12)

    # Print some statistics
    print("\n\nStatistics of energy histogram and fit (Compton events)")
    print("***********************************************************")
    print("Compton and pair events in energy histogram: %s (%s%%)" % ( len(energy_ComptonEvents[energySelection_fit1]), 100*len(energy_ComptonEvents[energySelection_fit1])/(numberOfComptonEvents + numberOfPairEvents)))
    print("")
    print("Mean of Guassian fit: %s keV" % mu_Guassian)
    print("Sigma of Gaussian fit: %s keV" % sigma)
    print("FWHM of Guassian fit: %s keV" % FWHM_Guassian)
    print("")
    print("Max of asymmetric Guassian fit: %s keV" % fitMax_skewedGuassian)
    print("FWHM of asymmetric Guassian fit: %s keV" % FWHM_skewedGuassian)
    
    if filename is not None:
        f=getDetailsFromFilename(filename)
        f1,f2 = f['MeV'], f['Cos']

        if onlyTrackedElectrons:
            titlestr = "Tracked"
        if onlyUntrackedElectrons:
            titlestr = "Untracked"
        if not onlyTrackedElectrons and not onlyUntrackedElectrons:
            titlestr = "Both"
       
        plot.savefig(f"{f1}MeV_Cos{f2}_energy_resolution_{titlestr}.png")

    if showPlots == True:
        plot.show()
    else:
        plot.close()

    return mu_Guassian, FWHM_Guassian, fitMax_skewedGuassian, FWHM_skewedGuassian, sigma


##########################################################################################

def plotDiagnostics(events, showPlots=True):

    # Get the angular resolution measurements
    FWHM, dphi = getARMForComptonEvents(events, showPlots=False)

    # Retrieve the event information
    energy_ComptonEvents = events['energy_ComptonEvents']
    energy_ComptonEvents_error = events['energy_ComptonEvents_error']
    energy_firstScatteredPhoton = events['energy_firstScatteredPhoton']
    energy_firstScatteredPhoton_error = events['energy_firstScatteredPhoton_error'] 
    energy_recoiledElectron = events['energy_recoiledElectron']
    energy_recoiledElectron_error = events['energy_recoiledElectron_error']

    index_untracked = events['index_untracked']
    index_tracked = events['index_tracked']


    # Plot the conversion coordinates
    fig = plot.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the geometry
    plotCube(shape=[50*2,50*2,35.75*2], position=[0,0,35.0], color='red', ax=ax)
    plotCube(shape=[40*2,40*2,30*2], position=[0,0,29.25], color='blue', ax=ax)
    plotCube(shape=[40*2,40*2,5.0*2], position=[0,0,-8.0], color='green', ax=ax)

    # Set the plot limits
    ax.set_xlim3d(-60,60)
    ax.set_ylim3d(-60,60)
    ax.set_zlim3d(-50,100)

    # Set the plot labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    x = []
    y = []
    z = []
    for coordinates in events['position_pairConversion']:
        x.append(coordinates[0])
        y.append(coordinates[1])
        z.append(coordinates[2])

    ax.scatter(x, y, z, marker='.', s=1)
    plot.show()


    # Photon Energy vs ARM (Compton events)
    plot.figure(figsize=[11,9])
    if len(index_tracked) != 0:
        plot.scatter(dphi[index_tracked], energy_ComptonEvents[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
    if len(index_untracked) != 0:
        plot.scatter(dphi[index_untracked], energy_ComptonEvents[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
    plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
    plot.ylabel('Reconstructed Photon Energy (keV)')
    plot.xlabel(r'ARM ($\Delta\theta$)')

    # Set the minor ticks
    ax = plot.gca()
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))


    # Incoming photon energy error vs ARM (Compton events)
    plot.figure(figsize=[11,9])
    percentError_ComptonEvents = 100. * energy_ComptonEvents_error/energy_ComptonEvents
    if len(index_tracked) != 0:
        plot.scatter(dphi[index_tracked], percentError_ComptonEvents[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
    if len(index_untracked) != 0:
        plot.scatter(dphi[index_untracked], percentError_ComptonEvents[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
    plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
    plot.ylabel('Reconstructed Photon Energy Error (%)')
    plot.xlabel(r'ARM ($\Delta\theta$)')

    # Set the minor ticks
    ax = plot.gca()
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))


    # First scattered photon energy vs ARM (Compton events)
    plot.figure(figsize=[11,9])
    if len(index_tracked) != 0:
        plot.scatter(dphi[index_tracked], energy_firstScatteredPhoton[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
    if len(index_untracked) != 0:
        plot.scatter(dphi[index_untracked], energy_firstScatteredPhoton[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
    plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
    plot.ylabel('Scattered Photon Energy (keV)')
    plot.xlabel(r'ARM ($\Delta\theta$)')

    # Set the minor ticks
    ax = plot.gca()
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))


    # Electron recoil energy vs ARM (Compton events)
    plot.figure(figsize=[11,9])
    if len(index_tracked) != 0:
        plot.scatter(dphi[index_tracked], energy_recoiledElectron[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
    if len(index_untracked) != 0:
        plot.scatter(dphi[index_untracked], energy_recoiledElectron[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
    plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
    plot.ylabel('Electron Recoil Energy (keV)')
    plot.xlabel(r'ARM ($\Delta\theta$)')

    # Set the minor ticks
    ax = plot.gca()
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))


    # First scattered photon energy error vs ARM (Compton events)
    plot.figure(figsize=[11,9])
    percentError_firstScatteredPhoton = 100. * energy_firstScatteredPhoton_error/energy_firstScatteredPhoton
    if len(index_tracked) != 0:
        plot.scatter(dphi[index_tracked], percentError_firstScatteredPhoton[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
    if len(index_untracked) != 0:
        plot.scatter(dphi[index_untracked], percentError_firstScatteredPhoton[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
    plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
    plot.ylabel('Scattered Photon Energy Error (%)')
    plot.xlabel(r'ARM ($\Delta\theta$)')

    # Set the minor ticks
    ax = plot.gca()
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))


    # Electron recoil energy error vs ARM (Compton events)
    plot.figure(figsize=[11,9])
    percentError_recoiledElectron = 100. * energy_recoiledElectron_error/energy_recoiledElectron
    if len(index_tracked) != 0:
        plot.scatter(dphi[index_tracked], percentError_recoiledElectron[index_tracked], color='darkred', marker='.', s=0.5, label='tracked')
    if len(index_untracked) != 0:
        plot.scatter(dphi[index_untracked], percentError_recoiledElectron[index_untracked], color='#3e4d8b', marker='.', s=0.5, label='untracked')
    plot.legend(numpoints=1, scatterpoints=1, fontsize='medium', frameon=True, loc='upper left')
    plot.ylabel('Electron Recoil Energy Error (%)')
    plot.xlabel(r'ARM ($\Delta\theta$)')

    # Set the minor ticks
    ax = plot.gca()
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))


    # Show the plot
    if showPlots == True:
        plot.show()
    else:
        plot.close()


##########################################################################################

def visualizePairs(events, sourceTheta=0, numberOfPlots=10):

    getARMForPairEvents(events, sourceTheta=sourceTheta, numberOfPlots=numberOfPlots, finishExtraction=False, showDiagnosticPlots=False)

    return


##########################################################################################

def visualizeCompton(events, showEvent=1, onlyShowTracked=True):

    # Open a new 3D plot
    fig = plot.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the geometry
    plotCube(shape=[50*2,50*2,35.75*2], position=[0,0,35.0], color='red', ax=ax)
    plotCube(shape=[40*2,40*2,30*2], position=[0,0,29.25], color='blue', ax=ax)
    plotCube(shape=[40*2,40*2,5.0*2], position=[0,0,-8.0], color='green', ax=ax)

    # Set the plot limits
    ax.set_xlim3d(-60,60)
    ax.set_ylim3d(-60,60)
    ax.set_zlim3d(-50,100)

    # Set the plot labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')


    if onlyShowTracked == True:
        index_type = events['index_tracked']
    else:
        index_type = numpy.arange(len(events['index_tracked'])) 

    index = showEvent       

    position_originalPhoton = events['position_originalPhoton'][index_type][index]  
    position_originalPhoton[2] = 100

    position_firstInteraction = events['position_firstInteraction'][index_type][index]
    position_secondInteraction = events['position_secondInteraction'][index_type][index]
    direction_recoilElectron = events['direction_recoilElectron'][index_type][index]

    direction_firstInteraction =  position_firstInteraction - position_originalPhoton
    direction_secondInteraction = position_secondInteraction - position_firstInteraction


    print(position_originalPhoton)
    print(position_firstInteraction)
    print(position_secondInteraction)

    print(direction_firstInteraction)
    print(direction_secondInteraction)

    # Plot the interactions
    ax.scatter(position_originalPhoton[0], position_originalPhoton[1], position_originalPhoton[2], color='black', marker=u'$\u2193$')
    ax.scatter(position_firstInteraction[0], position_firstInteraction[1], position_firstInteraction[2], marker='s', color='darkblue', label='Comp')
    # ax.scatter(position_secondInteraction[0], position_secondInteraction[1], position_secondInteraction[2], marker='o', label='BREM')

    ax.plot( [position_originalPhoton[0], position_firstInteraction[0]], [position_originalPhoton[1],position_firstInteraction[1]], zs=[position_originalPhoton[2],position_firstInteraction[2]], linestyle='--', color='purple', label=r"$\gamma$")
    # ax.plot( [position_firstInteraction[0], position_secondInteraction[0]], [position_firstInteraction[1],position_secondInteraction[1]], zs=[position_firstInteraction[2],position_secondInteraction[2]] )

    # Plot the electron and positron vectors                            
    # ax.quiver( position_originalPhoton[0], position_originalPhoton[1], position_originalPhoton[2], direction_firstInteraction[0], direction_firstInteraction[1], direction_firstInteraction[2], pivot='tail', arrow_length_ratio=0.05, color='darkblue', length=50)
    ax.quiver( position_firstInteraction[0], position_firstInteraction[1], position_firstInteraction[2], direction_secondInteraction[0], direction_secondInteraction[1], direction_secondInteraction[2], pivot='tail', arrow_length_ratio=0.10, color='purple',length=25, label=r"$\gamma$'")
    ax.quiver( position_firstInteraction[0], position_firstInteraction[1], position_firstInteraction[2], direction_recoilElectron[0], direction_recoilElectron[1], direction_recoilElectron[2], pivot='tail', arrow_length_ratio=0.10, color='darkblue', length=25, label=r'e$^{-}$')

    # Plot the legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), scatterpoints=1)

    plot.show()

    return

 
##########################################################################################


def performCompleteAnalysis(filename=None, directory=None, energies=None, angles=None, showPlots=False, energySearchUnit='MeV', maximumComptonEnergy=7, minimumPairEnergy=2, energyRangeCompton=None, phiRadiusCompton=15, openingAngleMax=60., energyHardCut = 5, parsing=True, events = None):

    """
    A function to plot the cosima output simulation file.
    Example Usage: 
    EventViewer.performCompleteAnalysis(filename='FarFieldPointSource_100MeV_Cos1.inc1.id1.tra')
    """

    if filename == None and directory == None:
        print("*** No filename or directory provide ***")
        print("Please provide a  filename, a list of filenames, or a directory name")
        return

    # Check to see if the user supplied a directory.  If so, include all .tra files in the directory
    if directory != None:
        old_dir = os.getcwd()
        os.chdir( directory )
        print ("Working in directory", directory )
        filenames = glob.glob( './*.tra') + glob.glob( './*.tra.gz')


    # Check if the user supplied a single file vs a list of files
    if isinstance(filename, list) == False and filename != None:
        filenames = [filename]


    # Try to get the energy from the filename
    if energies == None:
        energies = []
        for filename in filenames:
            try:
                filepath, filename = os.path.split(filename)                
                energy = float(filename.split('_')[1].replace(energySearchUnit,''))
                energies.append(energy)

            except:
                print("*** Unable to resolve energy from filename ***")
                print("Expected filename format: FarFieldPointSource_100MeV_Cos1.inc1.id1.tra")
                return

    # Try to get the angle from the filename
    if angles == None:
        angles = []
        for filename in filenames:
            try:
                filepath, filename = os.path.split(filename)                
                angle = float(filename.split('_')[2].split('.inc1')[0].replace('Cos',''))
                angles.append(angle)

            except:
                print("*** Unable to resolve angle from filename ***")
                print("Expected filename format: FarFieldPointSource_100MeV_Cos1.inc1.id1.tra")
                return

    # Print the identified files
    print("\nFiles identified for analysis:")
    for filename, energy, angle in zip(filenames, energies, angles):
        print("%s %s Cos %s %s" % (energy, energySearchUnit, angle, filename))
    print("")

    # Loop through the user specified filename(s) and extract the energy and angular resolution measurements
    for filename, energy, angle in zip(filenames, energies, angles):

        if parsing:
            print("Parsing: %s %s Cos %s %s" % (energy, energySearchUnit, angle, filename))
            # Parse the .tra file obtained from revan
            events = parse(filename, sourceTheta=angle)
            if not events:
                continue

        # Calculate the source theta in degrees
        source_theta = numpy.arccos(angle)*180./numpy.pi

        # Don't bother measuring the energy and angular resolutuon values for Compton events above the specified maximumComptonEnergy
        if energy <= maximumComptonEnergy:
            if energy >= 3:
                phiRadiusCompton = phiRadiusCompton/3.

            print("--------- All Compton Events ---------")
            # Calculate the energy resolution for Compton events
            print("Calculating the energy resolution for All Compton events...")
            print("EventAnalysis.getEnergyResolutionForComptonEvents(events, numberOfBins=100, energyPlotRange=None, energyFitRange=%s)" % (energyRangeCompton))
            mean, FWHM_energyComptonEvents, fitMax, FWHM_skewed_energyComptonEvents, sigma_Compton= getEnergyResolutionForComptonEvents(events, numberOfBins=100, onlyTrackedElectrons=False, onlyUntrackedElectrons=False, energyPlotRange=None, energyFitRange=energyRangeCompton, showPlots=showPlots, filename=filename, energyHardCut=energyHardCut, inputEnergy=energy)

            # Calculate the angular resolution measurement (ARM) for All Compton events
            print("\n\nCalculating the angular resolution measurement for Compton events...")
            print("EventAnalysis.getARMForComptonEvents(events, numberOfBins=100, phiRadius=%s)" % (phiRadiusCompton))
            FWHM_angleComptonEvents, dphi = getARMForComptonEvents(events, numberOfBins=100, phiRadius=phiRadiusCompton, onlyTrackedElectrons=False, onlyUntrackedElectrons=False, showPlots=showPlots, filename=filename, energyCutSelection = True, energyCut = [mean, sigma_Compton])

            # Calculate the energy resolution for tracked vs untracked Compton events
            print("--------- Untracked Compton Events ---------")
            print(" ")
            if energy <= 2:
                print("Calculating the energy resolution for Untracked Compton events...")
                print("EventAnalysis.getEnergyResolutionForComptonEvents(events, numberOfBins=100, energyPlotRange=None, energyFitRange=%s)" % (energyRangeCompton))
                mean_untracked, FWHM_energyUntrackedComptonEvents, UntrackedFitMax, FWHM_skewed_energyUntrackedComptonEvents, sigma_UntrackedCompton= getEnergyResolutionForComptonEvents(events, numberOfBins=100, onlyTrackedElectrons=False, onlyUntrackedElectrons=True, energyPlotRange=None, energyFitRange=energyRangeCompton, showPlots=showPlots, filename=filename, energyHardCut=energyHardCut, inputEnergy=energy)

                print("\n\nCalculating the angular resolution measurement for Untracked Compton events...")
                print("EventAnalysis.getARMForComptonEvents(events, numberOfBins=100, phiRadius=%s)" % (phiRadiusCompton))
                FWHM_angleUntrackedComptonEvents, dphi_untracked = getARMForComptonEvents(events, numberOfBins=100, phiRadius=phiRadiusCompton, onlyTrackedElectrons=False, onlyUntrackedElectrons=True, showPlots=showPlots, filename=filename, energyCutSelection = True, energyCut = [mean_untracked, sigma_UntrackedCompton])
            else:
                print("Energy too high (> 2 MeV) for untracked electrons. \n\n")
                mean_untracked = numpy.nan
                FWHM_energyUntrackedComptonEvents = numpy.nan
                UntrackedFitMax = numpy.nan
                FWHM_skewed_energyUntrackedComptonEvents = numpy.nan
                sigma_UntrackedCompton = numpy.nan
                FWHM_angleUntrackedComptonEvents = numpy.nan
                dphi_untracked = numpy.nan

            if energy>0.2:
                print("--------- Tracked Compton Events ---------")
                print("Calculating the energy resolution for Tracked Compton events...")
                print("EventAnalysis.getEnergyResolutionForComptonEvents(events, numberOfBins=100, energyPlotRange=None, energyFitRange=%s)" % (energyRangeCompton))
                mean_tracked, FWHM_energyTrackedComptonEvents, TrackedFitMax, FWHM_skewed_energyTrackedComptonEvents, sigma_TrackedCompton= getEnergyResolutionForComptonEvents(events, numberOfBins=100, onlyTrackedElectrons=True, onlyUntrackedElectrons=False, energyPlotRange=None, energyFitRange=energyRangeCompton, showPlots=showPlots, filename=filename, energyHardCut=energyHardCut, inputEnergy=energy)
                print("\n\nCalculating the angular resolution measurement for Tracked Compton events...")
                print("EventAnalysis.getARMForComptonEvents(events, numberOfBins=100, phiRadius=%s)" % (phiRadiusCompton))
                FWHM_angleTrackedComptonEvents, dphi_tracked = getARMForComptonEvents(events, numberOfBins=100, phiRadius=phiRadiusCompton, onlyTrackedElectrons=True, onlyUntrackedElectrons=False, showPlots=showPlots, filename=filename, energyCutSelection = True, energyCut = [mean_tracked, sigma_TrackedCompton])
            else:
                print("Energy too low (<0.2 MeV) for Tracked Compton events...")
                mean_tracked = numpy.nan
                TrackedFitMax = numpy.nan
                FWHM_energyTrackedComptonEvents = numpy.nan
                FWHM_skewed_energyTrackedComptonEvents = numpy.nan
                sigma_TrackedCompton = numpy.nan

        else:

            mean = numpy.nan
            FWHM_energyComptonEvents = numpy.nan
            FWHM_angleComptonEvents = numpy.nan
            sigma_Compton = numpy.nan
            fitMax = numpy.nan
            mean_tracked = numpy.nan
            FWHM_energyTrackedComptonEvents = numpy.nan
            FWHM_angleTrackedComptonEvents = numpy.nan
            sigma_TrackedCompton = numpy.nan
            mean_untracked = numpy.nan
            FWHM_energyUntrackedComptonEvents = numpy.nan
            FWHM_angleUntrackedComptonEvents = numpy.nan
            sigma_UntrackedCompton = numpy.nan
            fitMax = numpy.nan
            UntrackedFitMax = numpy.nan
            TrackedFitMax = numpy.nan

        # Don't bother measuring the energy and angular resolutuon values for pair events below the specified minimumPairEnergy
        if energy >= minimumPairEnergy:

            # Calculate the energy resolution for Pair events
            print("\n\nCalculating the energy resolution for pair events...")
            print("EventAnalysis.getEnergyResolutionForPairEvents(events, numberOfBins=100)")
            fileBase = "%sMeV_Cos%s" % (energy,angle)
            fitMax_pair, FWHM_pairComptonEvents, sigma_pair = getEnergyResolutionForPairEvents(events, numberOfBins=100, energyFitRange=True, showPlots=showPlots,fileBase=fileBase)

            # Calculate the angular resolution measurement (ARM) for pair events
            print("\n\nCalculating the angular resolution measurement for pair events...")
            print("EventAnalysis.getARMForPairEvents(events, numberOfBins=100, showDiagnosticPlots=False)")
            #angles, openingAngles, contaimentData_68, contaimentBinned_68, containmentFit_68, optimizedParameters = getARMForPairEvents(events, openingAngleMax=openingAngleMax, sourceTheta=source_theta, numberOfBins=500, showDiagnosticPlots=False, showPlots=showPlots, filename=filename, log=True, angleFitRange=[0,30], anglePlotRange=[30/500/10,30])
            angles, openingAngles, contaimentData_68, contaimentBinned_68, containmentFit_68, optimizedParameters = getARMForPairEvents(events, openingAngleMax=openingAngleMax, sourceTheta=source_theta, numberOfBins=500, showDiagnosticPlots=False, showPlots=showPlots, filename=filename, log=False, angleFitRange=[0,30], anglePlotRange=[-0.1,30])

        else:
            sigma_pair = numpy.nan
            fitMax_pair = numpy.nan
            contaimentData_68 = numpy.nan
            FWHM_pairComptonEvents = numpy.nan

        # Open the results filename for writing
        output_filename = filename.replace('.tra','.log')
        output = open(output_filename, 'w')

        # Write the results to disk
        if energy<0.2:
            sigma_TrackedCompton = numpy.nan
            FWHM_angleTrackedComptonEvents = numpy.nan
        if energy > 2:
            events['numberOfUntrackedElectronEvents'] = numpy.nan
        
        output.write("Results for simulation: %s %s Cos %s %s\n" % (energy, energySearchUnit, angle, filename))
        output.write("Compton Events Reconstructed: %s\n" % events['numberOfComptonEvents'])
        output.write("Compton Energy Resolution (keV): %s\n" % sigma_Compton) #FWHM_energyComptonEvents
        output.write("Compton Energy Mean (keV): %s\n" % mean)
        output.write("Compton Energy FitMax (keV): %s\n" % fitMax)
        output.write("Compton Angular Resolution (deg): %s\n" % FWHM_angleComptonEvents)
        output.write("Untracked Compton Events Reconstructed: %s\n" % events['numberOfUntrackedElectronEvents'])
        output.write("Untracked Compton Energy Mean (keV): %s\n" % mean_untracked)
        output.write("Untracked Compton Energy FitMax (keV): %s\n" % UntrackedFitMax)
        output.write("Untracked Compton Energy Resolution (keV): %s\n" % sigma_UntrackedCompton) #FWHM_energyUntrackedComptonEvents
        output.write("Untracked Compton Angular Resolution (deg): %s\n" % FWHM_angleUntrackedComptonEvents)
        output.write("Tracked Compton Events Reconstructed: %s\n" % events['numberOfTrackedElectronEvents'])
        output.write("Tracked Compton Energy Resolution (keV): %s\n" % sigma_TrackedCompton) #FWHM_energyTrackedComptonEvents
        output.write("Tracked Compton Energy Mean (keV): %s\n" % mean_tracked)
        output.write("Tracked Compton Energy FitMax (keV): %s\n" % TrackedFitMax)
        output.write("Tracked Compton Angular Resolution (deg): %s\n" % FWHM_angleTrackedComptonEvents)

        output.write("Pair Events Reconstructed: %s\n" % events['numberOfPairEvents'])
        output.write("Pair Energy Resolution (keV): %s\n" % sigma_pair) #FWHM_pairComptonEvents
        output.write("Pair Energy FitMax (keV): %s\n" % fitMax_pair)
        output.write("Pair Angular Containment (68%%): %s\n" % contaimentData_68)
        output.write("Pair Angular Resolution Fit Parameters: %s \n" % " ".join(map(str, optimizedParameters)))

        output.write("Events Not Reconstructed Flagged as Bad: %s\n" % events['numberOfBadEvents'])

        # Close the file
        output.close()

        print("Results saved to:\n%s" % output_filename)
        
        #try:
        #    clear_output(wait=True)
        #except:
        #    pass
    
    if directory is not None:
        os.chdir( old_dir )

    return events

##########################################################################################

def getTriggerEfficiency(filename=None, directory=None, save=True, savefile=None):
    """
    A function to extract the number of simulated events from a cosima .sim file
    Usage Examples: 
    EventViewer.getNumberOfSimulatedEvents(filename='FarFieldPointSource_100MeV_Cos1.inc1.id1.sim')
    EventViewer.getNumberOfSimulatedEvents(directory='./Simulations/MySimulations/')
    """

    if filename == None and directory == None:
        print("*** No filename or directory provide ***")
        print("Please provide a  filename, a list of filenames, or a directory name")
        return

    # Check to see if the user supplied a directory.  If so, include all .sim files in the directory
    if directory != None:
        print("\nSearching: %s\n" % directory)
        filenames = glob.glob(directory + '/*.sim' ) + glob.glob(directory + '/*.sim.gz' )

    # Check if the user supplied a single file vs a list of files
    if isinstance(filename, list) == False and filename != None:
        filenames = [filename]

    # Create a dictionary to return the results
    triggerEfficiency = {}      


    # Loop through each file
    for filename in filenames:

        print("Parsing: %s" % filename)

        # Use grep to find all lines with ID and pipe the results to the tail command to read only the last line
        #command = "grep 'ID ' %s | tail -n 1" % filename
        #output = os.popen(command).read()

        # Read backwards from the end of the file
        # Keep looking farther until you find an 'ID'
        # Much faster than the command above
        lookback = 1000
        IDs = []
        while len(IDs) == 0:
            if filename[-2:] == "gz": 
                command = "zcat %s | tail -n %d" % (filename, lookback) 
            else:
                command = "tail -n %d %s" % (lookback, filename)
            output = os.popen(command).read()
            IDs  = [line for line in output.split('\n') if "ID" in line]
            lookback += 1000
        output = IDs[-1]

        # Extract the number of triggers
        numberOfTriggers = float(output.split()[1])

        # Extract the number of simulated events
        numberOfSimulatedEvents = float(output.split()[2])

        # Add the values to the results dictionary
        triggerEfficiency[filename] = {
            'numberOfTriggers' : numberOfTriggers,
            'numberOfSimulatedEvents' : numberOfSimulatedEvents}

        # Add the details from the file.
        details = getDetailsFromFilename(filename)
        for key in details:
            triggerEfficiency[filename][key] = details[key]


    print("\nTrigger Efficiencies:")
    for filename in filenames:

        # Extract the values
        numberOfTriggers = triggerEfficiency[filename]['numberOfTriggers']
        numberOfSimulatedEvents = triggerEfficiency[filename]['numberOfSimulatedEvents']

        # Print the results
        print(filename, numberOfTriggers, numberOfSimulatedEvents, ' %.2f%%' % (100*numberOfTriggers/numberOfSimulatedEvents))

    # Check to see if the results should be written to disk
    if save == True:

        # Set a default filename if none was provided
        if savefile == None:
            savefile = 'TriggerEfficiency.txt'

        # Open the file for writing
        output = open(savefile, 'w')

        # Write the results to disk
        for filename in filenames:

            # Extract the values
            numberOfTriggers = triggerEfficiency[filename]['numberOfTriggers']
            numberOfSimulatedEvents = triggerEfficiency[filename]['numberOfSimulatedEvents']

            # Write out the values
            output.write("%s %s %s\n" % (filename, numberOfTriggers, numberOfSimulatedEvents))

        # Close the file
        output.close()

        print("\nResults saved to:\n./%s\n" % savefile)


    return triggerEfficiency

def getRevanTriggerEfficiency(filename=None, directory=None, save=True, savefile=None):

    """
    A function to extract the statistics from a revan tra file (used to check things)
    Usage Examples: 
    EventViewer.getNumberOfSimulatedEvents(filename='FarFieldPointSource_100MeV_Cos1.inc1.id1.sim')
    EventViewer.getNumberOfSimulatedEvents(directory='./Simulations/MySimulations/')
    """

    if filename == None and directory == None:
        print("*** No filename or directory provide ***")
        print("Please provide a  filename, a list of filenames, or a directory name")
        return

    # Check to see if the user supplied a directory.  If so, include all .tra files in the directory
    if directory != None:
        print("\nSearching: %s\n" % directory)
        filenames = glob.glob(directory + '/*.tra')

    # Check if the user supplied a single file vs a list of files
    if isinstance(filename, list) == False and filename != None:    
        filenames = [filename]

    # Create a dictionary to return the results
    revanStats = {}         

    # Loop through each file
    for filename in filenames:

        print("Parsing: %s" % filename)

        #Counter to look keep track of section dividers.
        counter = 0

        #Dictionary to save individual file results
        triggerStats = {}
        #Fill the details from the filname
        details = getDetailsFromFilename(filename)
        for key in details:
            triggerStats[key] = details[key]

            with open(filename) as infile:
                for line in infile:
                # Loop until you find these dividers.  Four dividers.
                    if line[0:5] == '-----':
                        counter += 1
            if counter == 2:
                split_line = line.split(':')
                if len(split_line) > 1:
                    if split_line[1] != '\n':
                        value = int(split_line[1].translate(None, '.'))
                        triggerStats[split_line[0].strip()] = value
            revanStats[filename] = triggerStats

    print("\nTrigger Efficiencies:")
    for filename in filenames:

        print(revanStats[filename])

    return revanStats





##########################################################################################


# if __name__ == "__main__":





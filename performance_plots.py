import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

from astropy.io import ascii
from operator import truediv

## The data file in ecsv format...
data = ascii.read('performance_plot_data.ecsv')

print "This is how many data points you have:", len(data)

def column(matrix, i):
    return [row[i] for row in matrix]

def make_plt(xdata,ydata,points, xmin, xmax, ymin, ymax, xlabel='',ylabel='',log=False,loglog=False):
    plt.figure()
    plt.plot(xdata,ydata)
    if log:
        plt.yscale('log')
    elif loglog:
        plt.yscale('log')
        plt.xscale('log')

    plt.plot(xdata, ydata, points)
    plt.gca().set_xlim(xmin,xmax)
    plt.gca().set_ylim(ymin,ymax)
    plt.gca().set_xlabel(xlabel)
    plt.gca().set_ylabel(ylabel)

#print column(data,'logEtrue_data')

logEtrue = column(data,'logEtrue_data')
Etrue = column(data,'Etrue_data')
Erec = column(data,'Ereco_data')
Eres = column(data,'Eres_data')
Ares = column(data,'Ares_data')
Aeff = column(data,'Aeff_data')

make_plt(logEtrue, Aeff,'ro-',2., 5.5, 1., 1000., r'log$_{10}$(Energy/keV)', r'A$_{\mathrm{eff}}$ [cm$^2$]', True, False)
plt.savefig('Aeff.png')
plt.show()
plt.close()

make_plt(logEtrue, Ares,'go-',2., 5.5, 0., 20., r'log$_{10}$(Energy/keV)', r'A$_{\mathrm{res}}$ [deg]', False, False)
plt.savefig('Ares.png')
plt.show()
plt.close()

make_plt(logEtrue, map(truediv, Eres, Etrue),'bo-',2., 5.5, 0., 0.2, r'log$_{10}$(Energy/keV)', r'E$_{\mathrm{res}}$ [E$_{\mathrm{width}}$/E$_{\mathrm{true}}$]', False, False)
plt.savefig('Eres.png')
plt.show()
plt.close()


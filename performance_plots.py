import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

from operator import truediv

## At some point should read in data files... 
logEtrue_data = [2.5, 3., 3.5, 4., 4.5, 5.]
Etrue_data = [316., 1000., 3160., 10000., 31600., 100000.]
Ereco_data = [314., 995., 3143., 9392., 26140., 68190.]
Eres_data = [12.24, 21.0, 48.0, 620.0, 3500., 9750.]
Ares_data = [8.0, 2.5, 1.9, 16.75, 6.25, 2.75]
Aeff_data = [25.4, 19.2, 7.2, 26.9, 133., 210.]

print "This is how many data points you have:", len(Aeff_data)

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

    plt.show()


make_plt(logEtrue_data, Aeff_data,'ro-',2., 5.5, 1., 1000., r'log$_{10}$(Energy/keV)', r'A$_{\mathrm{eff}}$ [cm$^2$]', True, False)
plt.savefig('Aeff.png')

make_plt(logEtrue_data, Ares_data,'go-',2., 5.5, 0., 20., r'log$_{10}$(Energy/keV)', r'A$_{\mathrm{res}}$ [deg]', False, False)
plt.savefig('Ares.png')

make_plt(logEtrue_data, map(truediv, Eres_data, Etrue_data),'bo-',2., 5.5, 0., 0.2, r'log$_{10}$(Energy/keV)', r'E$_{\mathrm{res}}$ [E$_width$/E$_true$]', False, False)
plt.savefig('Eres.png')


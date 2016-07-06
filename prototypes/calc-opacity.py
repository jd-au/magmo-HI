# Calculates the opacity for a specific location in a cube
# This uses the formula e^(-tau) = (T_on - T_off) / T_bg
# under the assumption that T_off is 0 (accurate for a long baseline)
# and that the target is a point source and thus T_on = flux and
# T_bg is the flux where there is no absorption

# Author James Dempsey
# Date 1 Jul 2016


patternClean='.*_clean$'
patternRestored='.*_restor$'

# for reading commandline args
import sys

#import pylab
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import os
import re
import subprocess
import math
import numpy as np
from astropy.table import Table, Column

class SourceStats:
    velocity = 0.0
    flux = 0.0

def getSpectrum(cube, region):
    """
    Retrieve the spectrum at a point in a miriad spectral cube

    :param cube: The Miriad forma spectral cube
    :param region: The definition on the region to be extracted
    :return: An array of data points, each with a plane id, velocity and flux
    """
    result = subprocess.check_output(["imspec", "in="+cube, "region="+region,
        "axes=ra,dec", "plot=mean", "options=guaranteespace,boxcar,6"])
    # imspec in=magmo-285.337-0.002_1420_sl_restor "region=boxes(267,198,267,198)"
    # axes=ra,dec plot=mean device=/xs options=boxcar,6,tb
    dataStarted = False
    spectrum = []
    zeroPoint = None
    for line in result.splitlines():
        # print line
        if dataStarted:
            # Skip first line and break once the data stops
            if line.strip().startswith("plane"):
                continue
            elif len(line.strip()) == 0:
                break
            # Store a data point
            stats = line.strip().split()
            #print stats
            point = (int(stats[0]), float(stats[1]), float(stats[2]))
            #print ">> " + str(point)
            if point[2] <= 0.0:
                # Trim leading and trailing zero values.
                zeroPoint = point
            else:
                if zeroPoint is not None and len(spectrum) > 0:
                    spectrum.append(zeroPoint)
                    zeroPoint = None
                spectrum.append(point)
        elif line is not None and line.find("VELO-") >= 0:
            dataStarted = True
    arr = np.array(spectrum, dtype=[('plane', int), ('velocity', float),
                                     ('flux', float)])
    spectrum_array = arr.view(np.recarray)
    return spectrum_array

def getMeanContinuum(spectrum):
    """
    Calulate the mean of the continuum values. Will divide the spectrum
    into continuum and line parts based on those values with a 1 sigma
    deviation from the mean of the first 5% of values.  This assumes that
    there is a leading continuum block in the spectrum.
    :param spectrum: The spectrum to be analysed, should be a numpy array of
        plane, velocity and flux values.
    :return: A single float whcih is the mean continuum flux.
    """
    fivePercent = len(spectrum) / 20
    firstFivePercRecs = spectrum.flux[:fivePercent-1]
    meanSeedFlux = np.mean(firstFivePercRecs)
    sdSeedFlux = np.std(firstFivePercRecs)
    #print str(firstFivePercRecs) + " " + str(meanSeedFlux) + " " + str(sdSeedFlux)

    continuumRecs = []
    lineRecs = []
    for point in spectrum:
        flux = point.flux
        if (abs(flux - meanSeedFlux) <= sdSeedFlux):
            continuumRecs.append(point[2])
        else:
            lineRecs.append(point[2])
    #print "Continuum [%d recs]: %s." % (len(continuumRecs), continuumRecs)
    #print "Line [%d recs]: %s." % (len(lineRecs), lineRecs)
    meanCont = np.mean(np.array(continuumRecs))
    #print meanCont
    return meanCont

def getOpacity(spectrum, mean):
    """
    Calculate the opacity profile for the spectrum. This simply divides the
    spectrum's flux by the mean.

    :param spectrum: The spectrum to be processed
    :param mean: The mean background flux, representing what the backlighting sources average flux.
    :return: The opacity (e^(-tau)) at each velocity step.
    """
    print spectrum.flux
    print spectrum.flux/mean
    return spectrum.flux/mean

def plotSpectrum(x, y, filename):
    """
    Output a plot of opacity vs LSR velocity to a specified file.

    :param x: The velocity data
    :param y: The opacity values for each velocity step
    :param filename: The file the plot should be written to. Should be an .eps or .pdf file.
    """
    fig = plt.figure()
    plt.plot(x, y)

    plt.axhline(1)

    plt.xlabel(r'Velocity relative to LSR (km/s)')
    plt.ylabel(r'$e^{(-\tau)}$')
    plt.grid(True)
    plt.savefig(filename)
    #plt.show()
    return

def outputSpectra(spectrum, opacity, filename):
    """
    Write the spectrum (velocity, flux and opacity) to a votable format file.

    :param spectrum: The spectrum to be output.
    :param opacity:  The opacity to be output.
    :param filename:  The filename to be created
    """
    table = Table([spectrum.plane, spectrum.velocity, opacity, spectrum.flux],
                  names=('plane', 'velocity', 'opacity', 'flux'),
                  meta={'name': 'Opacity', 'id' : 'opacity'})
    table.write(filename, format='votable', overwrite='True')



# Script begins here
if len(sys.argv) != 3:
    print("Incorrect number of parameters.")
    print("Usage: python calc-opacity.py cubename region")
    exit(1)

# spectrum = getSpectrum("magmo-285.337-0.002_1420_sl_restor", "boxes(267,198,267,198)")
spectrum = getSpectrum(sys.argv[1], sys.argv[2])
print spectrum
mean = getMeanContinuum(spectrum)
print mean
opacity = getOpacity(spectrum, mean)
print opacity
plotSpectrum(spectrum.velocity, opacity, "plot.pdf")
outputSpectra(spectrum, opacity, 'opacity.votable.xml')

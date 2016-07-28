# Plots a longitude velocity diagram for the MAGMO observations.
# This program will scan all votable files in the current directory
# extract the longitude, velocity and opacity data from each file
# and then plot it on a longitude-velocity diagram

# Author: James Dempsey

from astropy.io.votable import parse
import numpy as np
import matplotlib.pyplot as plt
import glob

def createSampleData():
    np.random.seed(0)
    n = 10000
    x = 360 * np.random.random(n) - 180
    y = np.random.random(n) * 100/180 * (180 - np.abs(x)) * (-x/abs(x))
    c = 1 - (np.random.random(n) * ((180 - np.abs(x)) / 180) )#* ((100 - np.abs(y))/100))
    return x, y, c

def readData():
    x = []
    y = []
    c = []

    voFiles = glob.glob('./*.votable.xml')
    for filename in voFiles:
        print 'Reading', filename
        votable = parse(filename, pedantic=False)
        results = next(resource for resource in votable.resources if resource.type == "results")
        if results is not None:
            for info in votable.infos:
                if info.name == 'longitude':
                    gal_info = info
            #gal_info = votable.get_info_by_id('longitude')
            if gal_info is None:
                print "No longitude provided for %s, skipping" , filename
                continue
            gal_long = int(float(gal_info.value))
            if gal_long > 180:
                gal_long -= 360
            results_array = results.tables[0].array
            for row in results_array:
                x.append(gal_long)
                y.append(row['velocity'])
                c.append(row['opacity'])
            #print results_array

    return x, y, c

def plot(x, y, c):
    xmin=-180
    xmax = 180
    ymin = -100
    ymax = 100

    # plt.subplots_adjust(hspace=0.5)
    plt.subplot(111, axisbg='black')
    plt.hexbin(x, y, c, cmap=plt.cm.gist_heat_r)
    # plt.scatter(x, y, cmap=plt.cm.YlOrRd_r)
    plt.axis([xmin, xmax, ymin, ymax])
    plt.title("Longitude-Velocity")
    plt.xlabel('Galactic longitude (l)')
    plt.ylabel('LSR Velocity (km/s)')
    cb = plt.colorbar()
    cb.set_label(r'$e^{(-\tau)}$')

    #plt.subplot(122)
    #plt.hexbin(x, y, c, bins='log', cmap=plt.cm.YlOrRd_r)
    #plt.axis([xmin, xmax, ymin, ymax])
    #plt.title("With a log color scale")
    #cb = plt.colorbar()
    #cb.set_label('log10(N)')

    plt.show()
    return

## Script starts here
#x, y, c = createSampleData()
x, y, c = readData()
plot(x, y, c)
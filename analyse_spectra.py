# Analyse the produced HI spectra and extract stats and produce diagrams.
#
# This program reads the previously generated spectra and stats files in the day
# folders and produces initial stage analysis data products. These include
# histograms of spectra produced vs used, and longitude velocity diagrams.

# Author James Dempsey
# Date 29 Sep 2016

from __future__ import print_function, division

from astropy.io.votable import parse
import numpy as np
import matplotlib.pyplot as plt
import glob
import time


class Field(object):

    def __init__(self, day, name):
        self.day = day
        self.name - name


class Spectrum(object):

    def __init__(self, day, field_name, src_id, longitude):
        self.day = day
        self.field_name = field_name
        self.src_id = src_id
        self.longitude = longitude


def read_data():
    """
    Read in the spectra produced in earlier pipeline stages, filter them for
     suitability and build up a set of longitude velocity data.

    :return: Three arrays, the x, y and opacity values for each data point in
     the usable spectra
    """
    x = []
    y = []
    c = []

    bad_spectra = 0
    prev_field = ''
    num_fields = 0
    vo_files = glob.glob('day*/*.votable.xml')
    for filename in sorted(vo_files):
        print ('Reading', filename)
        votable = parse(filename, pedantic=False)
        results = next(resource for resource in votable.resources if
                       resource.type == "results")
        if results is not None:
            gal_info = None
            for info in votable.infos:
                if info.name == 'longitude':
                    gal_info = info
            #gal_info = votable.get_info_by_id('longitude')
            if gal_info is None:
                print("No longitude provided for %s, skipping" % filename)
                continue
            gal_long = float(gal_info.value)
            if gal_long > 180:
                gal_long -= 360
            results_array = results.tables[0].array
            poor_sn = False
            #print gal_long
            for row in results_array:
                opacity = row['opacity']
                if opacity > 6 or opacity < -8:
                    poor_sn = True
            if not poor_sn:
                for row in results_array:
                    x.append(gal_long)
                    y.append(row['velocity']/1000.0) # Convert from m/s to km/s
                    opacity = row['opacity']
                    c.append(opacity)
                    if opacity > 10 or opacity < -15:
                        poor_sn = True
                field = filename.split('_')
                if field[0] != prev_field:
                    prev_field = field[0]
                    num_fields += 1
            if poor_sn:
                bad_spectra += 1

    print ("In %d fields read %d spectra of which %d had reasonable S/N. " % (
        num_fields, len(vo_files), len(vo_files)-bad_spectra))
    return x, y, c


def plot_lv(x, y, c, filename):
    """
    Produce a longitude-velocity diagram from the supplied data and write it
    out to a file.

    :param x: The x-coordinate of each point (galactic longitude in degrees)
    :param y: The y-coordinate of each point (LSR velocity in km/s)
    :param c: The opacity fraction (1= transparent, 0=completely opaque)
    :param filename: The file to write the plot to.
    :return: None
    """
    xmin =-180
    xmax = 180
    ymin = -300
    ymax = 300

    val = np.clip(c, -0.005, 1.05)
    #print val
    fig = plt.figure(1, (12,6))
    # plt.subplots_adjust(hspace=0.5)
    plt.subplot(111, axisbg='black')
    plt.hexbin(x, y, val, cmap=plt.cm.gist_heat_r)
    # plt.scatter(x, y, cmap=plt.cm.YlOrRd_r)
    plt.axis([xmax, xmin, ymin, ymax])
    plt.title("Longitude-Velocity")
    plt.xlabel('Galactic longitude (l)')
    plt.ylabel('LSR Velocity (km/s)')
    cb = plt.colorbar()
    cb.set_label(r'$e^{(-\tau)}$')

    plt.savefig(filename)
    return


def main():
    """
    Main script for analyse_spectra
    :return: The exit code
    """
    start = time.time()

    print("#### Started analysis on MAGMO spectra at %s ####" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    # Process
    x, y, c = read_data()
    plot_lv(x, y, c, 'magmo-lv.pdf')
    # also want
    # - Catalogue - Fields - day, field, peak, sn, coords, used
    # - Catalogue - Source - field, source, continuum, min, max, sn, used
    # - Histogram - fields observed/used per day
    # - Histogram - fields observed/used per day

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Processed spectra in %.02f s' %
          (end - start))
    return 0


# Run the script if it is called from the command line
if __name__ == "__main__":
    exit(main())


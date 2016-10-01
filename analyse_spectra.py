# Analyse the produced HI spectra and extract stats and produce diagrams.
#
# This program reads the previously generated spectra and stats files in the day
# folders and produces initial stage analysis data products. These include
# histograms of spectra produced vs used, and longitude velocity diagrams.

# Author James Dempsey
# Date 29 Sep 2016

from __future__ import print_function, division

from astropy.io.votable import parse, from_table, writeto
from astropy.table import Table

import csv
import glob
import matplotlib.pyplot as plt
import numpy as np
import time


class Field(object):
    """
    "name", "rms", "max", "sn", "strong"
    """

    def __init__(self, day, name, rms, max_flux, sn_ratio, used, longitude,
                 latitude):
        self.day = day
        self.name = name
        self.rms = rms
        self.max_flux = max_flux
        self.sn_ratio = sn_ratio
        self.used = used

        self.longitude = longitude
        self.latitude = latitude


class Spectrum(object):
    """

    """

    def __init__(self, day, field_name, src_id, longitude, latitude, velocities,
                 opacities, fluxes):
        self.day = day
        self.field_name = field_name
        self.src_id = src_id
        self.longitude = longitude
        self.latitude = latitude
        self.velocities = velocities
        self.opacities = opacities
        self.fluxes = fluxes
        self.low_sn = None

    def get_field_id(self):
        return str(self.day) + "-" + str(self.field_name)

    def __str__(self):
        return self.get_field_id() + ", src: " + str(self.src_id)


def read_spectra():
    """
    Read in the spectra produced in earlier pipeline stages.

    :return: An array of Spectrum objects
    """
    spectra = []

    vo_files = glob.glob('day*/*.votable.xml')
    for filename in sorted(vo_files):
        print ('Reading', filename)
        votable = parse(filename, pedantic=False)
        results = next(resource for resource in votable.resources if
                       resource.type == "results")
        if results is not None:
            gal_long = None
            gal_lat = None
            for info in votable.infos:
                if info.name == 'longitude':
                    gal_long = float(info.value)
                    if gal_long > 180:
                        gal_long -= 360
                if info.name == 'latitude':
                    gal_lat = float(info.value)
            if gal_long is None:
                print("No longitude provided for %s, skipping" % filename)
                continue
            results_array = results.tables[0].array

            velocities = np.zeros(len(results_array))
            opacities = np.zeros(len(results_array))
            fluxes = np.zeros(len(results_array))
            i = 0
            for row in results_array:
                opacity = row['opacity']
                velocities[i] = row['velocity']/1000.0
                opacities[i] = opacity
                fluxes[i] = row['flux']
                i += 1
            field = filename.split('_')
            parts = field[0].split('/')
            spectrum = Spectrum(str(parts[0][3:]), parts[1], field[1][3:],
                                gal_long, gal_lat, velocities, opacities,
                                fluxes)
            spectra.append(spectrum)

    return spectra


def read_field_stats():
    """
    Read in all of the stats csv files for the days processed and build a list
    of fields.

    :return: List of fields.
    """

    fields = []

    stats_files = glob.glob('day*/stats.csv')
    for filename in sorted(stats_files):
        print ('Reading', filename)
        day = int(filename.split('/')[0][3:])

        with open(filename, 'rb') as stats:
            reader = csv.reader(stats)
            first = True
            for row in reader:
                if first:
                    first = False
                else:
                    if "-" in row[0]:
                        coords = row[0].split('-')
                        coords[1] = "-"+coords[1]
                    else:
                        coords = row[0].split('+')

                    used = 'N'
                    if len(row) > 4:
                        used = row[4]
                    field = Field(day, row[0], row[1], row[2], row[3], used,
                                  float(coords[0]), float(coords[1]))
                    fields.append(field)

    return fields


def extract_lv(spectra):
    x = []
    y = []
    c = []
    bad_spectra = 0
    prev_field = ''
    num_fields = 0

    for spectrum in spectra:
        opacities = spectrum.opacities
        if np.max(opacities) > 6 or np.min(opacities) < -8:
            bad_spectra += 1
            spectrum.low_sn = True
            continue
        y = np.concatenate((y, spectrum.velocities))
        c = np.concatenate((c, opacities))
        for row in opacities:
            x.append(spectrum.longitude)
        field_id = spectrum.get_field_id()
        if field_id != prev_field:
            prev_field = field_id
            num_fields += 1

    print ("In %d fields read %d spectra of which %d had reasonable S/N. " % (
        num_fields, len(spectra), len(spectra)-bad_spectra))

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

    print("X: %d, Y: %d, data: %d" % (len(x), len(y), len(c) ))
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


def output_spectra_catalogue(spectra):
    """
    Output the list of spectrum stats to a VOTable file magmo-spectra.vot

    :param spectra: The list of Spectrum objects
    :return: None
    """
    rows = len(spectra)
    days = np.zeros(rows, dtype=int)
    fields = np.empty(rows, dtype=object)
    sources = np.empty(rows, dtype=object)
    longitudes = np.zeros(rows)
    latitudes = np.zeros(rows)
    max_flux = np.zeros(rows)
    max_opacity = np.zeros(rows)
    min_opacity = np.zeros(rows)
    used = np.empty(rows, dtype=bool)

    i = 0
    for spectrum in spectra:
        days[i] = int(spectrum.day)
        fields[i] = spectrum.field_name
        sources[i] = spectrum.src_id
        longitudes[i] = spectrum.longitude
        latitudes[i] = spectrum.latitude
        max_flux[i] = np.max(spectrum.fluxes)
        max_opacity[i] = np.max(spectrum.opacities)
        min_opacity[i] = np.min(spectrum.opacities)
        used[i] = True if not spectrum.low_sn else False
        i += 1

    spectra_table = Table(
        [days, fields, sources, longitudes, latitudes, max_flux, max_opacity,
         min_opacity, used],
        names=['Day', 'Field', 'Source', 'Longitude',
               'Latitude', 'Max_Flux', 'Max_Opacity', 'Min_Opacity', 'Used'])
    votable = from_table(spectra_table)
    filename = "magmo-spectra.vot"
    writeto(votable, filename)
    print("Wrote out", i, "spectra to", filename)


def output_field_catalogue(fields):
    rows = len(fields)
    days = np.zeros(rows, dtype=int)
    field_names = np.empty(rows, dtype=object)
    longitudes = np.zeros(rows)
    latitudes = np.zeros(rows)
    max_fluxes = np.zeros(rows)
    sn_ratios = np.zeros(rows)
    used = np.empty(rows, dtype=bool)

    i = 0
    for field in fields:
        days[i] = int(field.day)
        field_names[i] = field.name
        longitudes[i] = field.longitude
        latitudes[i] = field.latitude
        max_fluxes[i] = field.max_flux
        sn_ratios[i] = field.sn_ratio
        sn_ratios[i] = field.sn_ratio
        used[i] = True if not field.used == 'Y' else False
        i += 1

    fields_table = Table(
        [days, field_names, longitudes, latitudes, max_fluxes, sn_ratios, used],
        names=['Day', 'Field', 'Longitude',
               'Latitude', 'Max_Flux', 'SN_Ratio', 'Used'])
    votable = from_table(fields_table)
    filename = "magmo-fields.vot"
    writeto(votable, filename)

    print("Wrote out", i, "fields to", filename)


def main():
    """
    Main script for analyse_spectra
    :return: The exit code
    """
    start = time.time()

    print("#### Started analysis on MAGMO spectra at %s ####" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    # Process Spectra
    spectra = read_spectra()
    x, y, c = extract_lv(spectra)
    plot_lv(x, y, c, 'magmo-lv.pdf')
    output_spectra_catalogue(spectra)

    # Process Fields
    fields = read_field_stats()
    output_field_catalogue(fields)

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


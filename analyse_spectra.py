#!/usr/bin/env python -u

# Analyse the produced HI spectra and extract stats and produce diagrams.
#
# This program reads the previously generated spectra and stats files in the day
# folders and produces initial stage analysis data products. These include
# histograms of spectra produced vs used, and longitude velocity diagrams.

# Author James Dempsey
# Date 29 Sep 2016

from __future__ import print_function, division

from astropy.io import fits
from astropy.io.votable import parse, from_table, writeto
from astropy.table import Table
from scipy import ndimage
import analyse_data
import csv
import datetime
import glob
import magmo
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
import os
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

    def get_field_id(self):
        return str(self.day) + "-" + str(self.name)


class Spectrum(object):
    """

    """

    def __init__(self, day, field_name, src_id, longitude, latitude, velocity,
                 opacities, flux):
        self.day = day
        self.field_name = field_name
        self.src_id = src_id
        self.longitude = longitude
        self.latitude = latitude
        self.velocity = velocity
        self.opacities = opacities
        self.flux = flux
        self.low_sn = None
        self.range = 0
        self.opacity_range = 0
        self.max_s_max_n = 0
        self.continuum_sd = 0
        self.rating = 'A'

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

    vo_files = glob.glob('day*/*_opacity.votable.xml')
    print("Reading {} spectrum files.".format(len(vo_files)))
    for filename in sorted(vo_files):
        # print ('Reading', filename)
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
                velocities[i] = row['velocity'] / 1000.0
                opacities[i] = opacity
                fluxes[i] = row['flux']
                i += 1
            field = filename.split('_')
            parts = field[0].split('/')
            spectrum = Spectrum(str(parts[0][3:]), parts[1], field[1][3:],
                                gal_long, gal_lat, velocities, opacities,
                                fluxes)
            min_opacity = np.min(spectrum.opacities)
            max_opacity = np.max(spectrum.opacities)

            continuum_ranges = magmo.get_continuum_ranges()

            opacity_range = max_opacity - min_opacity
            max_s_max_n = (1 - min_opacity) / (max_opacity - 1)
            continuum_sd = calc_continuum_sd(spectrum, continuum_ranges)
            rating = calc_rating(opacity_range, max_s_max_n, continuum_sd)

            spectrum.opacity_range = opacity_range
            spectrum.max_s_max_n = max_s_max_n
            spectrum.continuum_sd = continuum_sd
            spectrum.rating = rating
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
    print("Reading {} day stats files.".format(len(stats_files)))
    for filename in sorted(stats_files):
        # print ('Reading', filename)
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
                        coords[1] = "-" + coords[1]
                    else:
                        coords = row[0].split('+')

                    used = 'N'
                    if len(row) > 4:
                        used = row[4]
                    field = Field(day, row[0], row[1], row[2], row[3], used,
                                  float(coords[0]), float(coords[1]))
                    fields.append(field)

    return fields


def extract_lv(spectra, min_rating='C'):
    x = []
    y = []
    c = []
    bad_spectra = 0
    prev_field = ''
    num_fields = 0
    used_fields = []

    for spectrum in spectra:
        opacities = spectrum.opacities
        if spectrum.rating > min_rating:
            bad_spectra += 1
            spectrum.low_sn = True
            continue
        #if np.max(opacities) > max_opacity or np.min(opacities) < min_opacity:
        #    bad_spectra += 1
        #   spectrum.low_sn = True
        #    continue
        y = np.concatenate((y, spectrum.velocity))
        c = np.concatenate((c, opacities))
        x = np.concatenate((x, np.full(len(opacities), spectrum.longitude)))
        field_id = spectrum.get_field_id()
        if field_id != prev_field:
            prev_field = field_id
            used_fields.append(field_id)
            num_fields += 1

    print("In %d fields read %d spectra of which %d had reasonable S/N. " % (
        num_fields, len(spectra), len(spectra) - bad_spectra))

    return x, y, c, used_fields


def plot_lv(x, y, c, filename, continuum_ranges, zoom):
    """
    Produce a longitude-velocity diagram from the supplied data and write it
    out to a file.

    :param x: The x-coordinate of each point (galactic longitude in degrees)
    :param y: The y-coordinate of each point (LSR velocity in km/s)
    :param c: The opacity fraction (1= transparent, 0=completely opaque)
    :param filename: The file to write the plot to.
    :param continuum_ranges: The file to write the plot to.
    :param zoom: Should the plot be zoomed in on the data region
    ,
    :return: None
    """
    xmin = -120 if zoom else -180
    xmax = 30 if zoom else 180
    ymin = -300
    ymax = 300

    # print("X: %d, Y: %d, data: %d" % (len(x), len(y), len(c) ))
    val = np.clip(c, -0.005, 1.05)

    fig_size = plt.rcParams["figure.figsize"]
    fig_size[1] = 4.5
    plt.rcParams["figure.figsize"] = fig_size

    # print val
    fig = plt.figure(1, (12, 6))
    # plt.subplots_adjust(hspace=0.5)
    plt.subplot(111, axisbg='black' if zoom else 'gray')
    plt.hexbin(x, y, val, cmap=plt.cm.gist_heat_r)
    # plt.scatter(x, y, cmap=plt.cm.YlOrRd_r)
    plt.axis([xmax, xmin, ymin, ymax])
    plt.title("Longitude-Velocity")
    plt.xlabel('Galactic longitude (deg)')
    plt.ylabel('LSR Velocity (km/s)')
    cb = plt.colorbar(orientation='horizontal')
    #cb = plt.colorbar()
    cb.set_label(r'$e^{(-\tau)}$')

    # Add bands for the continuum ranges
    if not zoom:
        for con_range in continuum_ranges:
            min_l = con_range['min_long']
            max_l = con_range['max_long']
            if min_l < 181:
                min_x = 0.5 - (min_l / 360.0)
                max_x = 0.5 - (max_l / 360.0)
            else:
                min_x = 0.5 + ((360 - min_l) / 360.0)
                max_x = 0.5 + ((360 - max_l) / 360.0)
            plt.axhline(con_range['min_con_vel'], xmin=min_x,
                        xmax=max_x, color='blue')
            # linestyle='dashed')
            plt.axhline(con_range['max_con_vel'], xmin=min_x,
                        xmax=max_x, color='blue')
            # linestyle='dashed')

    plt.grid(color='White')

    plt.savefig(filename)
    plt.close()
    print("Plotted ", len(c), "opacity points to", filename)

    return


def world_to_pixel(header, axis, value):
    """
    Calculate the pixel value for the provided world value using the WCS
    keywords on the specific axis. The axis must be linear.
    :param header: The FITS header describing the zxes.
    :param axis:  The number of the target axis.
    :param value: The world value to be converted.
    :return: The pixel value.
    """
    ax = str(axis)
    return int(header['CRPIX' + ax] + (value - header['CRVAL' + ax]) / header[
        'CDELT' + ax])


def get_lv_subset(data, header, l_min, l_max, v_min, v_max):
    """
    Extract a subset of velocity, longitude data based on physical bounds.

    :param data: The two dimensional array of data.
    :param header: The FITS header of the data with axes of longitude, velocity.
    :param l_min: The minimum of the desired longitude range.
    :param l_max: The maximum of the desired longitude range.
    :param v_min: The minimum of the desired velocity range.
    :param v_max: The maximum of the desired velocity range.
    :return: A numpy array with only the data from the requested range.
    """

    l_start = world_to_pixel(header, 1, l_max)
    l_end = world_to_pixel(header, 1, l_min)
    v_start = world_to_pixel(header, 2, v_min)
    v_end = world_to_pixel(header, 2, v_max)

    return data[v_start:v_end, l_start:l_end]


def plot_lv_image(x, y, c, filename):
    """
    Output a longitude-velocity plot of the provided data with the outline of
    emission from the GASS III dataset plotted over the data.

    :param x: The longitude value of each data point.
    :param y: The velocity value of each data point.
    :param c: The optical depth value of each data point.
    :param filename: The file name of the plot.
    :return: None
    """
    # Image dimensions
    l_max = 30
    l_min = -120
    l_dpd = 1 / 0.08
    l_size = int((l_max - l_min) * l_dpd)
    v_max = 300
    v_min = -300
    v_dpkms = 1 / 0.8245
    v_size = int((v_max - v_min) * v_dpkms)

    val = np.clip(c, -0.005, 1.05)

    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'

    dots_per_degree = l_dpd  # 4*3
    data = ma.array(np.ones((v_size, l_size)), mask=True)
    # print(data)
    xmax = data.shape[1]
    ymax = data.shape[0]
    for i in range(0, len(x)):
        x_idx = xmax - int((x[i] - l_min) * dots_per_degree)
        y_idx = ymax - int((y[i] - v_min) * v_dpkms)
        data[y_idx, x_idx - 3:x_idx + 4] = val[i]

    fig_size = plt.rcParams["figure.figsize"]
    fig_size[1] = 4.5
    plt.rcParams["figure.figsize"] = fig_size

    smoothed_data = ndimage.gaussian_filter(data, sigma=2, order=0)

    ax = plt.subplot(111)
    img = ax.imshow(smoothed_data, cmap=plt.cm.gist_heat_r)
    plt.title("Longitude-Velocity")
    plt.xlabel('Galactic longitude (deg)')
    plt.ylabel('LSR Velocity (km/s)')

    #cbaxes = plt.add_axes([0.05, 0.05, 0.9, 0.025])
    #cb = plt.colorbar(cax=cbaxes, mappable=mappable, orientation='horizontal')
    cb = plt.colorbar(img, ax=ax, orientation='horizontal')

    #cb = plt.colorbar(img, ax=ax)
    cb.set_label(r'$e^{(-\tau)}$')

    gass_lv = fits.open('gass-lv.fits')
    gass_subset = get_lv_subset(gass_lv[0].data, gass_lv[0].header, l_min,
                                l_max, v_min * 1000,
                                v_max * 1000)

    # Add an outline of the emission from GASS III
    plt.contour(np.log10(np.flipud(gass_subset)), 1, cmap='Pastel2')

    # Set the axis ticks and scales
    x_step = int(20 * (gass_subset.shape[1] / (l_max - l_min)))
    ax.set_xticks([i for i in range(gass_subset.shape[1], 0, -x_step)])
    ax.set_xticklabels([i for i in range(l_min, l_max + 1, 20)])

    y_step = int(100 * (gass_subset.shape[0] / (v_max - v_min)))
    ax.set_yticks([i for i in range(0, gass_subset.shape[0], y_step)])
    ax.set_yticklabels([i for i in range(v_max, v_min - 1, -100)])

    plt.grid(color='antiquewhite')

    plt.savefig(filename)
    plt.close()


def calc_continuum_sd(spectrum, continuum_ranges):
    """
    Calulate the standard deviaition of opacity in the velocity range
    designated as continuum for the spectrum's longitude. This gives a measure
    of the noise in wat should be an otherwise continuum only part of the
    spectrum.

    :param spectrum: The spectrum object being analysed.
    :param continuum_ranges: The defined contionuum ranges
    :return: The opacity standard deviation.
    """

    continuum_start_vel, continuum_end_vel = magmo.lookup_continuum_range(
        continuum_ranges, int(spectrum.longitude))
    continuum_range = np.where(
        continuum_start_vel < spectrum.velocity)
    bin_start = continuum_range[0][0]
    continuum_range = np.where(
        spectrum.velocity < continuum_end_vel)
    bin_end = continuum_range[0][-1]
    sd_cont = np.std(spectrum.opacities[bin_start:bin_end])
    return sd_cont


def calc_rating(opacity_range, max_s_max_n, continuum_sd):
    rating_codes = 'ABCDEF'
    rating = 0

    if opacity_range > 1.5:
        rating += 1
    if max_s_max_n < 3:
        rating += 1
    if continuum_sd*3 > 1:
        rating += 1

    return rating_codes[rating]


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
    max_velocity = np.zeros(rows)
    min_velocity = np.zeros(rows)
    rms_opacity = np.zeros(rows)
    opacity_range = np.zeros(rows)
    continuum_sd = np.zeros(rows)
    max_s_max_n = np.zeros(rows)
    rating = np.empty(rows, dtype=object)
    used = np.empty(rows, dtype=bool)
    filenames = np.empty(rows, dtype=object)
    local_paths = np.empty(rows, dtype=object)

    base_path = os.path.realpath('.')
    i = 0
    for spectrum in spectra:
        days[i] = int(spectrum.day)
        fields[i] = spectrum.field_name
        sources[i] = spectrum.src_id
        longitudes[i] = spectrum.longitude
        latitudes[i] = spectrum.latitude
        max_flux[i] = np.max(spectrum.flux)
        min_opacity[i] = np.min(spectrum.opacities)
        max_opacity[i] = np.max(spectrum.opacities)
        rms_opacity[i] = np.sqrt(np.mean(np.square(spectrum.opacities)))
        min_velocity[i] = np.min(spectrum.velocity)
        max_velocity[i] = np.max(spectrum.velocity)

        opacity_range[i] = spectrum.opacity_range
        max_s_max_n[i] = spectrum.max_s_max_n
        continuum_sd[i] = spectrum.continuum_sd
        rating[i] = spectrum.rating

        used[i] = not spectrum.low_sn
        filenames[i] = 'day' + spectrum.day + '/' + spectrum.field_name + \
                       "_src" + spectrum.src_id + "_plot.png"
        local_paths[i] = base_path + '/' + filenames[i]
        i += 1

    spectra_table = Table(
        [days, fields, sources, longitudes, latitudes, max_flux, min_opacity,
         max_opacity, rms_opacity, min_velocity, max_velocity, used,
         opacity_range, max_s_max_n, continuum_sd, rating,
         filenames, local_paths],
        names=['Day', 'Field', 'Source', 'Longitude', 'Latitude', 'Max_Flux',
               'Min_Opacity', 'Max_Opacity', 'RMS_Opacity', 'Min_Velocity',
               'Max_Velocity', 'Used', 'Opacity_Range', 'Max_S_Max_N',
               'Continuum_SD', 'Rating', 'Filename', 'Local_Path'],
        meta={'ID': 'magmo_spectra',
              'name': 'MAGMO Spectra ' + str(datetime.date.today())})
    votable = from_table(spectra_table)
    filename = "magmo-spectra.vot"
    writeto(votable, filename)
    print("Wrote out", i, "spectra to", filename)
    for grade in "ABCDEF":
        num_rating = len(np.where(rating == grade)[0])
        print ("%s: %3d" % (grade, num_rating))
    print ("Mean continuum sd %f" % np.mean(continuum_sd))


def output_field_catalogue(fields, used_fields):
    """
    Write out a catalogue of the fields observed under the MAGMO project
    with some basic stats for each field.

    :param fields: The fields to be written.
    :param used_fields: An aray of field ids which had spectra which were used.
    :return: None
    """
    rows = len(fields)
    days = np.zeros(rows, dtype=int)
    field_names = np.empty(rows, dtype=object)
    longitudes = np.zeros(rows)
    latitudes = np.zeros(rows)
    max_fluxes = np.zeros(rows)
    sn_ratios = np.zeros(rows)
    strong = np.empty(rows, dtype=bool)
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
        strong[i] = True if field.used == 'Y' else False
        used[i] = field.get_field_id() in used_fields
        i += 1

    fields_table = Table(
        [days, field_names, longitudes, latitudes, max_fluxes, sn_ratios,
         strong, used],
        names=['Day', 'Field', 'Longitude',
               'Latitude', 'Max_Flux', 'SN_Ratio', 'Strong', 'Used'],
        meta={'ID': 'magmo_fields',
              'name': 'MAGMO Fields ' + str(datetime.date.today())})
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

    print("#### Started analysis of MAGMO spectra at %s ####" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    # Process Spectra
    spectra = read_spectra()
    x, y, c, used_fields = extract_lv(spectra)
    continuum_ranges = magmo.get_continuum_ranges()
    plot_lv(x, y, c, 'magmo-lv.pdf', continuum_ranges, False)
    plot_lv(x, y, c, 'magmo-lv-zoom.pdf', continuum_ranges, True)
    plot_lv_image(x, y, c, 'magmo-lv-zoom-im.pdf')
    output_spectra_catalogue(spectra)

    # Output only the really good spectra
    x, y, c, temp = extract_lv(spectra, min_rating='B')
    plot_lv(x, y, c, 'magmo-lv_AB.pdf', continuum_ranges, False)
    plot_lv_image(x, y, c, 'magmo-lv-zoom-im_AB.pdf')

    # Process Fields
    fields = read_field_stats()
    output_field_catalogue(fields, used_fields)

    # also want
    # - Catalogue - Fields - day, field, peak, sn, coords, used
    # - Catalogue - Source - field, source, continuum, min, max, sn, used
    # - Histogram - fields observed/used per day
    # - Histogram - fields observed/used per day

    # Report
    end = time.time()
    print('#### Analysis completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Processed spectra in %.02f s' %
          (end - start))
    return 0


# Run the script if it is called from the command line
if __name__ == "__main__":
    exit(main())

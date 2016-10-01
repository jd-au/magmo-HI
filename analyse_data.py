# Find sources in the data and produce spectra for each suitable source.

# Author James Dempsey
# Date 28 Aug 2016

from __future__ import print_function, division

import magmo
import os
import sys
import time
import csv

from astropy.io import fits
from astropy.io import votable
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs import WCS
from astropy.table import Table, Column
from astropy.io.votable.tree import Param,Info
from astropy.io.votable import from_table, writeto
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
import numpy.core.records as rec

from string import Template


sn_min = 1.3
num_chan = 627


def get_high_signal_fields(day_dir_name):
    """
    Retrieve a list of fields observed in a particular day that have sufficient
    signal to noise to search for background sources.
    :param day_dir_name: The name of the day's directory.
    :return: A list of high signal fields.
    """
    field_list = []
    print ("Fields of interest:")
    with open(day_dir_name + '/stats.csv', 'rb') as stats:
        reader = csv.reader(stats)
        first = True
        for row in reader:
            if first:
                first = False
            else:
                if float(row[3]) > sn_min:
                    print (row)
                    field_list.append(row[0])

    return field_list


def find_sources(day_dir_name, field_name):
    """
    Search a contoinuum file for sources using the Aegean source finder. A
    VOTable file containing the list of discovered sources will be written out
    for the field. This function will use the Aegean source finder
    ( https://github.com/PaulHancock/Aegean ) to identify the sources.

    :param day_dir_name: The name of the day's directory.
    :param field_name:  The name fo the field to be searched for sources.
    :return: A list of error messages, if any
    """
    error_list = []
    cont_file = day_dir_name + "/1757/magmo-" + field_name + "_1757_restor.fits"
    table_file = day_dir_name + "/" + field_name + '_src.vot'
    try:
        print ("##--## Searching continuum image " + cont_file + " ##--##")
        magmo.run_os_cmd('bane ' + cont_file)
        aegean_cmd = 'aegean ' + cont_file + ' --autoload --telescope=ATCA ' \
                     '--cores=1 --table=' + table_file
        magmo.run_os_cmd(aegean_cmd)
    except magmo.CommandFailedError as e:
        error_list.append(str(e))
    return error_list


def read_sources(filename):
    print ("Extracting sources from " + filename)
    sources = []
    src_votable = votable.parse(filename, pedantic=False)
    results = src_votable.get_first_table().array
    for row in results:
        id = str(row['island']) + "-" + str(row['source'])
        ra = row['ra']
        dec = row['dec']
        rms = row['local_rms']
        flux = row['peak_flux']
        sn = flux / rms
        print ("Found source %s at %.4f, %.4f with flux %.4f and rms of %.4f "
               "giving S/N of %.4f" % (id, ra, dec, flux, rms, sn))
        if sn > 10 and flux > 0.02:
            sources.append([ra, dec, id, flux])
        else:
            print ("Ignoring source at %.4f, %.4f due to low S/N of %.4f or "
                   "flux of %.4f" % (ra, dec, sn, flux))

    return sources


def read_continuum_ranges():
    continuum_ranges = []
    with open('magmo-continuum.csv', 'rb') as con_def:
        reader = csv.reader(con_def)
        first = True
        for row in reader:
            if first:
                first = False
            else:
                continuum_ranges.append(
                    [int(row[0]), int(row[1]), int(row[2]), int(row[3])])

    print (continuum_ranges)
    return continuum_ranges


def find_edges(fluxes, num_edge_chan):
    """
    Seek from the edges to find where the data starts for this set of fluxes.
    This accounts for an optional number of channels in the data which have no
    data recorded.
    :param fluxes: The array of fluxes to be checked.
    :param num_edge_chan: The number of edge channels with data to be skipped
    :return: The index of the first and last cell to have data.
    """

    l_edge = 0
    r_edge = len(fluxes)-1

    while fluxes[l_edge] == 0 and l_edge < len(fluxes):
        l_edge += 1

    while fluxes[r_edge] == 0 and r_edge > 0:
        r_edge -= 1

    return l_edge + num_edge_chan, r_edge - num_edge_chan


def extract_spectra(daydirname, field):
    num_edge_chan = 10  # todo: make this a velocity relative value
    fits_filename = "{0}/1420/magmo-{1}_1420_sl_restor.fits".format(daydirname,
                                                                    field)
    src_filename = "{0}/{1}_src_comp.vot".format(daydirname, field)

    spectra = dict()
    source_ids = dict()
    if not os.path.exists(fits_filename):
        print ("Warning: File %s does not exist, skipping extraction." % \
              fits_filename)
        return spectra, source_ids

    sources = read_sources(src_filename)
    hdulist = fits.open(fits_filename)
    image = hdulist[0].data
    header = hdulist[0].header
    w = WCS(header)
    index = np.arange(header['NAXIS3'])
    velocities = w.wcs_pix2world(10,10,index[:],0,0)[2]
    for src in sources:
        pix = w.wcs_world2pix(src[0], src[1], 0, 0, 1)
        x_coord = int(round(pix[0])) - 1  # 266
        y_coord = int(round(pix[1])) - 1  # 197
        print ("Translated %.4f, %.4f to %d, %d" % (
            src[0], src[1], x_coord, y_coord))

        slice = image[0, :, y_coord, x_coord]
        l_edge, r_edge = find_edges(slice, num_edge_chan)
        print("Using data range %d - %d out of %d channels." % (
            l_edge, r_edge, len(slice)))

        # plotSpectrum(np.arange(slice.size), slice)
        spectrum_array = rec.fromarrays(
            [np.arange(slice.size)[l_edge:r_edge],
             velocities[l_edge:r_edge],
             slice[l_edge:r_edge]],
            names='plane,velocity,flux')
        c = SkyCoord(src[0], src[1], frame='icrs', unit="deg")
        spectra[c.galactic.l] = spectrum_array
        source_ids[c.galactic.l] = [src[2], src[3]]
    del image
    del header
    hdulist.close()

    return spectra, source_ids


def get_mean_continuum(spectrum, longitude, continuum_ranges):
    """
    Calculate the mean of the continuum values. Will divide the spectrum
    into continuum and line parts based on those values with a 1 sigma
    deviation from the mean of the first 5% of values.  This assumes that
    there is a leading continuum block in the spectrum.
    :param spectrum: The spectrum to be analysed, should be a numpy array of
        plane, velocity and flux values.
    :param longitude: The galactic longitude of the target object
    :param continuum_ranges: The predefined continuum blocks by longitude range
    :return: A single float which is the mean continuum flux.
    """
    int_l = int(longitude.degree)
    continuum_start_vel = -210
    continuum_end_vel = -150
    for row in continuum_ranges:
        if row[0] <= int_l <= row[1]:
            continuum_start_vel = row[2]
            continuum_end_vel = row[3]

    print(
        "Looking for velocity range %d to %d in data of %d to %d at longitude %f.3" %
        (continuum_start_vel, continuum_end_vel,
         np.min(spectrum.velocity) / 1000.0,
         np.max(spectrum.velocity) / 1000.0, longitude.degree))

    continuum_range = np.where(
        continuum_start_vel*1000 < spectrum.velocity)
    bin_start = continuum_range[0][0]
    continuum_range = np.where(
        spectrum.velocity < continuum_end_vel*1000)
    bin_end = continuum_range[0][-1]

    print("Using bins %d to %d (velocity range %d to %d) out of %d" % (
        bin_start, bin_end, continuum_start_vel, continuum_end_vel, len(spectrum.velocity)))
    continuum_sample = spectrum.flux[bin_start:bin_end]
    # print ("...gave sample of", continuum_sample)
    mean_cont = np.mean(continuum_sample)
    return mean_cont


def get_opacity(spectrum, mean):
    """
    Calculate the opacity profile for the spectrum. This simply divides the
    spectrum's flux by the mean.

    :param spectrum: The spectrum to be processed
    :param mean: The mean background flux, representing what the backlighting sources average flux.
    :return: The opacity (e^(-tau)) at each velocity step.
    """
    # print spectrum.flux
    # print spectrum.flux/mean
    return spectrum.flux/mean


def plot_spectrum(x, y, filename, title):
    """
    Output a plot of opacity vs LSR velocity to a specified file.

    :param x: The velocity data
    :param y: The opacity values for each velocity step
    :param filename: The file the plot should be written to. Should be
         an .eps or .pdf file.
    :param title: The title for the plot
    """
    fig = plt.figure()
    plt.plot(x/1000, y)

    plt.axhline(1, color='r')

    plt.xlabel(r'Velocity relative to LSR (km/s)')
    plt.ylabel(r'$e^{(-\tau)}$')
    plt.title(title)
    plt.grid(True)
    plt.savefig(filename)
    #plt.show()
    plt.close()
    return


def output_spectra(spectrum, opacity, filename, longitude, latitude):
    """
    Write the spectrum (velocity, flux and opacity) to a votable format file.

    :param spectrum: The spectrum to be output.
    :param opacity:  The opacity to be output.
    :param filename:  The filename to be created
    :param longitude: The galactic longitude of the target object
    :param latitude: The galactic latitude of the target object
    """
    table = Table([spectrum.plane, spectrum.velocity, opacity, spectrum.flux],
                  names=('plane', 'velocity', 'opacity', 'flux'),
                  meta={'name': 'Opacity', 'id' : 'opacity'})
    votable = from_table(table)
    votable.infos.append(Info('longitude', 'longitude', longitude.value))
    votable.infos.append(Info('latitude', 'latitude', latitude.value))
    writeto(votable, filename)


def produce_spectra(day_dir_name, day, field_list, continuum_ranges):
    with open(day_dir_name + '/spectra.html', 'w') as spectra_idx:
        t = Template(
            '<html>\n<head><title>D$day Spectra</title></head>\n'
            + '<body>\n<h1>Spectra previews for day $day</h1>\n<table>\n')
        spectra_idx.write(t.substitute(day=day))

        for field in field_list:
            spectra, source_ids = extract_spectra(day_dir_name, field)
            t = Template('<tr><td colspan=4><b>Field: ${field}</b></td></tr>\n' +
                         '<tr><td>Image Name</td><td> Source Longitude </td>' +
                         '<td>Peak Flux</td><td>Spectra</td></tr>\n')
            spectra_idx.write(t.substitute(field=field))

            idx = 0
            for longitude in sorted(spectra.keys()):
                spectrum = spectra.get(longitude)
                src_data = source_ids.get(longitude)
                name_prefix = field + '_src' + src_data[0]
                idx += 1
                mean = get_mean_continuum(spectrum, longitude, continuum_ranges)
                if mean < 0:
                    print(("WARNING: Skipped spectrum %s with negative " +
                          "mean: %.5f") % (name_prefix, mean))
                    continue;

                print ('Continuum mean of %s is %.5f Jy' % (name_prefix, mean))
                opacity = get_opacity(spectrum, mean)
                # print opacity
                dir_prefix = day_dir_name + "/"
                img_name = name_prefix + "_plot.png"
                plot_spectrum(spectrum.velocity, opacity, dir_prefix + img_name,
                              "Spectra for source {0} in field {1}".format(
                                  src_data[0], field))
                filename = dir_prefix + name_prefix + '_opacity.votable.xml'
                latitude = Angle(0 * u.deg) # todo: Read this from the spectra
                output_spectra(spectrum, opacity, filename, longitude, latitude)

                t = Template('<tr><td>${img}</td><td>${longitude}</td>' +
                             '<td>${peak_flux}</td><td><a href="${img}">' +
                             '<img src="${img}" width="500px"></a></td></tr>\n')
                spectra_idx.write(t.substitute(img=img_name, peak_flux=src_data[1],
                                               longitude=longitude))

        spectra_idx.write('</table></body></html>\n')


def main():
    """
    Main script for analyse_data
    :return: The exit code
    """
    # Read day parameter
    if len(sys.argv) != 2:
        print("Incorrect number of parameters.")
        print("Usage: python analyse_data.py day")
        return 1
    day = sys.argv[1]
    start = time.time()

    # Check metadata against file system
    day_dir_name = "day" + day
    if not os.path.isdir(day_dir_name):
        print ("Directory %s could not be found." % day_dir_name)
        return 1

    print ("#### Started source finding on MAGMO day %s at %s ####" % \
          (day, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))
    error_list = []

    # Read list of fields, filter for ones to be processed
    field_list = get_high_signal_fields(day_dir_name)

    # For each file, find the sources
    for field in field_list:
        error_list.extend(find_sources(day_dir_name, field))

    # For each file, extract spectra
    continuum_ranges = read_continuum_ranges()
    produce_spectra(day_dir_name, day, field_list, continuum_ranges)

    # Report
    end = time.time()
    print ('#### Processing completed at %s ####' \
          % (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))))
    print ('Searched %d images in %.02f s' % (len(field_list),
                                               end - start))
    if len(error_list) == 0:
        print ("Hooray! No errors found.")
    else:
        print ("%d errors were encountered:" % (len(error_list)))
        for err in error_list:
            print (err)
    return 0


# Run the script if it is called from the command line
if __name__ == "__main__":
    exit(main())
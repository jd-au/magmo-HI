#!/usr/bin/env python -u

# Find sources in the data and produce spectra for each suitable source.

# Author James Dempsey
# Date 28 Aug 2016

from __future__ import print_function, division

import magmo
import sgps
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
import math
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
    Search a continuum file for sources using the Aegean source finder. A
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
                     '--cores=1 --island --table=' + table_file
        magmo.run_os_cmd(aegean_cmd)
    except magmo.CommandFailedError as e:
        error_list.append(str(e))
    return error_list


def read_sources(filename):
    print ("Extracting sources from " + filename)
    sources = []

    if not os.path.exists(filename):
        print ("Warning: File %s does not exist, skipping source read." % \
               filename)
        return sources

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
            sources.append([ra, dec, id, flux, row['island']])
        else:
            print ("Ignoring source at %.4f, %.4f due to low S/N of %.4f or "
                   "flux of %.4f" % (ra, dec, sn, flux))

    return sources


def read_islands(filename):
    print ("Extracting islands from " + filename)
    islands = {}

    if not os.path.exists(filename):
        print ("Warning: File %s does not exist, skipping island read." % \
               filename)
        return {}

    isle_votable = votable.parse(filename, pedantic=False)
    results = isle_votable.get_first_table().array
    for row in results:
        islands[row['island']] = row
    return islands


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
    num_edge_chan = 10
    fits_filename = "{0}/1420/magmo-{1}_1420_sl_restor.fits".format(daydirname,
                                                                    field)
    src_filename = "{0}/{1}_src_comp.vot".format(daydirname, field)
    isle_filename = "{0}/{1}_src_isle.vot".format(daydirname, field)

    spectra = dict()
    source_ids = dict()
    if not os.path.exists(fits_filename):
        print ("Warning: File %s does not exist, skipping extraction." % \
              fits_filename)
        return spectra, source_ids

    sources = read_sources(src_filename)
    islands = read_islands(isle_filename)
    hdulist = fits.open(fits_filename)
    image = hdulist[0].data
    header = hdulist[0].header
    w = WCS(header)
    index = np.arange(header['NAXIS3'])
    beam_maj = header['BMAJ'] * 60 * 60
    beam_min = header['BMIN'] * 60 * 60
    beam_area = math.radians(header['BMAJ']) * math.radians(header['BMIN'])
    print ("Beam was %f x %f arcsec giving area of %f radians^2." % (beam_maj, beam_min, beam_area))
    velocities = w.wcs_pix2world(10,10,index[:],0,0)[2]
    for src in sources:
        pix = w.wcs_world2pix(src[0], src[1], 0, 0, 1)
        x_coord = int(round(pix[0])) - 1  # 266
        y_coord = int(round(pix[1])) - 1  # 197
        print ("Translated %.4f, %.4f to %d, %d" % (
            src[0], src[1], x_coord, y_coord))

        img_slice = image[0, :, y_coord, x_coord]
        l_edge, r_edge = find_edges(img_slice, num_edge_chan)
        print("Using data range %d - %d out of %d channels." % (
            l_edge, r_edge, len(img_slice)))

        # plotSpectrum(np.arange(slice.size), slice)
        spectrum_array = rec.fromarrays(
            [np.arange(img_slice.size)[l_edge:r_edge],
             velocities[l_edge:r_edge],
             img_slice[l_edge:r_edge]],
            names='plane,velocity,flux')
        c = SkyCoord(src[0], src[1], frame='icrs', unit="deg")
        spectra[c.galactic.l] = spectrum_array

        isle = islands.get(src[4], None)
        src_map = {'id': src[2], 'flux': src[3], 'pos': c, 'beam_area': beam_area}
        if isle:
            src_map['a'] = isle['a']
            src_map['b'] = isle['b']
            src_map['pa'] = isle['pa']
        # [ra, dec, id, flux])
        print (src_map)
        source_ids[c.galactic.l] = src_map
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
    continuum_start_vel, continuum_end_vel = magmo.lookup_continuum_range(
        continuum_ranges, int(longitude))

    print(
        "Looking for velocity range %d to %d in data of %d to %d at longitude %.3f" %
        (continuum_start_vel, continuum_end_vel,
         np.min(spectrum.velocity) / 1000.0,
         np.max(spectrum.velocity) / 1000.0, longitude))

    continuum_range = np.where(
        continuum_start_vel*1000 < spectrum.velocity)
    if len(continuum_range) ==0:
        return None, None, continuum_start_vel, continuum_end_vel

    bin_start = continuum_range[0][0]
    continuum_range = np.where(
        spectrum.velocity < continuum_end_vel*1000)
    bin_end = continuum_range[0][-1]

    print("Using bins %d to %d (velocity range %d to %d) out of %d" % (
        bin_start, bin_end, continuum_start_vel, continuum_end_vel, len(spectrum.velocity)))
    continuum_sample = spectrum.flux[bin_start:bin_end]
    # print ("...gave sample of", continuum_sample)
    mean_cont = np.mean(continuum_sample)
    sd_cont = np.std(continuum_sample)
    return mean_cont, sd_cont, continuum_start_vel, continuum_end_vel


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


def get_temp_bright(spectrum, beam_area, wavelen=0.210996048, ):
    """
    Calculate the brightness temperature (T_B) for the spectrum. This effectively converts the spectrum from units
    of Jy/beam to K.

    :param spectrum: The spectrum to be processed
    :param beam_area: The beam area in radian^2
    :return: The brightness temperature at each velocity step.
    """
    k = 1.3806503E-23  # boltzmann constant in J K^-1
    jy_to_si = 1E-26  # J s^-1 m^-2 Hz^-1

    factor = (wavelen**2 / (2*k)) * jy_to_si / (np.pi*beam_area/4)
    print (factor)
    return factor * spectrum.flux


def plot_spectrum(x, y, filename, title, con_start_vel, con_end_vel):
    """
    Output a plot of opacity vs LSR velocity to a specified file.

    :param x: The velocity data
    :param y: The opacity values for each velocity step
    :param filename: The file the plot should be written to. Should be
         an .eps or .pdf file.
    :param title: The title for the plot
    :param con_start_vel: The minimum velocity that the continuum was measured at.
    :param con_end_vel: The maximum velocity that the continuum was measured at.
    """
    fig = plt.figure()
    plt.plot(x/1000, y)

    plt.axhline(1, color='r')
    plt.axvline(con_start_vel, color='g', linestyle='dashed')
    plt.axvline(con_end_vel, color='g', linestyle='dashed')

    plt.xlabel(r'Velocity relative to LSR (km/s)')
    plt.ylabel(r'$e^{(-\tau)}$')
    plt.title(title)
    plt.grid(True)
    plt.savefig(filename)
    #plt.show()
    plt.close()
    return


def plot_emission_spectrum(velocity, em_mean, em_std, filename, title, con_start_vel, con_end_vel):
    """
    Output a plot of emission vs LSR velocity to a specified file.

    :param velocity: The velocity data
    :param em_mean: The mean temperature values for each velocity step
    :param em_std: The standard deviation in temperature values for each velocity step
    :param filename: The file the plot should be written to. Should be
         an .eps or .pdf file.
    :param title: The title for the plot
    :param con_start_vel: The minimum velocity that the continuum was measured at.
    :param con_end_vel: The maximum velocity that the continuum was measured at.
    """

    if len(em_mean) == 0:
        if os.path.exists(filename):
            os.remove(filename)
        return

    fig = plt.figure()
    plt.plot(velocity/1000, em_mean)

    em_max = em_mean + em_std
    em_min = em_mean - em_std
    plt.fill_between(velocity/1000, em_min, em_max, facecolor='lightgray', color='lightgray')

    plt.axvline(con_start_vel, color='g', linestyle='dashed')
    plt.axvline(con_end_vel, color='g', linestyle='dashed')

    plt.xlabel(r'Velocity relative to LSR (km/s)')
    plt.ylabel(r'$T_B$ (K)')
    plt.title(title)
    plt.grid(True)
    plt.savefig(filename)
    #plt.show()
    plt.close()
    return


def output_spectra(spectrum, opacity, filename, longitude, latitude, em_mean, em_std, temp_bright, beam_area):
    """
    Write the spectrum (velocity, flux and opacity) to a votable format file.

    :param spectrum: The spectrum to be output.
    :param opacity:  The opacity to be output.
    :param filename:  The filename to be created
    :param longitude: The galactic longitude of the target object
    :param latitude: The galactic latitude of the target object
    """
    table = Table(meta={'name': filename, 'id': 'opacity'})
    table.add_column(Column(name='plane', data=spectrum.plane))
    table.add_column(Column(name='velocity', data=spectrum.velocity, unit='m/s'))
    table.add_column(Column(name='opacity', data=opacity))
    table.add_column(Column(name='flux', data=spectrum.flux, unit='Jy', description='Flux per beam'))
    table.add_column(Column(name='temp_brightness', data=temp_bright, unit='K'))
    if len(em_mean) > 0:
        # The emission may not be available, so don't include it if not
        table.add_column(Column(name='em_mean', data=em_mean, unit='K'))
        table.add_column(Column(name='em_std', data=em_std, unit='K'))

    votable = from_table(table)
    votable.infos.append(Info('longitude', 'longitude', longitude.value))
    votable.infos.append(Info('latitude', 'latitude', latitude.value))
    votable.infos.append(Info('beam_area', 'beam_area', beam_area))
    writeto(votable, filename)


def output_emission_spectra(filename, longitude, latitude, velocity, em_mean,
                            em_std, ems):
    """
    Write the emission spectrum (velocity, flux and opacity) to a votable format
    file.

    :param filename:  The filename to be created
    :param longitude: The galactic longitude of the target object
    :param latitude: The galactic latitude of the target object
    :param velocity:
    :param em_mean:
    :param em_std:
    :param ems:
    """
    table = Table(meta={'name': filename, 'id': 'emission'})
    table.add_column(Column(name='velocity', data=velocity, unit='m/s'))
    table.add_column(Column(name='em_mean', data=em_mean, unit='K'))
    table.add_column(Column(name='em_std', data=em_std, unit='K'))
    for i in range(len(ems)):
        table.add_column(Column(name='em_'+str(i), data=ems[i].flux, unit='K'))

    votable = from_table(table)
    votable.infos.append(Info('longitude', 'longitude', longitude.value))
    votable.infos.append(Info('latitude', 'latitude', latitude.value))
    writeto(votable, filename)


def calc_offset_points(longitude, latitude, beam_size, num_points=6):
    spacing = 2.0 * math.pi / float(num_points)
    points = []
    for i in range(0, num_points):
        angle = spacing * i
        l = longitude + math.sin(angle)*beam_size
        b = latitude + math.cos(angle)*beam_size
        coord = SkyCoord(l, b, frame='galactic', unit="deg")
        print ("Point at angle %f is %s" % (math.degrees(angle), str(coord)))
        points.append(coord)

    return points


def get_emission_spectra(centre, velocities, file_list, filename_prefix):
    """
    Extract SGPS emission spectra around a central point.

    :param centre: A SkyCoord containing the location of the central point
    :param velocities: The velocities list sothat the emission data can be matched.
    :param file_list: A list of dictionaries describing the SGPS files.
    :return: An array fo the mean and standard deviation of emission at each velocity.
    """

    #file_list = sgps.get_hi_file_list()
    filename = filename_prefix + '_emission.votable.xml'
    coords = calc_offset_points(centre.galactic.l.value,
                                centre.galactic.b.value, 0.03611)
    ems = sgps.extract_spectra(coords, file_list)
    if ems:
        all_em = np.array([ems[i].flux for i in range(len(ems))])
        em_std = np.std(all_em, axis=0)
        em_mean = np.mean(all_em, axis=0)
        em_std_interp = np.interp(velocities, ems[0].velocity, em_std)
        em_mean_interp = np.interp(velocities, ems[0].velocity, em_mean)

        output_emission_spectra(filename,
                                centre.galactic.l, centre.galactic.b,
                                ems[0].velocity, em_mean, em_std, ems)

        return em_mean_interp, em_std_interp

    print("WARNING: Unable to find emission data for " + str(centre))
    if os.path.exists(filename):
        os.remove(filename)
    return [], []


def produce_spectra(day_dir_name, day, field_list, continuum_ranges):
    file_list = sgps.get_hi_file_list()
    with open(day_dir_name + '/spectra.html', 'w') as spectra_idx:
        t = Template(
            '<html>\n<head><title>D$day Spectra</title></head>\n'
            + '<body>\n<h1>Spectra previews for day $day</h1>\n<table>\n')
        spectra_idx.write(t.substitute(day=day))

        neg_mean = 0
        no_mean = 0
        all_cont_sd = []
        all_opacity = []
        for field in field_list:
            spectra, source_ids = extract_spectra(day_dir_name, field)
            t = Template('<tr><td colspan=4><b>Field: ${field}</b></td></tr>\n' +
                         '<tr><td>Image Name</td><td>Details</td>' +
                         '<td>Absorption</td><td>Emission</td></tr>\n')
            spectra_idx.write(t.substitute(field=field))

            idx = 0
            for longitude in sorted(spectra.keys()):
                spectrum = spectra.get(longitude)
                src_data = source_ids.get(longitude)
                name_prefix = field + '_src' + src_data['id']
                idx += 1
                mean, cont_sd, min_con_vel, max_con_vel = get_mean_continuum(
                    spectrum,
                    longitude.degree,
                    continuum_ranges)
                if mean is None:
                    print("WARNING: Skipped spectrum %s with no continuum data" % (name_prefix, mean))
                    no_mean += 1
                    continue;

                if mean < 0:
                    print(("WARNING: Skipped spectrum %s with negative " +
                          "mean: %.5f") % (name_prefix, mean))
                    neg_mean += 1
                    continue;

                print('Continuum mean of %s is %.5f Jy, sd %.5f' % (
                    name_prefix, mean, cont_sd))
                all_cont_sd.append(cont_sd)
                opacity = get_opacity(spectrum, mean)
                temp_bright = get_temp_bright(spectrum, src_data['beam_area'])
                # print opacity
                dir_prefix = day_dir_name + "/"
                img_name = name_prefix + "_plot.png"
                plot_spectrum(spectrum.velocity, opacity, dir_prefix + img_name,
                              "Spectra for source {0} in field {1}".format(
                                  src_data['id'], field), min_con_vel, max_con_vel)
                filename = dir_prefix + name_prefix + '_opacity.votable.xml'
                latitude = src_data[2].galactic.b

                em_mean, em_std = get_emission_spectra(src_data['pos'],
                                                       spectrum.velocity,
                                                       file_list, dir_prefix + name_prefix)
                em_img_name = name_prefix + "_emission.png"
                plot_emission_spectrum(spectrum.velocity, em_mean, em_std,
                                       dir_prefix + name_prefix + "_emission.png",
                                       "Emission around source {0} in field {1}".format(
                                           src_data['id'], field), min_con_vel,
                                       max_con_vel)
                output_spectra(spectrum, opacity, filename, longitude, latitude,
                               em_mean, em_std, temp_bright, src_data['beam_area'])
                all_opacity.append(opacity)

                t = Template('<tr><td>${img}</td><td>l:&nbsp;${longitude}<br/>' +
                             'Peak:&nbsp;${peak_flux}&nbsp;Jy<br/>Mean:&nbsp;${mean}&nbsp;Jy<br/>'
                             'Cont&nbsp;SD:&nbsp;${cont_sd}</td><td><a href="${img}">' +
                             '<img src="${img}" width="500px"></a></td><td><a href="${em_img}">' +
                             '<img src="${em_img}" width="500px"></a></td></tr>\n')
                spectra_idx.write(t.substitute(img=img_name, em_img=em_img_name, peak_flux=src_data['flux'],
                                               longitude=longitude, mean=mean, cont_sd=cont_sd))


        spectra_idx.write('</table></body></html>\n')

        if no_mean > 0:
            print("Skipped %d spectra with no continuum data." % no_mean)

        print("Skipped %d spectra with negative mean continuum." % neg_mean)
        print("Produced %d spectra with continuum sd of %.5f." % (
            len(all_cont_sd), np.mean(all_cont_sd)))
        return all_opacity


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
    continuum_ranges = magmo.get_continuum_ranges()
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
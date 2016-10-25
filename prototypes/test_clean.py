#!/usr/bin/env python -u

# Try out various clean strategies to find optimal approach for the magmo data

# Author James Dempsey
# Date 23 Oct 2016

from __future__ import print_function, division

import magmo
import analyse_data
import process_data

import numpy as np
import os
import time
import csv


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
                if str(row[4]) == "Y":
                    print (row)
                    field_list.append(row[0])

    return field_list


def find_freq_file(name_prefix, freq_list):
    prefixes = [name_prefix]
    if name_prefix.startswith("0"):
        subname = name_prefix.lstrip("0")
        if subname.startswith("."):
            subname = "0" + subname
        prefixes.append(subname)
    for prefix in prefixes:
        for src_freq in freq_list:
            bp_file = prefix + '.' + src_freq
            print ('Looking for', bp_file)
            if os.path.exists(bp_file):
                return bp_file, src_freq

    return None, None


def clean_spectral_data(day_dir_name, band, fields):

    error_list = []
    os.chdir(day_dir_name)
    # We only work with the HI spectral line data
    freq = band['main']
    magmo.ensure_dir_exists(freq)
    print ("Cleaning %s MHz spectral cubes for %d fields" % (freq, len(fields)))

    for src_name in fields:
        try:
            src_file, freq_suffix = find_freq_file(src_name, band['freqs'])
            if src_file is None:
                print ("No spectral line file found for field %s and band %s." % (src_name, freq))
                exit(1)
            name_prefix = freq + '/magmo-' + src_name + '_' + freq + '_'
            dirty_file = name_prefix + 'sl_dirty'
            clean_file = name_prefix + 'sl_clean'
            beam_file = name_prefix + 'sl_beam'

            rms_1420, max_1420 = process_data.get_signal_noise_ratio(
                dirty_file)
            cutoff = 3 * rms_1420

            magmo.run_os_cmd('rm -r ' + clean_file, failOnErr=False)
            #cmd = 'clean niters=2000 speed=+1 map=' + dirty_file + ' beam=' \
            #    + beam_file + ' out=' + clean_file
            cmd = 'clean niters=250 speed=+1 options=positive cutoff=' + str(cutoff) + \
                  ' map=' + dirty_file + ' beam=' + beam_file + ' out=' + \
                  clean_file
            magmo.run_os_cmd(cmd)
        except magmo.CommandFailedError as e:
            error_list.append(str(e))

    os.chdir('..')
    return error_list


def produce_cube(day_dir_name, band, fields):

    error_list = []
    os.chdir(day_dir_name)
    # We only work with the HI spectral line data
    freq = band['main']
    magmo.ensure_dir_exists(freq)
    print ("Producing %s MHz spectral cubes for %d fields" % (freq, len(fields)))

    for src_name in fields:
        try:
            src_file, freq_suffix = find_freq_file(src_name, band['freqs'])
            if src_file is None:
                print ("No spectral line file found for field %s and band %s." % (src_name, freq))
                exit(1)
            name_prefix = freq + '/magmo-' + src_name + '_' + freq + '_'
            dirty_file = name_prefix + 'sl_dirty'
            clean_file = name_prefix + 'sl_clean'
            beam_file = name_prefix + 'sl_beam'
            restored_file = name_prefix + 'sl_restor'
            fits_file = restored_file + '.fits'

            magmo.run_os_cmd('rm -r ' + restored_file + ' ' + fits_file, failOnErr=False)
            cmd = 'restor model=' + clean_file + ' beam=' \
                + beam_file + ' map=' + dirty_file + ' out=' + restored_file
            magmo.run_os_cmd(cmd)
            cmd = 'fits op=xyout in=' + restored_file + ' out=' + fits_file
            magmo.run_os_cmd(cmd)
        except magmo.CommandFailedError as e:
            error_list.append(str(e))

    os.chdir('..')
    return error_list


def report_cube_stats(day_dir_name, fields):

    rms_list = []
    max_list = []
    for src_name in fields:
        restored_img = day_dir_name+"/1420/magmo-" + src_name + "_1420_sl_restor"
        if os.path.exists(restored_img):
            rms_1420, max_1420 = process_data.get_signal_noise_ratio(
                restored_img)
            rms_list.append(rms_1420)
            max_list.append(max_1420)

    rms_arr = np.array(rms_list)
    max_arr = np.array(max_list)
    print(
        "Produced %d cubes\n Mean rms %.4f\n   sd rms %.4f\n mean max %.4f\n   sd max %.4f" % (
            len(rms_list), np.mean(rms_arr), np.max(rms_arr), np.mean(max_arr),
            np.max(max_arr)))

    return


def report_spectra_stats(all_spectra):

    # Compare min opacity and rms_opacity for each spectra,
    # likely using mean, median, min, and max stats
    ok_spectra = 0
    good_spectra = 0
    all_min = []
    all_max = []

    for spectrum in all_spectra:
        min_opacity = np.min(spectrum)
        max_opacity = np.max(spectrum)
        all_min.append(min_opacity)
        all_max.append(max_opacity)

        if min_opacity > -4 and max_opacity < 4:
            ok_spectra += 1
        if min_opacity > -2 and max_opacity < 2:
            good_spectra += 1

    print("Found %d ok and %d good spectra, mean min: %.4f, mean max: %.4f" % (
        ok_spectra, good_spectra, np.mean(all_min), np.mean(all_max)))
    return


def main():
    day = "21"
    start = time.time()

    print ("#### Started processing MAGMO day %s at %s ####" % \
          (day, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))
    error_list = []

    # Check metadata against file system
    day_dir_name = "day" + day
    if not os.path.isdir(day_dir_name):
        print ("Directory %s could not be found." % day_dir_name)
        exit(1)

    # set up map of parent/child map of frequencies
    line_band = {'main': '1420', 'freqs': ['1420', '1421', '1420.5'], 'line': True}
    cont_band = {'main': '1757', 'freqs': ['1757'], 'line': False}
    band_list = [line_band, cont_band]
    for band in band_list:
        freq = band['main']
        freq_dir = day_dir_name + "/" + freq
        magmo.ensure_dir_exists(freq_dir)

    fields = get_high_signal_fields(day_dir_name)
    print ("Found %d fields. First field: %s" % (len(fields), fields[0]))

    error_list.extend(clean_spectral_data(day_dir_name, line_band, fields))

    error_list.extend(produce_cube(day_dir_name, line_band, fields))

    # For each file, extract spectra
    magmo.run_os_cmd('rm -r ' + day_dir_name+"/*votable.xml", failOnErr=False)
    magmo.run_os_cmd('rm -r ' + day_dir_name+"/*plot.png", failOnErr=False)
    continuum_ranges = magmo.get_continuum_ranges()
    all_spectra = analyse_data.produce_spectra(day_dir_name, day, fields,
                                               continuum_ranges)

    # Report
    report_cube_stats(day_dir_name, fields)
    report_spectra_stats(all_spectra)

    end = time.time()
    print ('#### Processing completed at %s ####' \
          % (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))))
    print('Cleaned %d images in %.02f s' % (len(fields),
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
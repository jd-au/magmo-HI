# Flag, calibrate, produce 2 GHz continuum image, produce HI image cube
#
# This program extracts the HI spectra (1.42 GHz) and continuum (2 GHz) data
# from the day's RPFITS files into miriad format. It then backs up the metadata
# to allow any processing to be simply rolled back.

# TODO: Add dynamic flagging support - use an issues file
# TODO: Generate a calibration report with diagnostic plots
# TODO: Generate an image index page with previews of each image

# Author James Dempsey
# Date 30 Jul 2016

import csv
import sys
import os
import glob
import magmo
import time as timer


# functions to be made common
def get_day_obs_data(day):
    """
    Read the magmo-obs.csv file and find the rows for the requested day.
    :param day: The day to be found
    :return: The day's rows, or None if the day is not defined
    """

    sources = []
    with open('magmo-obs.csv', 'rb') as magmo_obs:
        reader = csv.reader(magmo_obs)
        for row in reader:
            if row[0] == day:
                src = dict()
                src['source'] = row[4]
                src['phase_cal'] = row[10]
                src['gal_l'] = row[2]
                sources.append(src)
    return sources


# functions specific to this step
def flag_data(dirname, day):
    """
    Flag out any bad data in the visibilities. This is a combinationof global flagging
    and using previously compiled flaggin lists for the day's data.

    TODO: Add read of file and flagging based on the file

    :param dirname: The name of the day directory
    :param day: The day number being processed.
    :return: None
    """
    # Need to add dynamic flagging e.g. day7: antenna 6,pol x and antenna 4
    uvDirs = glob.glob(dirname + '/*.[0-9][0-9][0-9][0-9]')
    for filename in uvDirs:
        uvflag_cmd = "uvflag flagval=f options=brief vis='"+filename+"' select='amplitude(500)'"
        magmo.run_os_cmd(uvflag_cmd)
        # todo: make these configurable
        uvflag_cmd = "uvflag flagval=f options=brief vis='" + filename + "' select='ant(4)'"
        magmo.run_os_cmd(uvflag_cmd)
        uvflag_cmd = "uvflag flagval=f options=brief vis='" + filename + "' select='ant(6),pol(xx)'"
        magmo.run_os_cmd(uvflag_cmd)

    return []


def find_bandpasscal(dirname):
    """
    Identify the bandpass calibrator that was used for the day.
    :param dirname: The name of the day directory
    :return: The source name of the calibrator, or None if no know calibrators were present.
    """
    potential_cals = ['1934-638']
    for cal in potential_cals:
        path = dirname + '/' + cal + '.*'
        uvdirs = glob.glob(path)
        for filename in uvdirs:
            #print filename
            if os.path.isdir(filename):
                return cal
    print "Unable to find any of the %s cals in %s " % (potential_cals, dirname)
    return None


def find_freq_file(name_prefix, freq_list):
    for src_freq in freq_list:
        bp_file = name_prefix + '.' + src_freq
        print 'Looking for', bp_file
        if os.path.exists(bp_file):
            return bp_file, src_freq

    return None, None


def calibrate(dirname, bandpass_cal, sources, band_list):
    """
    Prepare the bandpass, flux and phase calibration data.

    TODO: Consider outputing a calibration report with plots.

    :param dirname: The name of the day directory
    :param bandpass_cal: The source name of the bandpass & flux calibrator
    :param sources: The list of sources, which includes the pahse calibrator details
    :param band_list: The list of bands being processed, each band is a map
    :return: error list
    """
    error_list = []
    phase_cals = set()
    for src in sources:
        phase_cals.add(src["phase_cal"])
    print phase_cals

    for band in band_list:
        freq = band['main']
        print "###### Calibrating freq " + freq + " ######"
        bp_file, freq_suffix = find_freq_file(dirname + '/' + bandpass_cal, band['freqs'])
        if bp_file is None:
            print "No bandpass file found for band %s." % (freq)
            exit(1)

        magmo.run_os_cmd("mfcal vis="+bp_file+" options=interpolate")
        magmo.run_os_cmd("gpcal vis="+bp_file+" options=xyvary")

        for cal in phase_cals:
            for src_freq in band['freqs']:
                cal_file = dirname + '/' + cal + '.' + freq
                if os.path.exists(cal_file):
                    try:
                        print "##--## Processing phase cal " + cal_file + " ##--##"
                        magmo.run_os_cmd('gpcopy vis='+bp_file+' out='+cal_file)
                        magmo.run_os_cmd('gpcal vis='+cal_file+' options=xyvary,qusolv')
                        magmo.run_os_cmd('gpboot vis='+cal_file+' cal='+bp_file)
                        magmo.run_os_cmd('mfboot vis=' + cal_file + ',' + bp_file +
                                         ' "select=source(' + bandpass_cal + ')"')
                    except magmo.CommandFailedError as e:
                        error_list.append(str(e))

    return error_list


def build_images(day_dir_name, sources, band):
    """
    Generate continuum images of each field.

    :param day_dir_name: The name of the day directory
    :param sources: The list of sources, which includes the pahse calibrator details
    :param band: The map describing the continuum band to be processed
    :return: None
    """

    error_list = []
    os.chdir(day_dir_name)
    # We only work with the continuum data
    freq = band['main']
    magmo.ensure_dir_exists(freq)
    print "Producing %s MHz continuum images for %d sources" % (freq, len(sources))

    for src in sources:
        src_name = src['source']
        try:
            src_file, freq_suffix = find_freq_file(src_name, band['freqs'])
            if src_file is None:
                print "No continuum file found for source %s and band %s." % (src_name, freq)
                exit(1)

            phase_cal_file = src['phase_cal'] + "." + freq_suffix
            magmo.run_os_cmd('gpcopy vis=' + phase_cal_file + ' out=' + src_file)

            name_prefix = freq + '/magmo-' + src_name + '_' + freq + '_'
            dirty_file = name_prefix + 'dirty'
            clean_file = name_prefix + 'clean'
            beam_file = name_prefix + 'beam'
            restored_file = name_prefix + 'restor'
            fits_file = restored_file + '.fits'

            cmd = 'invert robust=0.5 options=systemp,mfs,double stokes=ii vis=' + src_file \
                + ' map=' + dirty_file + ' beam=' + beam_file
            magmo.run_os_cmd(cmd)
            cmd = 'clean niters=2000 speed=+1 map=' + dirty_file + ' beam=' \
                + beam_file + ' out=' + clean_file
            magmo.run_os_cmd(cmd)
            cmd = 'restor options=mfs model=' + clean_file + ' beam=' \
                + beam_file + ' map=' + dirty_file + ' out=' + restored_file
            magmo.run_os_cmd(cmd)
            cmd = 'fits op=xyout in=' + restored_file + ' out=' + fits_file
            magmo.run_os_cmd(cmd)
        except magmo.CommandFailedError as e:
            error_list.append(str(e))

    os.chdir('..')
    return error_list


def build_cubes(day_dir_name, sources, band):
    """
    Generate HI spectral cubes of each field.

    :param day_dir_name: The name of the day directory
    :param sources: The list of sources, which includes the pahse calibrator details
    :param freq_list: The list of short identifiers of the frequencies
    :return: None
    """

    error_list = []
    os.chdir(day_dir_name)
    # We only work with the HI spectral line data
    freq = band['main']
    magmo.ensure_dir_exists(freq)
    print "Producing %s MHz spectral cubes for %d sources" % (freq, len(sources))

    for src in sources:
        src_name = src['source']
        try:
            src_file, freq_suffix = find_freq_file(src_name, band['freqs'])
            if src_file is None:
                print "No spectral line file found for source %s and band %s." % (src_name, freq)
                exit(1)

            phase_cal_file = src['phase_cal'] + "." + freq_suffix
            magmo.run_os_cmd('gpcopy vis=' + phase_cal_file + ' out=' + src_file)

            name_prefix = freq + '/magmo-' + src_name + '_' + freq + '_'
            ave_file = name_prefix + 'ave'
            dirty_file = name_prefix + 'sl_dirty'
            clean_file = name_prefix + 'sl_clean'
            beam_file = name_prefix + 'sl_beam'
            restored_file = name_prefix + 'sl_restor'
            fits_file = restored_file + '.fits'
            #line = 'felocity,1053,-250.0,0.4,0.4'
            line = 'felocity,793,-225.0,0.4,0.4'

            cmd = 'uvaver line=' + line + ' vis=' + src_file + ' out=' + ave_file
            magmo.run_os_cmd(cmd)
            cmd = 'invert robust=0.5 cell=5 options=systemp,nopol,mosaic,double stokes=i '\
                + ' slop=0.5 line='+ line + ' vis=' + ave_file \
                + ' map=' + dirty_file + ' beam=' + beam_file
            magmo.run_os_cmd(cmd)
            cmd = 'clean niters=2000 speed=+1 map=' + dirty_file + ' beam=' \
                + beam_file + ' out=' + clean_file
            magmo.run_os_cmd(cmd)
            cmd = 'restor model=' + clean_file + ' beam=' \
                + beam_file + ' map=' + dirty_file + ' out=' + restored_file
            magmo.run_os_cmd(cmd)
            cmd = 'fits op=xyout in=' + restored_file + ' out=' + fits_file
            magmo.run_os_cmd(cmd)
        except magmo.CommandFailedError as e:
            error_list.append(str(e))

    os.chdir('..')
    return error_list


# ### Script starts here ###

# Read day parameter
if len(sys.argv) != 2:
    print("Incorrect number of parameters.")
    print("Usage: python process-data.py day")
    exit(1)
day = sys.argv[1]
start = timer.clock()

# metadata needed: flux/bandpass cal, phase cals for each source, extra flagging
# Read metadata for the day (file pattern fragments etc)
sources = get_day_obs_data(day)
if sources is None or len(sources) == 0:
    print "Day %s is not defined." % (day)
    exit(1)
print "Found %d sources. First source:\n%s" % (len(sources), sources[0])
error_list = []

# Check metadata against file system
dayDirName = "day" + day
if not os.path.isdir(dayDirName):
    print "Directory %s could not be found." % dayDirName
    exit(1)

# set up map of parent/child map of frequencies
line_band = {'main': '1420', 'freqs': ['1420', '1421']}
cont_band = {'main': '1720', 'freqs': ['1720', '1721']}
band_list = [line_band, cont_band]
for band in band_list:
    freq = band['main']
    freqDir = dayDirName + "/" + freq
    magmo.ensure_dir_exists(freqDir)

# Flag
error_list.extend(flag_data(dayDirName, day))

# Calibrate
bandpasscal = find_bandpasscal(dayDirName)
print "Bandpass cal:", bandpasscal
error_list.extend(calibrate(dayDirName, bandpasscal, sources, band_list))
print error_list

# Produce 2 GHz continuum image
error_list.extend(build_images(dayDirName, sources, cont_band))

# Produce HI image cube
error_list.extend(build_cubes(dayDirName, sources, line_band))

# Report
end = timer.clock()
print '#### Processing Completed ####'
print 'Processed %d sources in %.03f s at %s' % (len(sources), end - start, timer.time())
if len(error_list) == 0:
    print "Hooray! No errors found."
else:
    print "%d errors were encountered:" % (len(error_list))
    for err in error_list:
        print err

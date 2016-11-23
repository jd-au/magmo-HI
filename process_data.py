# Flag, calibrate, produce 2 GHz continuum image, produce HI image cube
#
# This program extracts the HI spectra (1.42 GHz) and continuum (2 GHz) data
# from the day's RPFITS files into miriad format. It then backs up the metadata
# to allow any processing to be simply rolled back.

# Author James Dempsey
# Date 30 Jul 2016

import aplpy
import csv
import glob
import magmo
import math
import os
import re
import subprocess
import sys
import time
import numpy as np
from collections import OrderedDict

from string import Template
from astropy.io import ascii

sn_min = 1.3
num_chan = 627


def get_day_flags(day):
    """
    Read the magmo-flagging.csv file and find the flagging needed for the requested day.
    :param day: The day to be found
    :return: The uvflag select clauses needed
    """

    flags = []
    with open('magmo-flagging.csv', 'rb') as magmo_obs:
        reader = csv.reader(magmo_obs)
        for row in reader:
            if row[0] == day:
                select = ''
                if len(row[1]) > 0:
                    select = 'time('+row[1]+','+row[2]+')'
                if len(row[3]) > 0:
                    if len(select) > 0:
                        select += ','
                    select += 'ant('+row[3]+')'
                    if len(row[4]) > 0:
                        select += '('+row[4]+')'
                if len(row[5]) > 0:
                    if len(select) > 0:
                        select += ','
                    select += 'pol('+row[5]+')'
                flag_details = {'uv_file':row[6], 'select':select, 'line':row[7]}
                flags.append(flag_details)
    return flags


def find_zero_hi_plane(dirname, bandpass_cal):
    """
    Determine the plane containing the 0 km/s HI line for the bandpass
    calibrator. This can then be used as the basis for flaging out the local
    HI absorbtion in the bandpass cal.

    :param dirname: The name of the day directory
    :param bandpass_cal: The source name of the bandpass calibrator
    :return: The integer plane (aka line)  number of the zero km/s plane.
    """
    cal_file = "{0}/{1}.1420".format(dirname, bandpass_cal)
    result = subprocess.check_output(
        ["uvlist", "vis=" + cal_file, "options=spectral"])
    opt_vel_next = False
    vel_start = 0.0
    vel_end = 0.0
    vel_step = 0.0
    for line in result.splitlines():
        print line
        if opt_vel_next:
            stats = line.strip().split(':')
            if len(stats) > 1:
                if stats[0].strip().startswith('Start velocity'):
                    vel_start = float(stats[1])
                if stats[0].strip().startswith('End velocity'):
                    vel_end = float(stats[1])
                if stats[0].strip().startswith('Velocity'):
                    vel_step = float(stats[1])
        elif re.match('^Optical Velocities:', line):
            opt_vel_next = True

    if vel_step != 0.0:
        # print vel_start, vel_end, vel_step
        plane = int(vel_start / vel_step * -1.0)
        return plane
    else:
        print "Unable to find Optical Velocity data in:"
        print result
        return None


def flag_data(dirname, day, bandpass_cal, sources):
    """
    Flag out any bad data in the visibilities. This is a combination of global
    flagging and using previously compiled flagging lists for the day's data.

    :param dirname: The name of the day directory
    :param day: The day number being processed.
    :param bandpass_cal: The source name of the bandpass calibrator
    :param sources: The list of sources, which includes the phase calibrator
           details
    :return: None
    """

    calibrators = set()
    # calibrators.add(bandpass_cal) # tvclip seems to make a mess of the bamdpass calibrator
    for src in sources:
        calibrators.add(src["phase_cal"])
    dynamic_flags = get_day_flags(day)
    uvDirs = glob.glob(dirname + '/*.[0-9][0-9][0-9][0-9]')
    for filename in uvDirs:
        uvflag_cmd = "uvflag flagval=f options=brief vis='"+filename+"' select='amplitude(500)'"
        magmo.run_os_cmd(uvflag_cmd)

        for flag in dynamic_flags:
            flag_file = flag['uv_file']
            if flag_file is None or len(flag_file) == 0 or filename.endswith(flag_file):
                uvflag_cmd = "uvflag flagval=f options=brief vis='" + filename + "' "
                if len(flag['select']) > 0:
                    uvflag_cmd += " select='" + flag['select'] + "'"
                if len(flag['line']) > 0:
                    uvflag_cmd += " line='" + flag['line'] + "'"
            magmo.run_os_cmd(uvflag_cmd)

        # Discard edge channels
        if filename.endswith('1757'):
            uvflag_cmd = "uvflag flagval=f options=brief vis='" + filename + "' edge=20,550"
            magmo.run_os_cmd(uvflag_cmd)
        else:
            uvflag_cmd = "uvflag flagval=f options=brief vis='" + filename + "' edge=200,200"
            magmo.run_os_cmd(uvflag_cmd)

    # Flag the 0km/s region in the bandpass calibrator
    zero_km_plane = find_zero_hi_plane(dirname, bandpass_cal)
    filename = "{0}/{1}.1420".format(dirname, bandpass_cal)
    line = "channel,110,{0},1,1".format(zero_km_plane - 35)
    uvflag_cmd = "uvflag flagval=f options=brief vis='" + filename + \
                 "' line=" + line
    magmo.run_os_cmd(uvflag_cmd)

    # Clip the data too far away for the mean in the calibrators
    for cal_src in calibrators:
        pattern = dirname + '/' + cal_src + '.[0-9][0-9][0-9][0-9]'
        print pattern
        uvDirs = glob.glob(pattern)
        for filename in uvDirs:
            tvclip_cmd = "tvclip options=notv commands=clip,diff,clip clip=3 vis='" + filename + "'"
            magmo.run_os_cmd(tvclip_cmd)

    return []


def find_bandpasscal(dirname):
    """
    Identify the bandpass calibrator that was used for the day.
    :param dirname: The name of the day directory
    :return: The source name of the calibrator, or None if no know calibrators were present.
    """
    potential_cals = ['1934-638', '0823-500']
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


def calibrate(dirname, bandpass_cal, sources, band_list, day):
    """
    Prepare the bandpass, flux and phase calibration data.

    :param dirname: The name of the day directory
    :param bandpass_cal: The source name of the bandpass & flux calibrator
    :param sources: The list of sources, which includes the phase calibrator details
    :param band_list: The list of bands being processed, each band is a map
    :param day: The day number being processed
    :return: error list
    """
    error_list = []
    phase_cals = set()
    for src in sources:
        phase_cals.add(src["phase_cal"])
    print phase_cals

    # Setup te index file
    cal_idx = open(dirname + '/calibration.html', 'w')
    t = Template('<html>\n<head><title>Day $day calibration</title></head>\n'
                 + '<body>\n<h1>Calibration results for day $day</h1>\n<table>')
    cal_idx.write(t.substitute(day=day))

    for band in band_list:
        freq = band['main']
        print "###### Calibrating freq " + freq + " ######"
        bp_file, freq_suffix = find_freq_file(dirname + '/' + bandpass_cal, band['freqs'])
        if bp_file is None:
            print "No bandpass file found for band %s." % (freq)
            exit(1)

        cal_interval = 0.25
        magmo.run_os_cmd(
            "mfcal vis=" + bp_file + " options=interpolate interval=" + str(
                cal_interval))
        magmo.run_os_cmd(
            "gpcal vis=" + bp_file + " options=xyvary interval=" + str(
                cal_interval))

        # Plot bandpass of bandpass cal
        magmo.run_os_cmd(
            "gpplt options=bandpass vis=" + bp_file + " device=" + bp_file +
            "-bandpass.png/png")
        magmo.run_os_cmd(
            "uvplt options=nobase vis=" + bp_file + " device=" + bp_file +
            "-amp.png/png")
        magmo.run_os_cmd(
            "uvplt options=equal,nobase stokes=i,q,u,v axis=real,imag vis=" +
            bp_file + " device=" + bp_file + "-scatter.png/png ")

        t = Template('<tr><td colspan="2"><br>Bandpass of $bp_cal at $src_freq MHz</td></tr>\n' +
                     '<tr>\n<td>' +
                     '<a href="${bp_file}-bandpass.png"><img src="${bp_file}-bandpass.png" width="400px"></a>' +
                     '</td><td>'
                     '<a href="${bp_file}-amp.png"><img src="${bp_file}-amp.png" width="400px"></a>' +
                     '</td><td>'
                     '<a href="${bp_file}-scatter.png"><img src="${bp_file}-scatter.png" width="400px"></a>' +
                     '</td>\n</tr>')
        cal_idx.write(t.substitute(bp_cal=bandpass_cal, src_freq=freq,
                                   bp_file=bp_file[len(dirname) + 1:]))

        phase_cal_interval = 5  # min
        for cal in phase_cals:
            for src_freq in band['freqs']:
                cal_file = dirname + '/' + cal + '.' + src_freq
                if os.path.exists(cal_file):
                    try:
                        print "##--## Processing phase cal " + cal_file + " ##--##"
                        magmo.run_os_cmd('gpcopy vis='+bp_file+' out='+cal_file)
                        magmo.run_os_cmd(
                            'gpcal vis=' + cal_file + ' options=xyvary,qusolv' +
                            ' interval=' + str(phase_cal_interval))
                        magmo.run_os_cmd('gpboot vis='+cal_file+' cal='+bp_file)
                        magmo.run_os_cmd('mfboot vis=' + cal_file + ',' + bp_file +
                                         ' "select=source(' + bandpass_cal + ')"')
                    except magmo.CommandFailedError as e:
                        error_list.append(str(e))

                    prefix = dirname + '/' + freq + '/cal-' + cal + "-" + src_freq
                    cmd = 'uvplt stokes=i,q,u,v axis=real,imag vis=' + cal_file \
                          + ' device=' + prefix + '.png/png options=equal'
                    magmo.run_os_cmd(cmd, False)
                    cmd = 'gpplt yaxis=phase options=xygains vis=' + cal_file \
                          + ' device=' + prefix + '-gain.png/png'
                    magmo.run_os_cmd(cmd, False)
                    cmd = 'uvplt axis=time,amplitude vis=' + cal_file \
                          + ' device=' + prefix + '-amp.png/png'
                    magmo.run_os_cmd(cmd, False)
                    t = Template('<tr><td colspan="3"><br>Source $cal at $src_freq MHz</td></tr>\n'
                                 + '<tr>\n'
                                 + '<td><a href="${prefix}.png"><img src="${prefix}.png" width="400px"></a></td></td>'
                                 + '<td><a href="${prefix}-amp.png"><img src="${prefix}-amp.png" width="400px"></a></td>'
                                 + '<td><a href="${prefix}-gain.png"><img src="${prefix}-gain.png" width="400px"></a></td>'
                                 + '</tr>'
                                + '<tr>\n'
                                + '<td><a href="${prefix}.png_2"><img src="${prefix}.png_2" width="400px"></a></td></td>'
                                + '<td><a href="${prefix}-amp.png_2"><img src="${prefix}-amp.png_2" width="400px"></a></td>'
                                + '<td></td>'
                                + '</tr>')
                    cal_idx.write(t.substitute(cal=cal, src_freq=src_freq, prefix=prefix[len(dirname) + 1:]))

    cal_idx.write('</table></body></html>\n')
    cal_idx.close()
    return error_list


def build_images(day_dir_name, sources, band, day):
    """
    Generate continuum images of each field.

    :param day_dir_name: The name of the day directory
    :param sources: The list of sources, which includes the pahse calibrator details
    :param band: The map describing the continuum band to be processed
    :param day: The day being processed
    :return: None
    """

    error_list = []
    os.chdir(day_dir_name)
    # We only work with the continuum data
    freq = band['main']
    magmo.ensure_dir_exists(freq)
    print "Producing %s MHz continuum images for %d sources" % (freq, len(sources))

    img_idx = open('images.html', 'w')
    t = Template('<html>\n<head><title>Image previews for day $day</title></head>\n'
                 + '<body>\n<h1>Image previews for day $day at $freq MHz</h1>\n<table>')
    img_idx.write(t.substitute(day=day, freq=freq))

    for src in sources:
        src_name = src['source']
        try:
            src_file, freq_suffix = find_freq_file(src_name, band['freqs'])
            if src_file is None:
                print "No continuum file found for source %s and band %s." % (src_name, freq)
                exit(1)

            phase_cal_file = src['phase_cal'] + "." + freq_suffix
            magmo.run_os_cmd('gpcopy vis=' + phase_cal_file + ' out=' + src_file)

            name_prefix = freq + '/magmo-' + src_name + '_' + freq
            dirty_file = name_prefix + '_dirty'
            clean_file = name_prefix + '_clean'
            beam_file = name_prefix + '_beam'
            restored_file = name_prefix + '_restor'
            fits_file = restored_file + '.fits'
            png_file = name_prefix + '.png'

            cmd = 'invert robust=0.5 options=systemp,mfs,double stokes=ii vis=' + src_file \
                + ' slop=1.0 ' \
                + ' map=' + dirty_file + ' beam=' + beam_file
            magmo.run_os_cmd(cmd)
            cmd = 'clean niters=2000 speed=+1 map=' + dirty_file \
                  + ' beam=' + beam_file + ' out=' + clean_file
            magmo.run_os_cmd(cmd)
            cmd = 'restor options=mfs model=' + clean_file + ' beam=' \
                + beam_file + ' map=' + dirty_file + ' out=' + restored_file
            magmo.run_os_cmd(cmd)
            cmd = 'fits op=xyout in=' + restored_file + ' out=' + fits_file
            magmo.run_os_cmd(cmd)

            fig = aplpy.FITSFigure(fits_file)
            fig.set_theme('publication')
            fig.show_grayscale(0, None, 0.25)
            fig.add_colorbar()
            fig.save(png_file)
            fig.close()
            t = Template('<tr><td><br>Source ${src_name}</td></tr>\n<tr>\n'
                         + '<td><a href="${png_file}"><img src="${png_file}" width="500px"></a></td></tr>')
            img_idx.write(t.substitute(png_file=png_file, src_name=src_name))

        except magmo.CommandFailedError as e:
            error_list.append(str(e))

    img_idx.write('</table></body></html>\n')
    img_idx.close()
    os.chdir('..')
    return error_list


def get_signal_noise_ratio(miriad_file):
    """
    Retrieve the signal to noise ratio for a miriad image file.
    :param miriad_file:
    :return: The root mean square of the data and the maximum flux density
    """
    result = subprocess.check_output(["imstat", "in="+miriad_file, "options=guaranteespaces"])
    total_next = False
    for line in result.splitlines():
        # print line
        if total_next:
            stats = line.strip().split()
            #print stats
            #mean = float(stats[1])
            rms = float(stats[2])
            max = float(stats[3])
            return rms, max
        if re.match('^[ \t]*Total', line):
            total_next = True
    print "Unable to find totals in:"
    print result


def find_strong_sources(day_dir_name, freq, sources, num_chan, min_sn):
    """
    Identify the source files with a maximum signal to noise which is greater
    than a threshold, that is those restored images which show a source
    sufficiently brighter than the background noise.

    :param day_dir_name: The name of the day directory
    :param freq: The primary identifier of the frequency to be checked.
    :param sources: The list of sources, which includes the pahse calibrator details
    :param num_chan: The number of channels in the line data
    :param min_sn: The minimum signal to noise ratio to be a 'strong' source.
    :return: A list of the source names filtered to just those with bright sources
    """
    strong_sources = []
    src_names = []
    src_rms = np.zeros(len(sources))
    src_max = np.zeros(len(sources))
    src_sn = np.zeros(len(sources))
    src_used = []

    i = 0
    for src in sources:
        src_name = src['source']
        restored_img = day_dir_name + "/" + freq + "/magmo-" + src_name + "_" + freq + "_restor"
        if os.path.exists(restored_img):
            rms, max = get_signal_noise_ratio(restored_img)
            sn = 0
            if rms > 0:
                sn = max / rms
            sn /= math.sqrt(num_chan)
            strong = sn > min_sn
            if strong:
                strong_sources.append(src)
            src_names.append(src_name)
            src_rms[i] = rms
            src_max[i] = max
            src_sn[i] = sn
            src_used.append("Y" if strong else "N")
            i += 1

    with open(day_dir_name+'/stats.csv', "wb") as stats:
        writer = csv.writer(stats, quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(["name", "rms", "max", "sn", "strong"])
        for i in range(0, len(src_names)):
            writer.writerow(
                [src_names[i], src_rms[i], src_max[i], src_sn[i], src_used[i]])
    #table = OrderedDict([('name', src_names), ('rms', src_rms), ('max', src_max), ('s/n', src_sn)])
    #ascii.write(table, day_dir_name+'/stats.dat', format='fixed_width', bookend=False, delimiter=None, quotechar='"')

    return strong_sources


def build_cubes(day_dir_name, sources, band):
    """
    Generate HI spectral cubes of each field.

    :param day_dir_name: The name of the day directory
    :param sources: The list of sources, which includes the pahse calibrator details
    :param band: The definition of the band being processed.
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
            line = 'felocity,627,-250.0,0.8,0.8'
            #line = 'felocity,1053,-250.0,0.4,0.4'

            cmd = 'uvaver line=' + line + ' vis=' + src_file + ' out=' + ave_file
            magmo.run_os_cmd(cmd)
            cmd = 'invert robust=0.5 cell=5 options=systemp,nopol,mosaic,double stokes=i '\
                + ' slop=1.0 line='+ line + ' vis=' + ave_file \
                + ' map=' + dirty_file + ' beam=' + beam_file
            magmo.run_os_cmd(cmd)
            cmd = 'clean niters=500 mode=steer speed=+1 map=' + dirty_file + ' beam=' \
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


def main():
    """
    Main script for process_data
    :return: None
    """
    # Read day parameter
    if len(sys.argv) != 2:
        print("Incorrect number of parameters.")
        print("Usage: python process_data.py day")
        exit(1)
    day = sys.argv[1]
    start = time.time()

    # metadata needed: flux/bandpass cal, phase cals for each source, extra flagging
    # Read metadata for the day (file pattern fragments etc)
    sources = magmo.get_day_obs_data(day)
    if sources is None or len(sources) == 0:
        print "Day %s is not defined." % (day)
        exit(1)
    print "#### Started processing MAGMO day %s at %s ####" % \
          (day, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))
    print "Found %d sources. First source:\n%s" % (len(sources), sources[0])
    error_list = []

    # Check metadata against file system
    dayDirName = "day" + day
    if not os.path.isdir(dayDirName):
        print "Directory %s could not be found." % dayDirName
        exit(1)

    # set up map of parent/child map of frequencies
    line_band = {'main': '1420', 'freqs': ['1420', '1421', '1420.5'], 'line': True}
    cont_band = {'main': '1757', 'freqs': ['1757'], 'line': False}
    band_list = [line_band, cont_band]
    for band in band_list:
        freq = band['main']
        freq_dir = dayDirName + "/" + freq
        magmo.ensure_dir_exists(freq_dir)

    # Flag
    bandpasscal = find_bandpasscal(dayDirName)
    print "Bandpass cal:", bandpasscal
    error_list.extend(flag_data(dayDirName, day, bandpasscal, sources))

    # Calibrate
    error_list.extend(calibrate(dayDirName, bandpasscal, sources, band_list, day))
    print error_list

    # Produce 2 GHz continuum image
    error_list.extend(build_images(dayDirName, sources, cont_band, day))

    # Produce HI image cube
    strong_sources = find_strong_sources(dayDirName, cont_band['main'], sources,
                                         num_chan, sn_min)
    print "### Found the following bright sources in the data ", strong_sources
    error_list.extend(build_cubes(dayDirName, strong_sources, line_band))

    # Report
    end = time.time()
    print '#### Processing Completed at %s ####' \
          % (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print 'Processed %d sources (%d strong enough to produce cubes) in %.02f s' % (len(sources),
                                                                                   len(strong_sources),
                                                                                   end - start)
    if len(error_list) == 0:
        print "Hooray! No errors found."
    else:
        print "%d errors were encountered:" % (len(error_list))
        for err in error_list:
            print err

# Run the script if it is called from the command line
if __name__ == "__main__":
    main()
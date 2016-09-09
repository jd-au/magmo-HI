# Quick extra program to proiduce images of the HI 1420 data and provide a side by side comparison with the
# 1757 MHz continuum data

import aplpy
import magmo
import math
import os
import process_data
import sys
import time

from string import Template

# Read day parameter
if len(sys.argv) != 2:
    print("Incorrect number of parameters.")
    print("Usage: python process_data.py day")
    exit(1)
day = sys.argv[1]

sources = magmo.get_day_obs_data(day)
if sources is None or len(sources) == 0:
    print "Day %s is not defined." % (day)
    exit(1)

dayDirName = "day" + day
start = time.time()
print "#### Started imaging MAGMO day %s at %s ####" % \
      (day, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

line_band = {'main': '1420', 'freqs': ['1420', '1421', '1420.5'], 'line': True}

error_list = []

error_list.extend(process_data.build_images(dayDirName, sources, line_band, day))

img_idx = open(dayDirName + '/images-comp.html', 'w')
t = Template('<html>\n<head><title>Day $day image previews</title></head>\n'
             + '<body>\n<h1>Image previews for day $day</h1>\n<table>\n'
             + '<tr><td>Source</td><td>1757 MHz</td><td>1420 MHz</td><td>1757 Beam</td></tr>')
img_idx.write(t.substitute(day=day))

os.chdir(dayDirName)
for src in sources:
    src_name = src['source']
    name_prefix = '/magmo-' + src_name + '_'
    restored_img = "1757/magmo-" + src_name + "_1757_restor"
    sn = 0.0
    if os.path.exists(restored_img):
        rms, max = process_data.get_signal_noise_ratio(restored_img)
        if rms > 0:
            sn = max / rms
        sn /= math.sqrt(1053)
    max_1420 = max
    restored_img = "1420/magmo-" + src_name + "_1420_restor"
    if os.path.exists(restored_img):
        rms_1420, max_1420 = process_data.get_signal_noise_ratio(restored_img)

    beam_file = "1757/magmo-" + src_name + "_1757_beam"
    beam_img = beam_file + ".png"
    if os.path.exists(beam_file):
        cmd = 'cgdisp region="percentage(15)" in=' + beam_file + ' device=' + beam_img\
              + '/png type=p'
        magmo.run_os_cmd(cmd)

    img_idx.write('<tr><td>' + src_name + '<br>S/N: ' + str(sn) + '<br>Max: ' + str(max)
                  + '<br>RMS: ' + str(rms) + '</td>\n')
    for freq in ['1757', '1420']:
        fits_file = freq + name_prefix + freq + '_restor.fits'
        img = freq + name_prefix + freq + '.png'
        #print 'File %s exists %s' % (fits_file, os.path.exists(fits_file))
        if os.path.exists(fits_file):
            fig = aplpy.FITSFigure(fits_file)
            fig.set_theme('publication')
            fig.show_grayscale(0, None, (max_1420 if freq == '1420' else max))
            # fig.show_grayscale(0, None, 0.1)
            fig.add_colorbar()
            fig.save(img)
            fig.close()

        t = Template('<td><a href="${img}"><img src="${img}" width="500px"></a></td>')
        img_idx.write(t.substitute(img=img))

    t = Template('<td><a href="${img}"><img src="${img}" width="500px"></a></td>')
    img_idx.write(t.substitute(img=beam_img))
    img_idx.write('</tr>\n')

img_idx.write('</table></body></html>\n')
img_idx.close()
os.chdir('..')

# Report
end = time.time()
print '#### Processing Completed at %s in %.02f s ####' \
      % (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)), end - start)

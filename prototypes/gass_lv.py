#!/usr/bin/env python -u

# Extract a longitude-velocity slice at galactic latitude = 0 from the
# GASS galactic projection cubes. This involves extracting the planes from each
# cube and stitching them together into a single image.

# Author James Dempsey
# Date 30 Oct 2016

from __future__ import print_function, division

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import glob
import numpy as np
import os


def copy_axis(in_header, in_axis, out_header, out_axis):
    for prefix in ['CTYPE', 'CRPIX', 'CDELT', 'CRVAL']:
        in_keyword = prefix+str(in_axis)
        out_keyword = prefix + str(out_axis)
        out_header[out_keyword] = in_header[in_keyword]


gass_folder = '/Volumes/Data-Dempsey/GASSIII/CAR_CUBES'

in_files = glob.glob(gass_folder + "/*.fits")
full_slice = None
max_vel = 0
gass_header = None
for input in in_files:
    print("Reading file", input)
    hdulist = fits.open(input, memmap=True)
    header = hdulist[0].header
    if gass_header is None:
        gass_header = header
    data = hdulist[0].data
    # print (data.shape)
    w_orig = WCS(header)
    pixcrd = np.array([[0, 0, 0], [0, 0, data.shape[0]]])
    #print (pixcrd.shape)
    #print (pixcrd)
    world = w_orig.wcs_pix2world(pixcrd, 0)
    print("Velocity range is %.4f to %.4f." % (world[0][2], world[1][2]))
    max_vel = int(world[1][2]/1000.0)
    lon, lat, vel = w_orig.wcs_world2pix(0, 0, 0, 0)
    #print (lat)
    #print (vel)
    # w_part = w_orig.slice()
    # w_lv = w_orig.swapaxes(2,1)
    # w_lv = w_lv.sub(2)

    #print (data.shape)
    #print (w_orig)
    lv_data = data[:, int(lat), :]
    # print(w_lv)
    #print(lv_data.shape)

    if full_slice is None:
        full_slice = lv_data
    else:
        full_slice = np.concatenate((full_slice, lv_data[1:]))
    print("Updated full slice:", full_slice.shape)

# Produce a fits cube
os.remove("gass-lv.fits")
hdu = fits.PrimaryHDU(full_slice)
copy_axis(gass_header, 1, hdu.header, 1)
copy_axis(gass_header, 2, hdu.header, 3)
copy_axis(gass_header, 3, hdu.header, 2)
hdu.header['CRPIX3']= 1
hdu.header['CRVAL3']= 0
for keyword in ['RESTFREQ', 'TELESCOP', 'SPECSYS', 'OBJECT', 'OBSERVER',
                'BUNIT', 'BSCALE', 'BZERO']:
    hdu.header[keyword] = gass_header[keyword]
hdu.writeto('gass-lv.fits')

# Plot the emission and its outline.
fig_filename = "full-lv.pdf"
print("Writing", fig_filename)
# fig = plt.figure()
# fig.add_subplot(111, projection=w_lv)
ax = plt.subplot(111)
#ax.imshow(np.log10(full_slice), origin='lower', cmap='viridis')
ax.imshow(full_slice, origin='lower', cmap='YlOrRd')
plt.xlabel('glon(deg)')
x_step = int(30*(full_slice.shape[1]/360))
print(x_step)
ax.set_xticks([i for i in range(0,full_slice.shape[1],x_step)])
ax.set_xticklabels([-i for i in range(-180,180,30)])

max_vel = 450
print(max_vel)
y_step = int(100*(full_slice.shape[0]/(500+max_vel)))
print(y_step)
ax.set_yticks([i for i in range(0,full_slice.shape[0],y_step)])
ax.set_yticklabels([i for i in range(-500, max_vel, 100)])
#plt.xlim([-180, 40])
# plt.axes.get_yaxis
# plt.x_scale()
plt.contour(np.log10(full_slice), 1, cmap='Blues')

plt.grid(color='darkgrey')

plt.ylabel('LSR velocity (km/s)')
plt.savefig(fig_filename)
plt.close()


# Plot just the contour of the GASS emission.
x_val = np.repeat(np.arange(0,full_slice.shape[0]), full_slice.shape[1])
y_val = np.tile(np.arange(0,full_slice.shape[1]), full_slice.shape[0])

data = np.reshape(full_slice, -1)

plt.figure()
#CS = plt.contour(x_val, y_val, data)
CS = plt.contour(np.log10(full_slice), 3)
# plt.clabel(CS, inline=1, fontsize=10)
plt.title('Contour plot of HI emission')
plt.savefig("contour.pdf")

# Prototype fucntion to extract a sp[ecyrum form a FITS format image cube
# from a specific spatial location
# Author: James Dempsey

from astropy.io import fits
from astropy.io import votable
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np


def readSources(filename):
    print "Extracting sources from " + filename
    sources = []
    src_votable = votable.parse(filename, pedantic=False)
    results = src_votable.get_first_table().array
    for row in results:
        ra = row['ra']
        dec = row['dec']
        rms = row['local_rms']
        flux = row['peak_flux']
        sn = flux / rms
        print "Found source at %.4f, %.4f with flux %.4f and rms of %.4f giving S/N of %.4f" % \
              (ra, dec, flux, rms, sn)
        if (sn > 10):
            sources.append([ra, dec])
        else:
            print "Ignoring source at %.4f, %.4f due to low S/N of %.4f" % \
                  (ra, dec, sn)

    return sources


def plotSpectrum(x, y):
    """
    Output a plot of opacity vs LSR velocity to a specified file.

    :param x: The velocity data
    :param y: The opacity values for each velocity step
    :param filename: The file the plot should be written to. Should be an .eps or .pdf file.
    """
    fig = plt.figure()
    #ax = fig.add_subplot()
    #pylab.figure(fig)

    #pylab.plot(x, y,
    #           color='Blue',linestyle='--',label='Separation distance')
    plt.plot(x, y)

    #plt.ylim([-0.01, 0.01])
    #plt.set_autoscaley_on(False)
    #plt.yscale(-0.1, 0.1)
    #plt.axhline(1)

    plt.xlabel(r'Velocity relative to LSR (km/s)')
    plt.ylabel(r'$e^{(-\tau)}$')
    plt.grid(True)
    #plt.savefig(filename)
    plt.show()
    return


#fits_filename = 'day21/1420/magmo-270.255+0.835_1420_sl_restor.fits'
#fits_filename = 'day21/1420/magmo-291.270-0.719_1420_sl_restor.fits'
field = '285.337-0.002'
daydirname = 'day21'
fits_filename = "{0}/1420/magmo-{1}_1420_sl_restor.fits".format(daydirname, field)
src_filename = "{0}/{1}_src_comp.vot".format(daydirname, field)
ra = 157.87033137328115
dec = -58.04001286446428


#sources = readSources('day21/285.337-0.002_src_comp.vot')
sources = readSources(src_filename)
hdulist_im = fits.open(fits_filename)
image = hdulist_im[0].data
header = hdulist_im[0].header
#header = fits.getheader(fits_filename)
w = WCS(header)
for src in sources:
    pix = w.wcs_world2pix(src[0], src[1], 0, 0, 1)
    x_coord = int(round(pix[0]))-1 # 266
    y_coord = int(round(pix[1]))-1 # 197
    print "Translated %.4f, %.4f to %d, %d" % (src[0], src[1], x_coord, y_coord)

    slice = image[0,:,y_coord, x_coord]

    plotSpectrum(np.arange(slice.size), slice)

x_size = image.shape[2]
y_size = image.shape[3]
continuum = np.zeros((x_size, y_size))

for i in range(0,x_size-1):
    for j in range(0,y_size-1):
        slice = image[0,28:228,i,j]
        continuum[i,j] = max(0,np.mean(slice))
        if continuum[i,j] > 0.1:
            print i, j, continuum[i,j]
            if i > 1 and j > 1:
                slice = image[0, :, i, j]
                #plotSpectrum(np.arange(slice.size), slice)

plt.imshow(continuum, cmap='gray')
plt.colorbar()
plt.show()



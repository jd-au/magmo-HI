#!/usr/bin/env python -u

# Decompose each of the MAGMO spectra into Gaussian components.
#
# This program reads the previously generated spectra summaries and uses GaussPy to fit Gaussian components to
# the spectra. Comparison diagrams are produced for each spectrum.

# Author James Dempsey
# Date 5 Jan 2017

from __future__ import print_function, division

from astropy.io import fits, votable
from astropy.io.votable import parse, from_table, writeto
from astropy.table import Table, Column
from string import Template

import argparse
import datetime
import gausspy.gp as gp
import numpy as np
import pickle
import time
import matplotlib.pyplot as plt


FILENAME_DATA_GAUSSPY = 'spectra.pickle'


def parseargs():
    """
    Parse the command line arguments
    :return: An args map with the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Decompose Galactic plane HI absorption spectra into Gaussian components")
    parser.add_argument("-i", "--input", help="The input spectra catalogue",
                        default='magmo-spectra.vot')
    parser.add_argument("--long_min", help="The smallest longitude to be decomposed",
                        type=int, default=-180)
    parser.add_argument("--long_max", help="The largest longitude to be decomposed",
                        type=int, default=180)
    parser.add_argument("-o", "--output", help="The file name of the decomposition result.",
                        default='magmo-decomp.pickle')
    parser.add_argument("--plot_only", help="Produce plots for the result of a previous decomposition", default=False,
                        action='store_true')
    parser.add_argument("--train", help="Train GaussPy using the selected spectra. The produced alpha value can then " +
                                        "be used for later decompostion.", default=False, action='store_true')
    parser.add_argument("--alpha", help="The value for the GaussPy smoothing parameter",
                        type=float, default=1.12)
    parser.add_argument("--snr_thresh", help="The signal to noise ratio threshold",
                        type=float, default=0.2)

    args = parser.parse_args()
    return args


def read_spectra(filename):
    votable = parse(filename, pedantic=False)
    results = next(resource for resource in votable.resources if
                   resource.type == "results")
    results_array = results.tables[0].array
    return results_array


def read_opacity(filename):
    votable = parse(filename, pedantic=False)
    results = next(resource for resource in votable.resources if
                   resource.type == "results")
    results_array = results.tables[0].array
    return results_array


def filter_spectra(spectra, min_long, max_long):
    filtered = spectra[spectra['Longitude'] >= min_long]
    filtered = filtered[filtered['Longitude'] <= max_long]
    filtered = filtered[filtered['Rating'] <= 'C']
    return filtered


def get_opacity_filename(day, field, source):
    t = Template('day${day}/${field}_src${source}_opacity.votable.xml')
    return t.substitute(day=day, field=field, source=source)


def prepare_spectra(spectra):
    data = {}

    # TODO: Convert this to use a real error value
    rms = 0.06

    # Convert to GaussPy format
    i = 0
    for spectrum in spectra:
        filename = get_opacity_filename(spectrum['Day'], spectrum['Field'], spectrum['Source'])
        print ("Reading", filename)
        opacity = read_opacity(filename)

        errors = np.ones(opacity.shape[0]) * rms
        location = np.array(spectrum['Longitude'], spectrum['Latitude'])
        inverted_opacity= 1- opacity['opacity']

        data['data_list'] = data.get('data_list', []) + [inverted_opacity]
        data['x_values'] = data.get('x_values', []) + [opacity['velocity']/1000]
        data['errors'] = data.get('errors', []) + [errors]
        data['location'] = data.get('location', []) + [location]
        i+= 1

    # Save the file to be used by GaussPy
    pickle.dump(data, open(FILENAME_DATA_GAUSSPY, 'w'))


def decompose(spectra, out_filename, alpha1, snr_thresh):
    start = time.time()
    print("## Commenced decomposition at %s ##" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    prepare_spectra(spectra)

    end_read = time.time()
    print("## Finished conversion of %d spectra at %s taking %.02f s ##" %
          (len(spectra), time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_read)), (end_read - start)))

    # Load GaussPy
    g = gp.GaussianDecomposer()

    # Setting AGD parameters
    g.set('phase', 'one')
    g.set('SNR_thresh', [snr_thresh, snr_thresh])
    g.set('alpha1', alpha1)

    # Run GaussPy
    decomposed_data = g.batch_decomposition(FILENAME_DATA_GAUSSPY)

    # Save decomposition information
    pickle.dump(decomposed_data, open(out_filename, 'w'))

    end = time.time()
    print("## Finished decomposition of %d spectra at %s taking %.02f s ##" %
          (len(spectra), time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)), (end - end_read)))


def train(spectra, alpha1_initial, snr_thresh):
    start = time.time()
    print("## Commenced training at %s ##" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    prepare_spectra(spectra)

    end_read = time.time()
    print("## Finished conversion of %d spectra at %s taking %.02f s ##" %
          (len(spectra), time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_read)), (end_read - start)))

    g = gp.GaussianDecomposer()
    # Next, load the training dataset for analysis:
    g.load_training_data(FILENAME_DATA_GAUSSPY)
    # Set GaussPy parameters
    g.set('phase', 'one')
    g.set('SNR_thresh', [snr_thresh, snr_thresh])

    # Train AGD starting with initial guess for alpha
    g.train(alpha1_initial=alpha1_initial)

    end = time.time()
    print("## Finished training using %d spectra at %s taking %.02f s ##" %
          (len(spectra), time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)), (end - end_read)))


def output_component_catalogue(spectra, data_decomposed):
    days = []
    field_names = []
    sources = []
    longitudes = []
    latitudes = []
    amps = []
    fwhms = []
    means = []

    num_no_comps = {}

    for i in range(len(data_decomposed['fwhms_fit'])):
        spectrum = spectra[i]
        fit_fwhms = data_decomposed['fwhms_fit'][i]
        fit_means = data_decomposed['means_fit'][i]
        fit_amps = data_decomposed['amplitudes_fit'][i]
        if len(fit_amps) > 0.:
            for j in range(len(fit_amps)):
                days.append(int(spectrum['Day']))
                field_names.append(spectrum['Field'])
                sources.append(spectrum['Source'])
                longitudes.append(spectrum['Longitude'])
                latitudes.append(spectrum['Latitude'])
                amps.append(fit_amps[j])
                fwhms.append(fit_fwhms[j])
                means.append(fit_means[j])
        else:
            rating = spectrum['Rating']
            num_no_comps[rating] = num_no_comps.get(rating, 0) + 1


    temp_table = Table(
        [days, field_names, sources, longitudes, latitudes, amps, fwhms, means],
        names=['Day', 'Field', 'Source', 'Longitude', 'Latitude', 'Amplitude', 'FWHM', 'Mean'],
        meta={'ID': 'magmo_components',
              'name': 'MAGMO Components ' + str(datetime.date.today())})
    votable = from_table(temp_table)
    filename = "magmo-components.vot"
    writeto(votable, filename)

    print("Wrote out", len(fwhms), "components to", filename, "No components generated for")
    for rating in sorted(num_no_comps.keys()):
        print ("%s: %3d" % (rating, num_no_comps[rating]))
    # , ".", num_no_comps, "spectra (of", len(spectra), ") had no components found")


def output_decomposition(spectra, out_filename):

    # For each spectra:
    # Output plots
    # Output component catalogue

    max_results = 6 if len(spectra) > 6 else len(spectra)
    print (max_results)

    data = pickle.load(open(FILENAME_DATA_GAUSSPY))
    data_decomposed = pickle.load(open(out_filename))
    index_values = np.argsort(np.random.randn(5000))  # plot random results
    fig = plt.figure(0, [12, 12])
    for i in range(max_results):
        ax = fig.add_subplot(4, 3, i+1+((i//3)*3)) #, sharex=True)
        index = i #index_values[i]
        x = data['x_values'][index]
        y = 1-data['data_list'][index]
        fit_fwhms = data_decomposed['fwhms_fit'][index]
        fit_means = data_decomposed['means_fit'][index]
        fit_amps = data_decomposed['amplitudes_fit'][index]
        g_sum = np.zeros(len(x))
        ax.plot(x, y, color='grey')
        # Plot individual components
        if len(fit_amps) > 0.:
            for j in range(len(fit_amps)):
                amp, fwhm, mean = fit_amps[j], fit_fwhms[j], fit_means[j]
                yy = amp * np.exp(-4. * np.log(2) * (x - mean) ** 2 / fwhm ** 2)
                g_sum += yy
                yy = 1- yy
                ax.plot(x, yy, '-', lw=1.5, color='purple')

        ax.plot(x, 1-g_sum, '--', lw=1.0, color='blue')

        if i % 3 == 0:
            ax.set_ylabel('$1 - e^{-\\tau}$')

        # Residual plot
        residual = y - (1-g_sum)
        ax = fig.add_subplot(4, 3, i+4+((i//3)*3))
        #frame2 = ax.add_axes((.1, .1, .8, .2))
        ax.plot(x, residual, 'or', markerfacecolor='None', markersize=2, markeredgecolor='blue')
        ax.grid()

        #ax.set_xlim(400, 600)
        if i >= 3:
            ax.set_xlabel('LSR Velocity (km/s)')
    plt.savefig("magmo-decomp.pdf")

    output_component_catalogue(spectra, data_decomposed)


def main():
    # Parse command line options
    args = parseargs()

    start = time.time()
    print("#### Started processing MAGMO spectra at %s ####" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    # Read in spectra
    spectra = read_spectra(args.input)
    spectra = filter_spectra(spectra, args.long_min, args.long_max)

    # Decompose all spectra
    if args.train:
        train(spectra, args.alpha, args.snr_thresh)
    elif not args.plot_only:
        decompose(spectra, args.output, args.alpha, args.snr_thresh)

    # Read in result
    if not args.train:
        output_decomposition(spectra, args.output)

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Processed spectra in %.02f s' %
          (end - start))
    return 0


if __name__ == '__main__':
    exit(main())

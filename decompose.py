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
from matplotlib import gridspec
from string import Template

import argparse
import datetime
import gausspy.gp as gp
import magmo
import numpy as np
import pickle
import time
import matplotlib.pyplot as plt
import aplpy

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
    parser.add_argument("-q", "--quality", help="The minimum quality level to include",
                        default='C')
    parser.add_argument("-o", "--output", help="The file name of the decomposition result.",
                        default='magmo-decomp.pickle')
    parser.add_argument("--plot_only", help="Produce plots for the result of a previous decomposition", default=False,
                        action='store_true')
    parser.add_argument("--train", help="Train GaussPy using the selected spectra. The produced alpha value can then " +
                                        "be used for later decomposition.", default=False, action='store_true')
    parser.add_argument("--alpha1", help="The value for the first GaussPy smoothing parameter",
                        type=float, default=1.12)
    parser.add_argument("--alpha2", help="The value for the second GaussPy smoothing parameter",
                        type=float, default=5)
    parser.add_argument("--snr_thresh", help="The signal to noise ratio threshold",
                        type=float, default=5)

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


def filter_spectra(spectra, min_long, max_long, min_quality):
    filtered = spectra[spectra['Longitude'] >= min_long]
    filtered = filtered[filtered['Longitude'] <= max_long]
    filtered = filtered[filtered['Rating'] <= min_quality]
    filtered = filtered[filtered['Duplicate'] == False]
    return filtered


def get_opacity_filename(day, field, source):
    t = Template('day${day}/${field}_src${source}_opacity.votable.xml')
    return t.substitute(day=day, field=field, source=source)


def convert_from_ratio(absorption):
    """
    Convert an array of absorption values (I/I_0) to opacity values (tau).
    Note that this currently clips any negative values to avoid nan in the
    result.
    :param absorption: The array of absorption values
    :return: The equivalent tau values
    """
    #return -1 * np.log(np.maximum(absorption, 1e-16))
    return 1 - absorption


def convert_to_ratio(opacity):
    """
    Convert an array of opacity values (tau) to absorption values (I/I_0).
    :param opacity: The array of tau values
    :return: The equivalent absorption ratio values
    """
    #return np.exp(-1 * opacity)
    return 1 - opacity


def prepare_spectra(spectra, data_filename):
    data = {}

    # Convert to GaussPy format
    i = 0
    for spectrum in spectra:
        filename = get_opacity_filename(spectrum['Day'], spectrum['Field'], spectrum['Source'])
        #print("Reading", filename)
        opacity = read_opacity(filename)
        #if spectrum['Day'] == 43 and spectrum['Rating'] == 'A' and spectrum['Field'] == '019.612-0.120': # and spectrum['Source'] == '25-0':
        #    print("Skipping ", spectrum['Rating'], spectrum['Day'], spectrum['Longitude'], spectrum['Field'], spectrum['Source'])
        #    continue

        longitude = spectrum['Longitude']
        if longitude < 0:
            longitude += 360
        print (i, spectrum['Rating'], spectrum['Day'], longitude, spectrum['Field'], spectrum['Source'])
        rms = spectrum['Continuum_SD']

        errors = np.ones(opacity.shape[0]) * rms
        location = np.array(spectrum['Longitude'], spectrum['Latitude'])
        tau = convert_from_ratio(opacity['opacity'])
        #print (tau)
        #inverted_opacity = 1 - opacity['opacity']

        data['data_list'] = data.get('data_list', []) + [tau]
        data['x_values'] = data.get('x_values', []) + [opacity['velocity'] / 1000]
        data['errors'] = data.get('errors', []) + [errors]
        data['location'] = data.get('location', []) + [location]
        data['spectrum_idx'] = data.get('spectrum_idx', []) + [i]
        data['rating'] = data.get('rating', []) + [spectrum['Rating']]
        i += 1

    # Save the file to be used by GaussPy
    pickle.dump(data, open(data_filename, 'w'))


def decompose(spectra, out_filename, alpha1, alpha2, snr_thresh, data_filename):
    start = time.time()
    print("## Commenced decomposition at %s ##" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    prepare_spectra(spectra, data_filename)

    end_read = time.time()
    print("## Finished conversion of %d spectra at %s taking %.02f s ##" %
          (len(spectra), time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_read)), (end_read - start)))

    # Load GaussPy
    g = gp.GaussianDecomposer()

    # Setting AGD parameters
    g.set('mode','conv')
    g.set('phase', 'two' if alpha2 else 'one')
    g.set('SNR_thresh', [snr_thresh, snr_thresh])
    g.set('alpha1', alpha1)
    #g.set('verbose', True)
    if alpha2:
        g.set('alpha2', alpha2)

    # Run GaussPy
    decomposed_data = g.batch_decomposition(data_filename)

    # Save decomposition information
    pickle.dump(decomposed_data, open(out_filename, 'w'))

    end = time.time()

    if len(decomposed_data['means_fit']) != len(spectra):
        print(
            '###! WARNING: Original %d and decomposed spectra %d counts differ!' % (
            len(spectra), len(decomposed_data['means_fit'])))

    print("## Finished decomposition of %d spectra at %s taking %.02f s ##" %
          (len(spectra), time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)), (end - end_read)))


def output_component_catalogue(spectra, data, data_decomposed):
    days = []
    field_names = []
    sources = []
    longitudes = []
    latitudes = []
    amps = []
    fwhms = []
    means = []
    best_fit_rchi2s = []
    amps_fit_errs = []
    fwhms_fit_errs = []
    means_fit_errs = []


    num_no_comps = {}

    for i in range(len(data_decomposed['fwhms_fit'])):
        if i >= len(data['spectrum_idx']):
            print("Error: data index of %d is invalid for data array of len %d" % (
                i, len(data['spectrum_idx'])))
        spectrum_idx = data['spectrum_idx'][i]
        if spectrum_idx >= len(spectra):
            print("Error: spectra index of %d at row %d is invalid for spectra array of len %d" % (
                spectrum_idx, i, len(spectra)))
        spectrum = spectra[spectrum_idx]
        fit_fwhms = data_decomposed['fwhms_fit'][i]
        fit_means = data_decomposed['means_fit'][i]
        fit_amps = data_decomposed['amplitudes_fit'][i]
        best_fit_rchi2 = data_decomposed['best_fit_rchi2'][i]
        means_fit_err = data_decomposed['means_fit_err'][i]
        fwhms_fit_err = data_decomposed['fwhms_fit_err'][i]
        amplitudes_fit_err = data_decomposed['amplitudes_fit_err'][i]

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
                best_fit_rchi2s.append(best_fit_rchi2[0])
                amps_fit_errs.append(means_fit_err[j])
                fwhms_fit_errs.append(fwhms_fit_err[j])
                means_fit_errs.append(amplitudes_fit_err[j])
        else:
            rating = spectrum['Rating']
            num_no_comps[rating] = num_no_comps.get(rating, 0) + 1

    temp_table = Table(
        [days, field_names, sources, longitudes, latitudes, amps, fwhms, means, best_fit_rchi2s, amps_fit_errs,
         fwhms_fit_errs, means_fit_errs],
        names=['Day', 'Field', 'Source', 'Longitude', 'Latitude', 'Amplitude', 'FWHM', 'Mean', 'Best Fit Rchi2',
               'Amplitude Fit Err', 'FWHM Fit Err', 'Mean Fit Err'],
        meta={'ID': 'magmo_components',
              'name': 'MAGMO Components ' + str(datetime.date.today())})
    votable = from_table(temp_table)
    filename = "magmo-components.vot"
    writeto(votable, filename)

    total_nc = 0
    for rating, count in num_no_comps.items():
        total_nc += count

    print("Wrote out", len(fwhms), "components to", filename, "No components generated for", total_nc)
    for rating in sorted(num_no_comps.keys()):
        print("%s: %3d" % (rating, num_no_comps[rating]))
        # , ".", num_no_comps, "spectra (of", len(spectra), ") had no components found")


def calc_residual(velo, opacity, fit_amps, fit_fwhms, fit_means):
    g_sum = np.zeros(len(velo))
    # Plot individual components
    if len(fit_amps) > 0.:
        for j in range(len(fit_amps)):
            amp, fwhm, mean = fit_amps[j], fit_fwhms[j], fit_means[j]
            yy = amp * np.exp(-4. * np.log(2) * (velo - mean) ** 2 / fwhm ** 2)
            g_sum += yy
    residual = opacity - g_sum
    return residual


def output_decomposition_catalogue(folder, spectra, data, data_decomposed):
    days = []
    field_names = []
    sources = []
    longitudes = []
    latitudes = []
    cont_sd = []
    residual_rms = []
    ratings = []
    num_comps = []

    for i in range(len(data_decomposed['fwhms_fit'])):
        spectrum = spectra[data['spectrum_idx'][i]]

        days.append(int(spectrum['Day']))
        field_names.append(spectrum['Field'])
        sources.append(spectrum['Source'])
        longitudes.append(spectrum['Longitude'])
        latitudes.append(spectrum['Latitude'])
        ratings.append(spectrum['Rating'])
        cont_sd.append(spectrum['Continuum_SD'])

        fit_fwhms = data_decomposed['fwhms_fit'][i]
        fit_means = data_decomposed['means_fit'][i]
        fit_amps = data_decomposed['amplitudes_fit'][i]
        num_comps.append(len(fit_fwhms))

        velo = data['x_values'][i]
        # opacity = convert_to_ratio(data['data_list'][i])
        residual = calc_residual(velo, data['data_list'][i], fit_amps, fit_fwhms, fit_means)
        residual_rms.append(np.sqrt(np.mean(np.square(residual))))

    temp_table = Table(
        [days, field_names, sources, longitudes, latitudes, residual_rms, ratings, num_comps, cont_sd],
        names=['Day', 'Field', 'Source', 'Longitude', 'Latitude', 'Residual_RMS', 'Rating', 'Num_Comp', 'Continuum_SD'],
        meta={'ID': 'magmo_decomposition',
              'name': 'MAGMO Decomposition ' + str(datetime.date.today())})
    votable = from_table(temp_table)
    filename = folder + "/magmo-decomposition.vot"
    writeto(votable, filename)


def plot_single_spectrum(ax, velo, opacity, fit_amps, fit_fwhms, fit_means, name):
    g_sum = np.zeros(len(velo))
    ax.plot(velo, opacity, color='grey')
    # Plot individual components
    if len(fit_amps) > 0.:
        for j in range(len(fit_amps)):
            amp, fwhm, mean = fit_amps[j], fit_fwhms[j], fit_means[j]
            yy = amp * np.exp(-4. * np.log(2) * (velo - mean) ** 2 / fwhm ** 2)
            g_sum += yy
            #yy = convert_to_ratio(yy)
            ax.plot(velo, yy, '--', lw=0.5, color='purple')
    #g_sum = convert_to_ratio(g_sum)
    ax.plot(velo, g_sum, '-', lw=1.0, color='blue')
    plt.title(name)
    residual = opacity - g_sum
    return residual


def plot_spectrum(velo, opacity, fit_amps, fit_fwhms, fit_means, name, filename, formats):
    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax = fig.add_subplot(gs[0])

    y = convert_from_ratio(opacity)
    residual = plot_single_spectrum(ax, velo, y, fit_amps, fit_fwhms, fit_means, name)
    residual_rms = np.sqrt(np.mean(np.square(residual)))
    ax.set_ylabel('$1 - e^{-\\tau}$')

    # Residual plot
    ax = fig.add_subplot(gs[1])
    ax.plot(velo, residual, 'or', markerfacecolor='None', markersize=2, markeredgecolor='blue')
    ax.grid()

    plt.xlabel('LSR Velocity (km/s) n=%d rms=%.4f' % (len(fit_amps), residual_rms))

    for fmt in formats:
        plt.savefig(filename + "." + fmt)
    plt.close()
    return residual_rms


def plot_spectra(spectra, data, data_decomposed, alpha1, alpha2, folder='.'):
    plot_folder = folder + "/plots"
    magmo.ensure_dir_exists(plot_folder)
    residuals = np.zeros(len(data_decomposed['fwhms_fit']))
    for rating in 'ABCDEF':
        magmo.ensure_dir_exists(plot_folder + "/" + rating)

    for i in range(len(data_decomposed['fwhms_fit'])):
        spectrum = spectra[data['spectrum_idx'][i]]
        velo = data['x_values'][i]
        #opacity = 1 - data['data_list'][i]
        opacity = convert_to_ratio(data['data_list'][i])
        fit_fwhms = data_decomposed['fwhms_fit'][i]
        fit_means = data_decomposed['means_fit'][i]
        fit_amps = data_decomposed['amplitudes_fit'][i]
        rating = spectrum['Rating']
        field_name = spectrum['Field']
        day = str(spectrum['Day'])
        src_id = spectrum['Source']
        name = field_name + " src " + src_id + " on day " + day + "(" + rating + ")"
        filename = plot_folder + "/" + rating + "/"
        filename += field_name + "_" + day + "_src" + src_id + "_fit"
        residuals[i] = plot_spectrum(velo, opacity, fit_amps, fit_fwhms, fit_means, name, filename, ('pdf', 'png'))
    print("Residual RMS mean {:0.4} median {:0.4} sd {:0.4} for alphas {},{}".format(np.mean(residuals), np.median(residuals),
                                                                           np.std(residuals), alpha1, alpha2))


def plot_components(spectra, data_decomposed):
    fits_filename = 'hi4pi-total-car.fits'
    hdulist = fits.open(fits_filename, memmap=True)
    image = hdulist[0].data
    header = hdulist[0].header
    fig = aplpy.FITSFigure(image, header)
    fig.set_theme('publication')
    # fig.show_colorscale(vmin=vmin, vmax=np.max(image), cmap="jet")
    fig.add_colorbar()
    fig.save('magmo-map-comp.pdf')
    fig.close()


def get_samples(data, rating_count=[4, 1, 1], ratings='ABC'):
    indexes = []
    found_spectra = [0,0,0]
    print (len(data['rating']))
    for i in range(len(data['rating'])):
        print (data['rating'][i])
        key = ratings.index(data['rating'][i])
        if key >= 0 and found_spectra[key] < rating_count[key]:
            indexes.append(i)
            found_spectra[key] += 1
    return indexes


def output_decomposition(spectra, out_filename, folder, data_filename, alpha1, alpha2):
    # For each spectra:
    # Output plots
    # Output component catalogue

    max_results = 6 if len(spectra) > 6 else len(spectra)
    print(max_results)

    data = pickle.load(open(data_filename))
    data_decomposed = pickle.load(open(out_filename))

    index_values = get_samples(data)
    print (index_values, len(data_decomposed['fwhms_fit']), len(spectra))

    fig = plt.figure(0, [12, 12])
    i = 0
    for index in index_values:
        ax = fig.add_subplot(4, 3, i + 1 + ((i // 3) * 3))  # , sharex=True)
        #index = i  # index_values[i]
        x = data['x_values'][index]
        #y = 1 - data['data_list'][index]
        #y = convert_to_ratio(data['data_list'][index])
        y = data['data_list'][index]
        fit_fwhms = data_decomposed['fwhms_fit'][index]
        fit_means = data_decomposed['means_fit'][index]
        fit_amps = data_decomposed['amplitudes_fit'][index]

        spectrum = spectra[data['spectrum_idx'][index]]
        rating = spectrum['Rating']
        field_name = spectrum['Field']
        day = str(spectrum['Day'])
        src_id = spectrum['Source']

        name = field_name + " src " + src_id + " d " + day + "(" + rating + ")"

        residual = plot_single_spectrum(ax, x, y, fit_amps, fit_fwhms, fit_means, name)
        residual_rms = np.sqrt(np.mean(np.square(residual)))
        print (name, "has residual RMS of", residual_rms)

        if i % 3 == 0:
            ax.set_ylabel('$1 - e^{-\\tau}$')

        # Residual plot
        ax = fig.add_subplot(4, 3, i + 4 + ((i // 3) * 3))
        # frame2 = ax.add_axes((.1, .1, .8, .2))
        ax.plot(x, residual, 'or', markerfacecolor='None', markersize=2, markeredgecolor='blue')
        ax.grid()

        # ax.set_xlim(400, 600)
        if i >= 3:
            ax.set_xlabel('LSR Velocity (km/s)')
        i += 1

    plt.tight_layout()
    magmo.ensure_dir_exists(folder)
    plt.savefig(folder + "/magmo-decomp.pdf")
    plt.close()

    output_component_catalogue(spectra, data, data_decomposed)
    output_decomposition_catalogue(folder, spectra, data, data_decomposed)
    plot_spectra(spectra, data, data_decomposed, alpha1, alpha2, folder=folder)


def main():
    # Parse command line options
    args = parseargs()

    start = time.time()
    print("#### Started processing MAGMO spectra at %s ####" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    # Read in spectra
    spectra = read_spectra(args.input)
    spectra = filter_spectra(spectra, args.long_min, args.long_max, args.quality)

    alpha1_range = (4.36, 3)
    alpha2_range = (9.37, 7)

    #for i in range(0,1):  #len(alpha1_range)):
    for i in range(len(alpha1_range)):
        a1 = alpha1_range[i]
        a2 = alpha2_range[i]
        folder = "run"+ str(i+1)
        magmo.ensure_dir_exists(folder)
        data_filename = folder + "/" + FILENAME_DATA_GAUSSPY

        # Decompose all spectra
        if not args.plot_only:
            decompose(spectra, folder+"/"+args.output, a1, a2, args.snr_thresh, data_filename)

        # Read in result
        output_decomposition(spectra, folder+"/"+args.output, folder, data_filename, a1, a2)

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Processed spectra in %.02f s' %
          (end - start))
    return 0


if __name__ == '__main__':
    exit(main())

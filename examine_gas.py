#!/usr/bin/env python -u

# Take the Gaussian components and use them to examine the gas characteristics at each position.
#

# Author James Dempsey
# Date 26 Mar 2017


from __future__ import print_function, division

import argparse
import csv
import datetime
import glob
import math
import os
import re
from string import Template
import time

from astropy.coordinates import SkyCoord, matching
from astropy.io.votable import parse, from_table, writeto
from astropy.table import Table, Column, hstack, vstack
from astropy.io import ascii
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter
import numpy as np
import numpy.ma as ma
from scipy import interpolate
import scipy.stats as stats
import seaborn as sns

import magmo


class Gas(object):
    def __init__(self, day, field, src):
        self.day = day
        self.field = field
        self.src = src


def parseargs():
    """
    Parse the command line arguments
    :return: An args map with the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Examine the HI gas represented by each Gaussian component")

    parser.add_argument("--compare_only", help="Only run the comparison steps", default=False,
                        action='store_true')
    parser.add_argument("--no_compare", help="Only run the analysis steps", default=False,
                        action='store_true')

    args = parser.parse_args()
    return args


def read_votable_results(filename):
    votable = parse(filename, pedantic=False)
    results = next(resource for resource in votable.resources if
                   resource.type == "results")
    results_array = results.tables[0].array
    return results_array


def read_components(filename):
    return read_votable_results(filename)


def get_spectra_key(day, field, source):
    return '{}|{}|{}'.format(day, field, source)


def read_spectra(filename):
    spectra_array = read_votable_results(filename)
    spectra_array = spectra_array[spectra_array['Duplicate'] == False]
    spectra_map = {}
    for row in spectra_array:
        key = get_spectra_key(row['Day'], row['Field'], row['Source'])
        spectra_map[key] = row
    return spectra_map, spectra_array


def read_emission(filename):
    return read_votable_results(filename)


def get_field_key(name):
    key = name.strip()
    if re.search('^[0-9]\.', key):
        key = '00' + key
    elif re.search('^[0-9][0-9]\.', key):
        key = '0' + key
    return key


def read_mmb_cat(filename):
    """
    Read in the catalogue for the 6 GHz Methanol Multibeam Maser catalogue (2010MNRAS.404.1029C)
    :param filename: The filename of the votable catalogue.
    :return: A map of catalogue rows indexed by their field name keys
    """
    mmb_cat = read_votable_results(filename)
    maser_map = {}
    for row in mmb_cat:
        key = get_field_key(row['Name'])
        maser_map[key] = row
    return maser_map


def get_emission_filename(day, field, source):
    t = Template('day${day}/${field}_src${source}_emission.votable.xml')
    return t.substitute(day=day, field=field, source=source)


def get_temp(emission, comp_vel):
    velocities = emission['velocity']
    for i in range(0, len(velocities)):
        if velocities[i] >= comp_vel:
            temp = emission['em_mean'][i]
            temp_err = emission['em_std'][i]
            return temp, temp_err, velocities[i]
    return 0, 0, 0


def analyse_components(components, spectra_map, mmb_map):
    all_gas = []
    for component in components:

        comp_vel = component['Mean']
        comp_width = component['FWHM']
        comp_amp = component['Amplitude']
        optical_depth = 1 - comp_amp  # From 1-e^-tau to e^-tau

        gas = Gas(component['Day'], component['Field'], component['Source'])
        gas.comp_vel = comp_vel
        gas.comp_width = math.fabs(comp_width)
        gas.optical_depth = optical_depth
        gas.longitude = component['Longitude']
        gas.latitude = component['Latitude']
        gas.tau = -1 * np.log(np.maximum(optical_depth, 1e-16))
        gas.t_kmax = 21.866 * comp_width ** 2
        print("Gas T_kmax=", gas.t_kmax, "comp_width=", comp_width)
        gas.t_off = None
        gas.t_s = None
        gas.em_vel = None
        gas.name = component['Comp_Name']
        gas.spectra_name = component['Spectra_Name']

        spectrum = spectra_map[get_spectra_key(component['Day'], component['Field'], component['Source'])]
        gas.rating = spectrum['Rating']
        gas.continuum_sd = spectrum['Continuum_SD']

        loc = SkyCoord(gas.longitude, gas.latitude, frame='galactic', unit="deg")
        gas.loc = loc
        gas.ra = loc.icrs.ra.degree
        gas.dec = loc.icrs.dec.degree

        maser = mmb_map.get(component['Field'])
        if maser is None:
            print("unable to find maser for " + component['Field'])
        else:
            gas.maser_vel_low = maser['VL']
            gas.maser_vel_high = maser['VH']
            gas.maser_loc = SkyCoord(maser['RAJ2000'], maser['DEJ2000'], frame='fk5', unit="deg")

        # Load emission data
        t_off = 0
        emission_filename = get_emission_filename(component['Day'], component['Field'], component['Source'])
        if os.path.exists(emission_filename):
            emission = read_emission(emission_filename)
            t_off, t_off_err, em_vel = get_temp(emission, comp_vel * 1000)

            gas.em_vel = em_vel


            # Validate the component velocity
            if not spectrum['Min_Velocity'] <= comp_vel <= spectrum['Max_Velocity']:
                print("WARNING: Ignoring gas component outside of spectrum. Min: {} Max: {} Component: {}".format(
                    spectrum['Min_Velocity'], spectrum['Max_Velocity'], comp_vel))
                continue

        all_gas.append(gas)

        if t_off > 0 and optical_depth < 0.9:
            # Calculate spin temperature and column density
            #denominator = math.min(1-comp_amp)
            t_s = t_off / (1-np.exp(-gas.tau)) # (tb / (1-e^-tau)
            delta_t_s = t_s * ((t_off_err ** 2 / t_off ** 2) +
                               (gas.continuum_sd ** 2 / (1 - np.exp(-gas.tau)) ** 2))**0.5

            # Record
            gas.t_off = t_off
            gas.delta_t_off = t_off_err
            gas.t_s = t_s
            gas.delta_t_s = delta_t_s
            # component['t_s '] = t_s
            print("Component %s at velocity %.4f has t_s %.3f (%.3f/%.3f)" % (
                gas.name, comp_vel, t_s, t_off, comp_amp))
    return all_gas


def is_gas_near_maser(gas):
    if not hasattr(gas, 'maser_loc'):
        return False
    if gas.loc.separation(gas.maser_loc).value > (2 / 60):
        return False
    return gas.maser_vel_low - 10 <= gas.comp_vel <= gas.maser_vel_high + 10


def set_field_metadata(field, ucd, unit, description):
    if ucd:
        field.ucd = ucd
    if unit:
        field.unit = unit
    if description:
        field.description = description


def output_gas_catalogue(all_gas):
    num_gas = len(all_gas)
    names = []
    days = []
    field_names = []
    sources = []
    longitudes = []
    latitudes = []
    ras = []
    decs = []
    velocities = np.zeros(num_gas)
    em_velocities = np.ma.array(np.zeros(num_gas))
    optical_depths = np.zeros(num_gas)
    comp_widths = np.zeros(num_gas)
    temps_off = np.ma.array(np.zeros(num_gas))
    delta_temps_off = np.ma.array(np.zeros(num_gas))
    temps_spin = np.ma.array(np.zeros(num_gas))
    delta_temps_spin = np.ma.array(np.zeros(num_gas))
    temps_kmax = np.zeros(num_gas)
    tau = np.zeros(num_gas)
    continuum_sd = np.zeros(num_gas)
    maser_region = np.empty(num_gas, dtype=bool)
    ratings = np.empty(num_gas, dtype=object)
    filenames = np.empty(num_gas, dtype=object)
    local_paths = np.empty(num_gas, dtype=object)
    local_emission_paths = np.empty(num_gas, dtype=object)
    local_spectra_paths = np.empty(num_gas, dtype=object)

    base_path = os.path.realpath('.')

    for i in range(len(all_gas)):
        gas = all_gas[i]
        names.append(gas.name)
        days.append(gas.day)
        field_names.append(gas.field)
        sources.append(gas.src)
        longitudes.append(gas.longitude)
        latitudes.append(gas.latitude)
        ras.append(gas.ra)
        decs.append(gas.dec)
        velocities[i] = gas.comp_vel
        em_velocities[i] = gas.em_vel / 1000 if gas.em_vel else np.ma.masked
        optical_depths[i] = gas.optical_depth
        comp_widths[i] = gas.comp_width
        temps_kmax[i] = gas.t_kmax
        if gas.t_off is None:
            temps_off[i] = np.ma.masked
            delta_temps_off[i] = np.ma.masked
        else:
            temps_off[i] = gas.t_off
            delta_temps_off[i] = gas.delta_t_off
        if gas.t_s is None:
            temps_spin[i] = np.ma.masked
            delta_temps_spin[i] = np.ma.masked
        else:
            temps_spin[i] = gas.t_s
            delta_temps_spin[i] = gas.delta_t_s
        tau[i] = gas.tau
        maser_region[i] = is_gas_near_maser(gas)
        ratings[i] = gas.rating
        continuum_sd[i] = gas.continuum_sd
        # Need to read in spectra to get rating and include it in the catalogue and
        # link to the fit preview: e.g. plots/A/012.909-0.260_19_src4-1_fit
        prefix = 'day' + str(gas.day) + '/' + gas.field + \
                 "_src" + gas.src
        filenames[i] = prefix + "_plot.png"
        em_filename = prefix + "_emission.png"
        spectra_path = 'run2/plots/{}/{}_{}_src{}_fit.png'.format(gas.rating, gas.field, gas.day, gas.src)
        local_paths[i] = base_path + '/' + filenames[i]
        local_emission_paths[i] = base_path + '/' + em_filename
        local_spectra_paths[i] = base_path + '/' + spectra_path

    # bulk calc fields
    vel_diff = np.abs(velocities - em_velocities)
    equiv_width = np.abs((1-optical_depths) * comp_widths)
    fwhm = np.abs(comp_widths)
    column_density = tau * fwhm * temps_spin * 1.823E18 * 1.064
    sigma = fwhm / (2 * math.sqrt(2 * math.log(2)))
    mach_num = np.sqrt(4.2*((temps_kmax/temps_spin)-1))
    n_wnm = ((column_density * temps_spin) / 50 ) - column_density

    temp_table = Table(
        [names, days, field_names, sources, velocities, em_velocities, optical_depths, temps_off, temps_spin, temps_kmax,
         longitudes, latitudes, ras, decs, fwhm, sigma, vel_diff, equiv_width, tau, maser_region, column_density,
         mach_num, n_wnm, ratings, delta_temps_spin, delta_temps_off, continuum_sd,
         filenames, local_paths, local_emission_paths, local_spectra_paths],
        names=['Comp_Name', 'Day', 'Field', 'Source', 'Velocity', 'em_velocity', 'optical_depth', 'temp_off',
               'temp_spin', 'temp_kmax', 'longitude', 'latitude', 'ra', 'dec', 'fwhm', 'sigma', 'vel_diff',
               'equiv_width', 'tau', 'near_maser', 'column_density', 'mach', 'n_wnm', 'Rating',
               'delta_temp_spin', 'delta_temp_off', 'delta_optical_depth',
               'Filename', 'Local_Path', 'Local_Emission_Path', 'Local_Spectrum_Path'],
        meta={'ID': 'magmo_gas',
              'name': 'MAGMO Gas ' + str(datetime.date.today())})
    votable = from_table(temp_table)
    table = votable.get_first_table()
    set_field_metadata(table.get_field_by_id('longitude'), 'pos.galactic.lon', 'deg',
                       'Galactic longitude of the background source')
    set_field_metadata(table.get_field_by_id('latitude'), 'pos.galactic.lat', 'deg',
                       'Galactic latitude of the background source')
    set_field_metadata(table.get_field_by_id('ra'), 'pos.eq.ra;meta.main', 'deg',
                       'Right ascension of the background source (J2000)')
    set_field_metadata(table.get_field_by_id('dec'), 'pos.eq.dec;meta.main', 'deg',
                       'Declination of the background source (J2000)')
    set_field_metadata(table.get_field_by_id('fwhm'), '', 'km/s', 'Full width at half maximum of the Gaussian component')
    set_field_metadata(table.get_field_by_id('sigma'), '', 'km/s', 'Sigma value of the Gaussian component')
    set_field_metadata(table.get_field_by_id('temp_off'), 'phys.temperature;stat.mean', 'K',
                       'The mean temperature for the gas immediately adjacent to the source')
    set_field_metadata(table.get_field_by_id('temp_spin'), 'phys.temperature;stat.mean', 'K',
                       'The excitation or spin temperature of the gas')
    set_field_metadata(table.get_field_by_id('optical_depth'), '', '',
                       'The peak optical depth of the component ($e^(-\\tau)$)')
    set_field_metadata(table.get_field_by_id('column_density'), '', 'cm-2',
                       'The density of Cold HI gas in the Gaussian component')
    set_field_metadata(table.get_field_by_id('n_wnm'), '', 'cm-2',
                       'The density of Warm HI gas in the Gaussian component')
    set_field_metadata(table.get_field_by_id('mach'), '', '',
                       'The turbulent mach number of the gas in the Gaussian component')
    filename = "magmo-gas.vot"
    writeto(votable, filename)
    return table


def plot_equiv_width_lv(values, ax, key, label, max_clip=None):
    scale_vals = values[key]
    if max_clip:
        scale_vals = np.clip(scale_vals, 0, max_clip)
    base = np.array([values['longitude'], values['Velocity'], scale_vals],
                    dtype=[('l', float), ('b', float), ('scale', float)])

    sample = base # np.sort(base, axis=0, order='ew')
    #vmax = np.max(sample[2,:])
    #cm = plt.cm.get_cmap('RdYlBu_r')
    cm = plt.cm.get_cmap('viridis')
    x = sample[0,:]
    y = sample[1,:]
    z = sample[2,:]
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    # equiv width plot needs log, min 3
    sc = ax.scatter(x, y, c=z, s=25, cmap=cm, norm=matplotlib.colors.LogNorm(vmin=3))
    #sc = ax.scatter(sample[0,:], sample[1,:], c=sample[2,:], s=25, cmap=cm, norm=matplotlib.colors.PowerNorm(0.5, vmin=3))
    #cb = plt.colorbar(sc, norm=matplotlib.colors.LogNorm())
    formatter = LogFormatter(10, labelOnlyBase=False)
    ticks = [0, 5, 10, 20, 30, 40, 60, 80]
    cb = plt.colorbar(sc, ticks=ticks, format=formatter)

    #cb.set_ticklabels([0, 10, 20, 30, 40 , 50 ,60 , 70, 80])
    #cb.ax.get_yaxis().set_ticks([])
    #for j in range(0, 80, 10):
    #    cb.ax.text(.5, (2 * j + 1) / 8.0, str(j), ha='center', va='center')

    ax.set_xlim(values['longitude'].max() + 5, values['longitude'].min() - 5)

    ax.set_xlabel('Galactic longitude (deg)')
    ax.set_ylabel('LSR Velocity (km/s)')
    cb.set_label(label)

    return None


def plot_spin_temp_lv(values, ax, key, label, max_clip=None):
    scale_vals = values[key]
    dataset = values[~ma.getmask(scale_vals)]
    if max_clip:
        scale_vals = np.clip(dataset[key], 0, max_clip)
    base = np.array([dataset['longitude'], dataset['Velocity'], scale_vals],
                    dtype=[('l', float), ('b', float), ('scale', float)])

    sample = base
    #vmax = np.max(sample[2,:])
    #cm = plt.cm.get_cmap('RdYlBu_r')
    cm = plt.cm.get_cmap('viridis')

    x = sample[0,:]
    y = sample[1,:]
    z = sample[2,:]
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    # spin temp plot needs log, min 3
    sc = ax.scatter(x, y, c=z, s=25, cmap=cm, norm=matplotlib.colors.PowerNorm(0.5, vmin=3))
    formatter = LogFormatter(10, labelOnlyBase=False)
    ticks = [0, 5, 10, 20, 40, 70, 100, 150, 200, 300]
    cb = plt.colorbar(sc, ticks=ticks, format=formatter)

    ax.set_xlim(values['longitude'].max() + 5, values['longitude'].min() - 5)

    ax.set_xlabel('Galactic longitude (deg)')
    ax.set_ylabel('LSR Velocity (km/s)')
    cb.set_label(label)

    return None


def plot_spectra(spectra):
    print("## Plotting spectra stats ##")

    fig = plt.figure(figsize=(4, 3))
    gs = matplotlib.gridspec.GridSpec(2, 2)
    labels = ('A', 'B', 'C', 'D')
    ratings = 'ABCD'
    ind = np.arange(len(ratings))

    # Baseline noise box chart
    ax1 = fig.add_subplot(gs[0, 0])
    noise = []
    for rating in ratings:
        noise.append(spectra['Continuum_SD'][spectra['Rating'] == rating].compressed())
    ax1.boxplot(noise, labels=labels, sym='k.')
    ax1.axhline(y=1 / 3, linestyle='--', color='darkgray')
    ax1.semilogy()
    ax1.set_ylabel('$\sigma_{continuum}$')
    ax1.set_title('Baseline Noise', fontsize=10)

    ax2 = fig.add_subplot(gs[0, 1])
    maxsn = []
    for rating in ratings:
        maxsn.append(spectra['Max_S_Max_N'][spectra['Rating'] == rating].compressed())
    ax2.boxplot(maxsn, labels=labels, sym='k.')
    ax2.axhline(y=3, linestyle='--', color='darkgray')
    ax2.semilogy()
    ax2.set_title('Max Signal:Max Noise', fontsize=10)

    # Quality histogram
    ax3 = fig.add_subplot(gs[1, :])
    counts = []
    for rating in ratings:
        counts.append(len(spectra['Rating'][spectra['Rating'] == rating]))
    ax3.bar(ind, counts, width=1, align='center')
    ax3.set_xticks(ind)
    ax3.set_xticklabels(labels)
    ax3.set_xlim(-0.5, len(ind) - 0.5)
    ax3.set_xlabel('Quality Rating')
    ax3.set_ylabel('Number')

    gs.update(wspace=0.5, hspace=0.5)
    filename = 'magmo-spectra-stats.pdf'
    plt.savefig(filename, bbox_inches="tight")
    plt.close()


def plot_maser_comparison(gas_table):
    print("## Comparing gas near and away from masers ##")

    #df = pd.DataFrame(np.ma.filled(gas_array))
    gas_array = gas_table.array
    #df = gas_table.to_table().to_pandas()
    #grid = sns.FacetGrid(df, col="near_maser", margin_titles=True, sharey=False)
    #bins = np.logspace(0.01, 3, 21)
    #grid.map(plt.hist, "temp_spin", color="steelblue", lw=0, bins=bins)
    #plt.savefig('near_maser_temp.pdf')
    #plt.close()

    ts = gas_array['temp_spin']
    indexes = [(0 < ts)]
    ts_gas = gas_array[indexes]
    near_maser = ts_gas[ts_gas['near_maser']]
    away_maser = ts_gas[ts_gas['near_maser'] == False]
    #print(near_maser[0:10], away_maser[0:10])

    fig = plt.figure(figsize=(7.5, 3))
    gs = matplotlib.gridspec.GridSpec(1, 2)

    # Spin Temperature
    ax1 = fig.add_subplot(gs[0, 0])
    away_sample = away_maser['temp_spin']
    bins = np.linspace(0, 450, 19)
    hist, edges = build_hist_fraction(away_sample, bins, 450)
    ax1.step(edges, hist)

    near_sample = near_maser['temp_spin']
    hist, edges = build_hist_fraction(near_sample, bins, 450)
    ax1.step(edges, hist, color='black', ls='--')

    ax1.set_xlabel('Spin Temperature (K)')
    ax1.set_ylabel('Fraction of components')

    # Optical Depth
    ax2 = fig.add_subplot(gs[0, 1])
    away_sample = away_maser['tau']
    bins = np.linspace(0, 4, 21)
    hist, edges = build_hist_fraction(away_sample, bins, 4)
    ax2.step(edges, hist)

    near_sample = near_maser['tau']
    hist, edges = build_hist_fraction(near_sample, bins, 4)
    ax2.step(edges, hist, color='black', ls='--')

    ax2.set_xlabel('Optical depth $(\\tau)$')
    ax2.set_ylabel('Fraction of components')

    gs.update(wspace=0.5, hspace=0.5)
    plt.savefig('near_maser_temp.pdf', bbox_inches="tight")
    plt.close()


    print("Median near:{} v away:{} , sd n:{} v a:{}, count n:{} v a:{}".format(np.median(near_maser['temp_spin']),
                                                             np.median(away_maser['temp_spin']),
                                                             np.std(near_maser['temp_spin']),
                                                             np.std(away_maser['temp_spin']),
                                                             len(near_maser['temp_spin']),
                                                             len(away_maser['temp_spin'])))
    statistic, p_value = stats.ks_2samp(np.ma.filled(near_maser['temp_spin']), np.ma.filled(away_maser['temp_spin']))
    print ('Population T_S similarity K-s={} p_value={}'.format(statistic, p_value))
    statistic, p_value = stats.ks_2samp(np.ma.filled(near_maser['tau']), np.ma.filled(away_maser['tau']))
    print ('Population tau similarity K-s={} p_value={}'.format(statistic, p_value))


def plot_single_hist(gas_array, field, min, max, label=None, sample=None):
    if sample is None:
        sample = gas_array[field]
    hist, edges = np.histogram(sample, bins='auto', range=(min, max))
    outliers = len(sample[sample > max])
    hist[-1] += outliers
    #print ("tau edges",edges)
    plt.bar(edges[:-1], hist, width=edges[1])
    plt.xlabel(label if label else field)
    plt.ylabel('Number of components')
    filename = 'magmo-hist-{}.pdf'.format(field)
    plt.savefig(filename)
    plt.close()
    sample_no_outliers = sample[sample <= max]
    print("{} had {} values, mean {:0.3f}, median {:0.3f}, sd {:0.3f} outliers {} > max {}".format(field, len(sample),
                                                                                          np.mean(sample_no_outliers),
                                                                                          np.median(sample_no_outliers),
                                                                                          np.std(sample_no_outliers),
                                                                                          outliers, max))


def plot_observed(gas_array):
    fig = plt.figure(figsize=(7.5, 3))
    gs = matplotlib.gridspec.GridSpec(1, 2)

    # Optical depth
    ax1 = fig.add_subplot(gs[0, 0])
    sample = gas_array['tau']
    bins = np.linspace(0,4,21)
    hist, edges = np.histogram(sample, bins=bins, range=(0, 4))
    outliers = len(sample[sample > 4])
    hist[-1] += outliers
    hist = np.append(hist, hist[-1])
    #print ("tau edges",edges)
    ax1.step(edges, hist, where='post') # , width=edges[1])
    ax1.set_xlabel('Optical depth $(\\tau)$')
    ax1.set_ylabel('Number of components')

    # FWHM
    ax2 = fig.add_subplot(gs[0, 1])
    sample = gas_array['fwhm']
    bins = np.linspace(0,60,25)
    hist, edges = np.histogram(sample, bins=bins, range=(0, 60))
    outliers = len(sample[sample > 60])
    hist[-1] += outliers
    hist = np.append(hist[0], hist)
    ax2.step(edges, hist) # , width=edges[1])
    ax2.set_xlabel('FWHM (km/s)')
    ax2.set_ylabel('Number of components')

    gs.update(wspace=0.5, hspace=0.5)
    filename = 'magmo-hist-observed.pdf'
    plt.savefig(filename, bbox_inches="tight")
    plt.close()


def plot_equiv_width(gas_array):
    fig = plt.figure(figsize=(7.5, 3))
    gs = matplotlib.gridspec.GridSpec(1, 2)

    # Histogram - all
    ax1 = fig.add_subplot(gs[0, 0])
    sample = gas_array['equiv_width']
    bins = np.linspace(0,40,41)
    hist, edges = np.histogram(sample, bins=bins)
    outliers = len(sample[sample > 40])
    hist[-1] += outliers
    hist = np.append(hist[0], hist)
    ax1.step(edges, hist) # , width=edges[1])
    ax1.set_xlabel('Equivalent Width (km/s)')
    ax1.set_ylabel('Number of components')

    # Histogram - FWHM < 50
    ax2 = fig.add_subplot(gs[0, 1])
    sample = sample[gas_array['fwhm'] < 50]
    hist, edges = np.histogram(sample, bins=bins)
    outliers = len(sample[sample > 40])
    hist[-1] += outliers
    hist = np.append(hist[0], hist)
    ax2.step(edges, hist) # , width=edges[1])
    ax2.set_xlabel('Equivalent Width (km/s)')
    ax2.set_ylabel('Number of components')

    gs.update(wspace=0.5, hspace=0.5)
    filename = 'magmo-equiv-width.pdf'
    plt.savefig(filename, bbox_inches="tight")
    plt.close()


def plot_derived(gas_array):
    fig = plt.figure(figsize=(7.5, 3))
    gs = matplotlib.gridspec.GridSpec(1, 2)

    # Spin Temperature
    ax1 = fig.add_subplot(gs[0, 0])
    sample = np.ma.array(gas_array['temp_spin']).compressed()
    bins = np.linspace(0,450,19)
    hist, edges = np.histogram(sample, bins=bins)
    outliers = len(sample[sample > 450])
    hist[-1] += outliers
    hist = np.append(hist[0], hist)
    ax1.step(edges, hist) # , width=edges[1])
    ax1.set_xlabel('Spin Temperature (K)')
    ax1.set_ylabel('Number of components')

    # Column Density
    ax2 = fig.add_subplot(gs[0, 1])
    sample = np.log10(np.array(gas_array['column_density']))
    bins = np.linspace(19, 24, 21)
    hist, edges = np.histogram(sample, bins=bins)
    outliers = len(sample[sample > 24])
    hist[-1] += outliers
    hist = np.append(hist[0], hist)

    ## Call out the saturated values
    sat_index = gas_array['tau'] > 3
    sat_sample = np.log10(np.array(gas_array['column_density'][sat_index]))
    sat_hist, sat_edges = np.histogram(sat_sample, bins=bins)
    outliers = len(sat_sample[sat_sample > 24])
    sat_hist[-1] += outliers
    sat_hist = np.append(sat_hist[0], sat_hist)

    ax2.step(edges, hist) #, width=edges[1]-edges[0])
    ax2.step(edges, sat_hist, color='r') # , width=edges[1]-edges[0]
    label = 'Column Density $\\log_{10}(N_{H}$) (cm$^{-2}$)'
    ax2.set_xlabel(label)
    ax2.set_ylabel('Number of components')

    gs.update(wspace=0.5, hspace=0.5)
    filename = 'magmo-hist-derived.pdf'
    plt.savefig(filename, bbox_inches="tight")
    plt.close()


def plot_derived2(gas_array):
    fig = plt.figure(figsize=(7.5, 3))
    gs = matplotlib.gridspec.GridSpec(1, 2)

    # Mach number
    ax1 = fig.add_subplot(gs[0, 0])
    sample = np.ma.array(gas_array['mach']).compressed()
    bins = np.linspace(0,50,21)
    hist, edges = np.histogram(sample, bins=bins)
    outliers = len(sample[sample > 50])
    hist[-1] += outliers
    hist = np.append(hist[0], hist)
    ax1.step(edges, hist) # , width=edges[1])
    ax1.set_xlabel('Turbulent Mach number')
    ax1.set_ylabel('Number of components')

    # Cold fraction
    ax2 = fig.add_subplot(gs[0, 1])
    sample = 48 / np.ma.array(gas_array['temp_spin']).compressed()
    bins = np.linspace(0, 1, 21)
    hist, edges = np.histogram(sample, bins=bins)
    outliers = len(sample[sample > 1])
    hist[-1] += outliers
    hist = np.append(hist[0], hist)
    ax2.step(edges, hist) # , width=edges[1]-edges[0])
    label = 'Fraction of cold gas'
    ax2.set_xlabel(label)
    ax2.set_ylabel('Number of components')

    gs.update(wspace=0.5, hspace=0.5)
    filename = 'magmo-hist-derived2.pdf'
    plt.savefig(filename, bbox_inches="tight")
    plt.close()


def plot_histograms(gas_table):
    gas_array = gas_table.array

    plot_observed(gas_array)
    plot_equiv_width(gas_array)
    plot_derived(gas_array)
    plot_derived2(gas_array)

    fig = plt.figure(figsize=(7.5, 3))
    ax = plt.gca()
    plot_equiv_width_lv(gas_array, ax, 'equiv_width', 'Equivalent Width (km/s)')
    plt.savefig('magmo-equiv-width-lv.pdf', bbox_inches="tight")
    plt.close()

    fig = plt.figure(figsize=(7.5, 3))
    ax = plt.gca()
    plot_spin_temp_lv(gas_array, ax, 'temp_spin', 'Spin Temperature (K)', max_clip=300)
    plt.savefig('magmo-spin-temp-lv.pdf', bbox_inches="tight")
    plt.close()

    plot_single_hist(gas_array, 'optical_depth', 0, 1.2, label='$e^{-\\tau}$')
    plot_single_hist(gas_array, 'tau', 0, 8, label='$\\tau$')
    plot_single_hist(gas_array, 'fwhm', 0, 60, label='FWHM (km/s)')
    #plot_single_hist(gas_array, 'mach', 0, 50, label='Mach Number (M$_t$)')
    plot_single_hist(gas_array, 'temp_spin', 0, 450, label='Spin Temperature (K)', sample = np.ma.array(gas_array['temp_spin']).compressed())

    # Column density is special
    sample = np.log10(np.array(gas_array['column_density'])) # / 20.0
    bins = np.linspace(19, 24, 21)
    #bins = 10 ** np.linspace(np.log10(np.min(sample)), np.log10(np.max(sample)), 50)
    hist, edges = np.histogram(sample, bins=bins)
    outliers = len(sample[sample > 24])
    hist[-1] += outliers

    plt.bar(edges[:-1], hist, width=edges[1]-edges[0])
    label = 'Column Density $\\log_{10}(N_{H}$) (g/cm^3)'
    plt.xlabel(label)
    plt.ylabel('Number of components')
    filename = 'magmo-hist-column_density.pdf'
    plt.savefig(filename, bbox_inches="tight")
    plt.close()
    #sample_no_outliers = sample[sample <= max]
    print("{} had {} values, mean {:0.3f}, median {:0.3f}, sd {:0.3f} outliers".format('column_density', len(sample),
                                                                                          np.mean(sample),
                                                                                          np.median(sample),
                                                                                          np.std(sample))

    #plot_single_hist(gas_array, 'column_density', 0, 1e34, label='Column Density $N_{H20}$ (1E20 g/cm^3)',
                     )
    #tau = gas_array['tau']
    #hist, edges = np.histogram(tau, bins='fd', range=(0,8))
    #outliers = len(tau[tau > 8])
    #hist[-1] += outliers
    #print ("tau edges",edges)
    #plt.bar(edges[:-1], hist, width=edges[1]) # align='edge',
    #plt.xlabel('$\\tau$')
    #plt.ylabel('Number of components')
    #filename = 'magmo-tau.pdf'
    #plt.savefig(filename)
    #plt.close()
    #print("tau had {} values, mean {:0.3f}, median {:0.3f}, sd {:0.3f} outliers {}".format(len(tau), np.mean(tau),
    #                                                                                       np.median(tau),
    #                                                                                       np.std(tau), outliers))


def plot_uncertainty(gas_table):
    fig = plt.figure(figsize=(7.5, 3))
    gs = matplotlib.gridspec.GridSpec(1, 2)

    # Uncertainty spread
    gas_array = gas_table.array
    x = np.ma.array(gas_array['temp_spin']).compressed()
    y = np.ma.array(gas_array['delta_temp_spin']).compressed()

    # Calculate the point density
    xy = np.vstack([x, y])
    z = stats.gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    ax1 = fig.add_subplot(gs[0, 0])
    ax1.scatter(x, y, c=z, s=20, edgecolor='', cmap='viridis')
    ax1.set_ylim((0, 400))
    ax1.set_xlim(0)
    ax1.set_xlabel("Spin Temperature (K)")
    ax1.set_ylabel("$\Delta$ Spin Temperature (K)")

    # Uncertainty by fraction
    x = np.ma.array(gas_array['temp_spin']).compressed()
    y = np.ma.array(gas_array['delta_temp_spin']).compressed() / x
    print("Median T_S uncertainty fracton is", np.median(y))
    # Calculate the point density
    xy = np.vstack([x, y])
    z = stats.gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    #fig, ax = plt.subplots()
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.scatter(x, y, c=z, s=20, edgecolor='', cmap='viridis')
    ax2.set_ylim((0,2))
    ax2.set_xlim(0)
    # plt.yscale('log')
    ax2.set_xlabel("Spin Temperature (K)")
    ax2.set_ylabel("$\Delta T_S / T_S$")
    #plt.savefig('magmo-spin-temp-uncertainty-fraction.pdf', bbox_inches="tight")
    #plt.close()
    plt.savefig('magmo-spin-temp-uncertainty.pdf', bbox_inches="tight")
    plt.close()


def find_best_matches(magmo_coords, other_coords, max_dist, magmo_table):
    """
    Match the best MAGMO spectra within a defined distance to each source in a list of locations.
    Best here is best quality and then the lowest continuum standard deviation.

    :param magmo_coords: An array of SkyCoord objects representing the location of each MAGMO spectrum.
    :param other_coords: An array of SkyCoord objects representing the location of the sources to be matched.
    :param max_dist: The maximum allowed separation for matches, a quantity with angle units.
    :param magmo_table: The table of MAGMO spectra details in the same order as the magmo_coords.
    :return: An array of MAGMO index, other index and distance (in decimal degrees) for each match.
    """

    idx_other, idx_magmo, d2d, d3d = magmo_coords.search_around_sky(other_coords, max_dist)

    match_map = dict()
    dist_map = dict()
    cat2_match = []
    magmo_cat2_dist = []
    for val in np.unique(idx_other):
        match_map[val] = []
        dist_map[val] = []

    for i in range(len(idx_other)):
        match_map[idx_other[i]].append(idx_magmo[i])
        dist_map[idx_other[i]].append(d2d[i])

    for key in match_map.keys():
        magmo_rows = magmo_table[match_map[key]]
        distances = dist_map[key]
        prev_idx = -1
        prev_row = None
        dist = 0
        for i in range(len(magmo_rows)):
            if prev_idx < 0 or magmo_rows[i]['Rating'] < prev_row['Rating'] \
                    or magmo_rows[i]['Continuum_SD'] < prev_row['Continuum_SD']:
                prev_idx = i
                prev_row = magmo_rows[i]
                dist = distances[i].to(u.degree).value
        cat2_match.append(match_map[key][prev_idx])
        magmo_cat2_dist.append(dist)
    magmo_match = match_map.keys()
    matches = np.stack((magmo_match, cat2_match, magmo_cat2_dist), axis=-1)
    return matches


def get_opacity_filename(day, field, src_id):
    name_prefix = field + '_src' + src_id
    filename = 'day' + str(day) + '/' + name_prefix + '_opacity.votable.xml'
    return filename


def get_brown_filename(lat, lon):
    sep = '+' if lon >= 0 else ''
    filename = '../brown_data/{:.3f}{}{:.3f}.dat'.format(lat, sep, lon)
    return filename


def get_brown_graph_filename(lat, lon):
    sep = '+' if lon >= 0 else ''
    name = 'sgps_comp_{:.3f}{}{:.3f}'.format(lat, sep, lon)
    name = name.replace('.', '_')
    filename = 'figures/{}.pdf'.format(name)
    return filename


def resample(values, old_grid, other_grid, fill_value=1.0):
    """ resample values from the old_grid to match the new_grid, with linear interpolation """

    old_step = old_grid[1] - old_grid[0]
    new_step = other_grid[1] - other_grid[0]
    increment = new_step / old_step
    start = 0
    if old_grid[0] < other_grid[0]:
        start = np.searchsorted(old_grid, other_grid[0]) - 1
    num_steps = min(len(old_grid) / increment, len(other_grid)) - start - 1

    # Then resample data to match new_grid
    start = 0
    if old_grid[0] > other_grid[0]:
        start = np.searchsorted(other_grid, old_grid[0]) - 1

    f = interpolate.interp1d(old_grid, values, fill_value=fill_value, bounds_error=False)
    resampled = f(other_grid[start:int(start + num_steps)])
    padded = np.full(other_grid.shape, fill_value)
    padded[start:start + len(resampled)] = resampled
    return padded


def plot_spectrum(x, y, ax, title, line=True, **kwargs):
    """
    Output a plot of opacity vs LSR velocity to a specified file.

    :param x: The velocity data
    :param y: The opacity values for each velocity step
    :param filename: The file the plot should be written to. Should be
         an .eps or .pdf file.
    :param title: The title for the plot
    """
    fig = plt.figure()
    ax.plot(x / 1000, y, **kwargs)

    if line:
        ax.axhline(1, color='r')

    ax.set_xlabel(r'Velocity relative to LSR (km/s)', fontsize=11)
    ax.set_ylabel(r'$e^{(-\tau)}$')
    ax.set_title(title, fontsize=11)
    ax.set_xlim([-250, 250])
    ax.grid(True)
    return


def plot_source(title1, title2, brown_data, resampled_spec, resampled_sigma, filename):
    """
    Produce a plot of the two spectra (left) and the residual (right) and output it to a file.
    :param title1:  The first line of the title of the left panel.
    :param title2: The second line of the title of the left panel.
    :param brown_data: The spectrum data from the Brown et al 2014 table 2 for the source being plotted
    :param resampled_spec: The MAGMO spectrum data for the source at velocity steps matching the Brown et al data.
    :param cont_sd: The measured optical depth noise level of the continuum section of the MAGMO spectrum
    :param filename: The name of the file to write out the plot to.
    :return: None
    """
    # Plot spectra
    fig = plt.figure(figsize=(8.27, 3.8))
    ax = fig.add_subplot(1, 2, 1)
    title = title1 + '\n' + title2
    velocity = brown_data['col1'] * 1000
    plot_spectrum(velocity, resampled_spec, ax, title, color='mediumturquoise')
    plot_spectrum(velocity, brown_data['col5'], ax, title, line=False, color='blue')

    # Plot Residual
    ax = fig.add_subplot(1, 2, 2)
    plot_spectrum(velocity, resampled_spec - brown_data['col5'], ax,
                  'Difference between MAGMO and SGPS', line=False, color='gray')
    ax.set_ylim([-1.5, 1.5])
    plot_spectrum(velocity, (resampled_sigma * 3), ax, '', line=False, color='lightsteelblue', linewidth=1)
    plot_spectrum(velocity, -(resampled_sigma * 3), ax, '', line=False, color='lightsteelblue', linewidth=1)
    fig.tight_layout()
    fig.savefig(filename, bbox_inches="tight")
    plt.close()
    return


def assess_match(brown_data, resampled_spec, resampled_sigma):
    noise_threshold = 3*resampled_sigma
    residual = resampled_spec - brown_data['col5']
    abs_residual = np.abs(residual)
    noisy = residual[abs_residual > noise_threshold]
    #print ("Found {} out {} outside noise threshold.".format(len(noisy), len(residual)))
    return len(noisy)


def compare_brown_2014(magmo_coords, magmo_table):
    """
    Compare the MAGMO spectra with the Brown et al 2014 SGPS spectra for HII regions.
    This will produce a catalogue of sources in both datasets and comparative plots
    showing the two spectra and the difference for each source.
    :return: The number of matched spectra
    """
    print("## Comparing with Brown et al 2014 ##")

    magmo.ensure_dir_exists('figures')

    brown_table = ascii.read('../Brown-2014-HII-Regions.dat', format='cds')
    brown_coords = SkyCoord(brown_table['GLON'], brown_table['GLAT'], frame='galactic', unit="deg")

    matches = find_best_matches(magmo_coords, brown_coords, 2.16 * u.arcmin, magmo_table)
    t1 = brown_table[matches[:, 0].astype(int)]
    t2 = magmo_table[matches[:, 1].astype(int)]
    combined = hstack([t1, t2], join_type='exact')
    dist_col = Column(name='Separation', data=matches[:, 2], unit=u.degree)
    combined.add_column(dist_col)
    combined.sort('GLON')
    noisy_count = np.zeros(len(combined))
    noisy_col = Column(name='Noisy_Count', data=noisy_count,
                      description='Count of optical depth measurements more than 4 sigma difference between MAGMO and SGPS')
    combined.add_column(noisy_col)

    # Produce plots and a latex template
    i = 1
    with open('brown_2014_comp.tex', 'w') as latex_doc:
        for match_row in combined:
            brown_filename = get_brown_filename(match_row['GLON'], match_row['GLAT'])
            if not os.path.exists(brown_filename):
                print("Unable to find {}, skipping".format(brown_filename))
                continue

            brown_data = ascii.read(brown_filename)
            spectrum = read_votable_results(
                get_opacity_filename(match_row['Day'], match_row['Field'], match_row['Source']))
            resampled_spec = resample(spectrum['opacity'], spectrum['velocity'], brown_data['col1'] * 1000)

            cont_sd = match_row['Continuum_SD']
            resampled_sigma = resample(spectrum['sigma_tau'], spectrum['velocity'], brown_data['col1'] * 1000,
                                       fill_value=cont_sd)

            match_row['Noisy_Count'] = assess_match(brown_data, resampled_spec, resampled_sigma)

            title1 = match_row['SIMBAD']
            #title2 = 'Day {} Field {} Src {}'.format(match_row['Day'], match_row['Field'], match_row['Source'])
            title2 = match_row['Name']
            filename = get_brown_graph_filename(match_row['GLON'], match_row['GLAT'])
            plot_source(title1, title2, brown_data, resampled_spec, resampled_sigma, filename)
            latex_doc.write('\n\\begin{figure}[ht]\n')
            latex_doc.write('\\includegraphics[scale=0.75]{' + filename + '}\n')
            latex_doc.write(
                '\\caption{{Comparison of MAGMO and SGPS absorption spectra for {} }}\n'.format(match_row['SIMBAD']))
            latex_doc.write('\\label{{fig:Brown-2014-{} }}\n'.format(i))
            latex_doc.write('\\end{figure}\n')
            i += 1

    brown_table.write('sgps-hii.csv', format='csv', overwrite=True)
    combined.write('bm-combined.vot', format='votable', overwrite=True)
    #combined.write('bm-combined.csv', format='csv', overwrite=True)
    with open('bm-combined.csv', "wb") as stats:
        writer = csv.writer(stats) #, quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(["SIMBAD", "MAGMO Name", "Separation", "Noise Count", "Qual", "Rating", "sigma_tau", "l", "b"])
        for i in range(0, len(combined)):
            writer.writerow(
                [combined[i]['SIMBAD'], combined[i]['Name'], "{:0.0f}".format(combined[i]['Separation']*3600),
                 "{:0.0f}".format(combined[i]['Noisy_Count']), combined[i]['Qual'], combined[i]['Rating'],
                 "{:0.2f}".format(math.log(1 - combined[i]['Continuum_SD']) * -1),
                 "{:0.4f}".format(combined[i]['GLON']), "{:0.4f}".format(combined[i]['GLAT'])])

    return len(matches)


def get_magmo_table():
    votable = parse("magmo-spectra.vot", pedantic=False)
    table = votable.get_first_table()
    table = table.to_table()
    magmo_table = table[table['Duplicate'] == False]
    magmo_table = magmo_table[magmo_table['Rating'] <= 'C']
    magmo_coords = SkyCoord(magmo_table['Longitude'], magmo_table['Latitude'], frame='galactic', unit="deg")
    return magmo_coords, magmo_table


def compare_dickey_2003(magmo_coords, magmo_table):
    print("## Comparing with Dickey et al 2003 ##")

    dickey_table = ascii.read('../Dickey_2003_sources.txt')
    dickey_coords = SkyCoord(dickey_table['RA'], dickey_table['Dec'], frame='fk5', unit="deg")

    matches = find_best_matches(magmo_coords, dickey_coords, 2.16 * u.arcmin, magmo_table)
    t1 = dickey_table[matches[:, 0].astype(int)]
    t2 = magmo_table[matches[:, 1].astype(int)]
    combined = hstack([t1, t2], join_type='exact')
    dist_col = Column(name='Separation', data=matches[:, 2], unit=u.degree)
    combined.add_column(dist_col)
    combined.sort('Name_1')
    combined.write('dm-combined.vot', format='votable', overwrite=True)


def build_hist_fraction(sample, bins, clip_val):
    hist, edges = np.histogram(sample, bins=bins)
    outliers = len(sample[sample > clip_val])
    hist[-1] += outliers
    hist = hist / np.sum(hist)
    hist = np.append(hist[0], hist)
    return hist, edges


def compare_heiles_2003(magmo_gas):
    print("## Comparing with Heiles & Troland 2003 ##")

    # Read in ht03 data
    rdr = ascii.get_reader(Reader=ascii.Csv)
    ht03_table = rdr.read('../Millennium_data.csv')
    # filter for just the CNM data |b| < 10
    cnm = np.array(ht03_table['CNM'])
    sample = ht03_table[cnm >= '0']
    abs_lat = np.absolute(np.array(sample['GLAT']))
    ht03_low_cnm = sample[abs_lat <= 10]
    spin_temp = np.array(ht03_low_cnm['Ts'])
    print("Sample had {} values, mean {:0.3f}, median {:0.3f}, sd {:0.3f}".format(len(spin_temp),
                                                                                  np.mean(spin_temp),
                                                                                  np.median(spin_temp),
                                                                                  np.std(spin_temp)))

    votable = from_table(ht03_low_cnm)
    writeto(votable, 'millenium_spin.vot')

    # comparative histogram of the two CNM spin temp sets
    fig = plt.figure(figsize=(7.5, 3))
    gs = matplotlib.gridspec.GridSpec(1, 2)
    gas_array = magmo_gas

    # Spin Temperature
    ax1 = fig.add_subplot(gs[0, 0])
    sample = np.ma.array(gas_array['temp_spin']).compressed()
    bins = np.linspace(0,450,19)
    hist, edges = build_hist_fraction(sample, bins, 450)
    ax1.step(edges, hist)

    ht03_sample = np.array(ht03_low_cnm['Ts'])
    hist, edges = build_hist_fraction(ht03_sample, bins, 450)
    ax1.step(edges, hist, color='black', ls='--')

    ax1.set_xlabel('Spin Temperature (K)')
    ax1.set_ylabel('Fraction of components')

    statistic, p_value = stats.ks_2samp(np.ma.filled(sample), np.ma.filled(ht03_sample))
    print ('Spin temp population similarity p_value={}'.format(p_value))


    # Column Density
    ax2 = fig.add_subplot(gs[0, 1])
    sample = np.ma.array(gas_array['column_density']).compressed()
    sample = np.log10(sample)
    bins = np.linspace(19, 24, 21)
    hist, edges = build_hist_fraction(sample, bins, 24)

    sample = np.array(ht03_low_cnm['NHI']) * 1E20
    sample = np.log10(sample[sample > 0])
    ht03_hist, edges = build_hist_fraction(sample, bins, 24)

    ax2.step(edges, hist) #, width=edges[1]-edges[0])
    ax2.step(edges, ht03_hist, color='black', ls=':') # , width=edges[1]-edges[0]
    label = 'Column Density $\\log_{10}(N_{H}$) (cm$^{-2}$)'
    ax2.set_xlabel(label)
    ax2.set_ylabel('Fraction of components')

    statistic, p_value = stats.ks_2samp(np.ma.filled(gas_array['column_density']), np.ma.filled(np.array(ht03_low_cnm['NHI']) * 1E20))
    print ('Column density population similarity p_value={}'.format(p_value))

    gs.update(wspace=0.5, hspace=0.5)
    filename = 'magmo-heiles_2003_comp.pdf'
    plt.savefig(filename, bbox_inches="tight")
    plt.close()

    return


def compare_strasser_2007(magmo_gas):
    print("## Comparing with Strasser et al 2007 ##")

    # Read in s07 data
    s07_files = glob.glob('../strasser2007/*.dat')
    s07_data = None
    for file_name in s07_files:
        data = ascii.read(file_name)
        spin_temp = np.array(data['Ts'])
        subset = data[spin_temp > 0]
        if s07_data is None:
            s07_data = subset
        else:
            s07_data = vstack([s07_data, subset])

    s07_spin_temp = s07_data['Ts']
    s07_spin_temp = s07_spin_temp[s07_spin_temp<1E4]
    print("Sample had {} values, mean {:0.3f}, median {:0.3f}, sd {:0.3f}".format(len(s07_spin_temp),
                                                                                  np.mean(s07_spin_temp),
                                                                                  np.median(s07_spin_temp),
                                                                                  np.std(s07_spin_temp)))
    votable = from_table(s07_data)
    writeto(votable, 'strasser_2007_spin.vot')

    fig = plt.figure(figsize=(7.5, 3))
    gs = matplotlib.gridspec.GridSpec(1, 2)
    gas_array = magmo_gas

    # Spin Temperature
    ax1 = fig.add_subplot(gs[0, 0])
    sample = np.ma.array(gas_array['temp_spin']).compressed()
    bins = np.linspace(0,450,19)
    hist, edges = build_hist_fraction(sample, bins, 450)
    ax1.step(edges, hist)

    s07_sample = np.array(s07_data['Ts'])
    s07_sample = s07_sample[s07_sample<1E4]
    hist, edges = build_hist_fraction(s07_sample, bins, 450)
    ax1.step(edges, hist, color='black', ls='--')

    ax1.set_xlabel('Spin Temperature (K)')
    ax1.set_ylabel('Fraction of components')

    statistic, p_value = stats.ks_2samp(np.ma.filled(sample), np.ma.filled(s07_sample))
    print ('Spin temp population similarity p_value={}'.format(p_value))

    gs.update(wspace=0.5, hspace=0.5)
    filename = 'magmo-strasser_2007_comp.pdf'
    plt.savefig(filename, bbox_inches="tight")
    plt.close()
    return


def report_close_neighbours(magmo_coords, magmo_table):
    idx, sep2d, dist3d = matching.match_coordinates_sky(magmo_coords, magmo_coords, 2)
    idx_orig = np.arange(len(magmo_coords))
    idx_close = sep2d < 3 * u.arcsec
    idx_match1 = idx_orig[idx_close]
    idx_match2 = idx[idx_close]
    sep_12 = sep2d[idx_close]
    for i in range(len(idx_match1)):
        match1 = magmo_table[idx_match1[i]]
        match2 = magmo_table[idx_match2[i]]
        print("{} is only {:.2f} arcsec from {} Rating {} v {} ContSD {:.3f} v {:.3f}".format(match1['Name'],
                                                                                      sep_12[i].to(u.arcsec),
                                                                                      match2['Name'], match1['Rating'],
                                                                                      match2['Rating'],
                                                                                      match1['Continuum_SD'],
                                                                                      match2['Continuum_SD']))


def main():
    # Parse command line options
    args = parseargs()

    start = time.time()
    print("#### Started examining MAGMO components at %s ####" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    components = []
    all_gas = []
    # Load component catalogue
    if not args.compare_only:
        components = read_components("magmo-components.vot")
        spectra_map, spectra_array = read_spectra("magmo-spectra.vot")
        mmb_map = read_mmb_cat('methanol_multibeam_catalogue.vot')

        all_gas = analyse_components(components, spectra_map, mmb_map)

        # Output a catalogue
        gas_table = output_gas_catalogue(all_gas)

        # Plots of MAGMO specific
        plot_maser_comparison(gas_table)
        plot_histograms(gas_table)
        plot_spectra(spectra_array)
        plot_uncertainty(gas_table)

    # Comparisons
    num_matches_brown = 0
    if not args.no_compare:
        magmo_coords, magmo_table = get_magmo_table()
        votable = parse("magmo-gas.vot", pedantic=False)
        magmo_gas = votable.get_first_table().to_table()
        num_matches_brown = compare_brown_2014(magmo_coords, magmo_table)
        compare_dickey_2003(magmo_coords, magmo_table)
        report_close_neighbours(magmo_coords, magmo_table)
        compare_heiles_2003(magmo_gas)
        compare_strasser_2007(magmo_gas)

    # Report
    end = time.time()
    print('#### Examination completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Processed %d components and output %d gas stats, %d matches in %.02f s' %
          (len(components), len(all_gas), num_matches_brown, end - start))
    return 0


if __name__ == '__main__':
    exit(main())

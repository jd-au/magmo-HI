#!/usr/bin/env python -u

# Take the Gaussian components and use them to examine the gas characteristics at each position.
#

# Author James Dempsey
# Date 26 Mar 2017


from __future__ import print_function, division

import argparse
import datetime
import os
import re
from string import Template
import time

from astropy.coordinates import SkyCoord, matching
from astropy.io.votable import parse, from_table, writeto
from astropy.table import Table, Column, hstack
from astropy.io import ascii
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import scipy.stats as stats
import seaborn as sns


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
    temp = 0
    for i in range(0, len(velocities)):
        if velocities[i] >= comp_vel:
            temp = emission['em_mean'][i]
            return temp, velocities[i]
    return 0, 0


def analyse_components(components, spectra_map, mmb_map):
    all_gas = []
    for component in components:
        # Load emission data
        emission_filename = get_emission_filename(component['Day'], component['Field'], component['Source'])
        if os.path.exists(emission_filename):
            emission = read_emission(emission_filename)

            comp_vel = component['Mean']
            comp_width = component['FWHM']
            comp_amp = component['Amplitude']
            t_off, em_vel = get_temp(emission, comp_vel * 1000)
            optical_depth = 1 - comp_amp

            gas = Gas(component['Day'], component['Field'], component['Source'])
            gas.comp_vel = comp_vel
            gas.comp_amp = comp_amp
            gas.comp_width = comp_width
            gas.optical_depth = optical_depth
            gas.em_vel = em_vel
            gas.longitude = component['Longitude']
            gas.latitude = component['Latitude']
            gas.tau = -1 * np.log(np.maximum(optical_depth, 1e-16))
            gas.t_off = None
            gas.t_s = None
            gas.name = component['Comp_Name']
            gas.spectra_name = component['Spectra_Name']

            spectrum = spectra_map[get_spectra_key(component['Day'], component['Field'], component['Source'])]
            gas.rating = spectrum['Rating']

            # Validate the component velocity
            if not spectrum['Min_Velocity'] <= comp_vel <= spectrum['Max_Velocity']:
                print("WARNING: Ignoring gas component outside of spectrum. Min: {} Max: {} Component: {}".format(
                    spectrum['Min_Velocity'], spectrum['Max_Velocity'], comp_vel))
                continue

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

            all_gas.append(gas)

            if t_off > 0 and comp_amp > 0.02:
                # Calculate spin temperature and column density
                t_s = t_off / comp_amp

                # Record
                gas.t_off = t_off
                gas.t_s = t_s
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
    em_velocities = np.zeros(num_gas)
    optical_depths = np.zeros(num_gas)
    comp_widths = np.zeros(num_gas)
    temps_off = np.ma.array(np.zeros(num_gas))
    temps_spin = np.ma.array(np.zeros(num_gas))
    tau = np.zeros(num_gas)
    maser_region = np.empty(num_gas, dtype=bool)
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
        em_velocities[i] = gas.em_vel / 1000
        optical_depths[i] = gas.comp_amp
        comp_widths[i] = gas.comp_width
        if gas.t_off is None:
            temps_off[i] = np.ma.masked
        else:
            temps_off[i] = gas.t_off
        if gas.t_s is None:
            temps_spin[i] = np.ma.masked
        else:
            temps_spin[i] = gas.t_s
        tau[i] = gas.tau
        maser_region[i] = is_gas_near_maser(gas)
        # Need to read in spectra to get rating and include it in the catalogue and
        # link to the fit preview: e.g. plots/A/012.909-0.260_19_src4-1_fit
        prefix = 'day' + str(gas.day) + '/' + gas.field + \
                 "_src" + gas.src
        filenames[i] = prefix + "_plot.png"
        em_filename = prefix + "_emission.png"
        spectra_path = 'plots/{}/{}_{}_src{}_fit.png'.format(gas.rating, gas.field, gas.day, gas.src)
        local_paths[i] = base_path + '/' + filenames[i]
        local_emission_paths[i] = base_path + '/' + em_filename
        local_spectra_paths[i] = base_path + '/' + spectra_path

    # bulk calc fields
    vel_diff = np.abs(velocities - em_velocities)
    equiv_width = np.abs((1 - optical_depths) * comp_widths)
    column_density = equiv_width * temps_spin * 1.8E18

    temp_table = Table(
        [names, days, field_names, sources, velocities, em_velocities, optical_depths, temps_off, temps_spin,
         longitudes, latitudes, ras, decs, comp_widths, vel_diff, equiv_width, tau, maser_region, column_density,
         filenames, local_paths, local_emission_paths, local_spectra_paths],
        names=['Comp_Name', 'Day', 'Field', 'Source', 'Velocity', 'em_velocity', 'Optical_Depth', 'temp_off',
               'temp_spin', 'longitude', 'latitude', 'ra', 'dec', 'fwhm', 'vel_diff', 'equiv_width', 'tau',
               'near_maser', 'Column_density', 'Filename', 'Local_Path', 'Local_Emission_Path', 'Local_Spectrum_Path'],
        meta={'ID': 'magmo_gas',
              'name': 'MAGMO Gas ' + str(datetime.date.today())})
    votable = from_table(temp_table)
    table = votable.get_first_table()
    table.get_field_by_id('ra').ucd = 'pos.eq.ra;meta.main'
    table.get_field_by_id('dec').ucd = 'pos.eq.dec;meta.main'
    filename = "magmo-gas.vot"
    writeto(votable, filename)
    return table


def plot_equiv_width_lv(gas_table):
    values = gas_table.array
    cm = plt.cm.get_cmap('RdYlBu_r')
    sc = plt.scatter(values['longitude'], values['Velocity'], c=values['equiv_width'], s=35, cmap=cm)
    cb = plt.colorbar(sc, norm=matplotlib.colors.LogNorm())

    ax = plt.gca()
    ax.set_xlim(values['longitude'].max() + 5, values['longitude'].min() - 5)

    plt.title("Equivalent Width of Fitted Gas Components")
    plt.xlabel('Galactic longitude (deg)')
    plt.ylabel('LSR Velocity (km/s)')
    cb.set_label('Equivalent Width (km/s)')

    filename = 'magmo-equiv-width-lv.pdf'
    # plt.show()
    plt.savefig(filename)
    return None


def plot_maser_comparison(gas_table):
    print("## Comparing gas near and away from masers ##")

    #df = pd.DataFrame(np.ma.filled(gas_array))
    gas_array = gas_table.array
    df = gas_table.to_table().to_pandas()
    grid = sns.FacetGrid(df, col="near_maser", margin_titles=True, sharey=False)
    bins = np.logspace(0.01, 3.1, 21)
    grid.map(plt.hist, "temp_spin", color="steelblue", lw=0, bins=bins)
    plt.savefig('near_maser_temp.pdf')

    ts = gas_array['temp_spin']
    indexes = [~np.isnan(ts)]
    ts_gas = gas_array[indexes]

    print(ts[0:10], ts_gas[0:10])
    near_maser = ts_gas[ts_gas['near_maser']]
    away_maser = ts_gas[ts_gas['near_maser'] == False]
    print(near_maser[0:10], away_maser[0:10])
    print("Median near:{} v away:{} , sd n:{} v a:{}".format(np.median(near_maser['temp_spin']),
                                                             np.median(away_maser['temp_spin']),
                                                             np.std(near_maser['temp_spin']),
                                                             np.std(away_maser['temp_spin'])))
    statistic, p_value = stats.ks_2samp(np.ma.filled(near_maser['temp_spin']), np.ma.filled(away_maser['temp_spin']))
    print ('Population similarity p_value={}'.format(p_value))


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
    print ("Found {} out {} outside noise threshold.".format(len(noisy), len(residual)))
    return len(noisy)


def compare_brown_2014(magmo_coords, magmo_table):
    """
    Compare the MAGMO spectra with the Brown et al 2014 SGPS spectra for HII regions.
    This will produce a catalogue of sources in both datasets and comparative plots
    showing the two spectra and the difference for each source.
    :return: The number of matched spectra
    """
    print("## Comparing with Brown et al 2014 ##")

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
    dist_col = Column(name='Noisy_Count', data=noisy_count,
                      description='Count of optical depth measurements more than 4 sigma difference between MAGMO and SGPS')
    combined.add_column(dist_col)

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
    combined.write('bm-combined.csv', format='csv', overwrite=True)
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
        plot_equiv_width_lv(gas_table)

        # Plots of MAGMO specific
        plot_maser_comparison(gas_table)

    # Comparisons
    magmo_coords, magmo_table = get_magmo_table()
    num_matches_brown = compare_brown_2014(magmo_coords, magmo_table)
    compare_dickey_2003(magmo_coords, magmo_table)
    report_close_neighbours(magmo_coords, magmo_table)

    # Report
    end = time.time()
    print('#### Examination completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Processed %d components and output %d gas stats, %d matches in %.02f s' %
          (len(components), len(all_gas), num_matches_brown, end - start))
    return 0


if __name__ == '__main__':
    exit(main())

#!/usr/bin/env python -u

# Take the Gaussian components and use them to examine the gas characteristics at each position.
#

# Author James Dempsey
# Date 26 Mar 2017


from __future__ import print_function, division

from astropy.coordinates import SkyCoord
from astropy.io.votable import parse, from_table, writeto
from astropy.table import Table, Column
from numpy import ma
from string import Template

import argparse
import datetime
import numpy as np
import os
import time


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
    return spectra_map


def read_emission(filename):
    return read_votable_results(filename)


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


def analyse_components(components, spectra_map):
    all_gas = []
    for component in components:
        # Load emission data
        emission_filename = get_emission_filename(component['Day'], component['Field'], component['Source'])
        if os.path.exists(emission_filename):
            emission = read_emission(emission_filename)

            comp_vel = component['Mean']
            comp_width = component['FWHM']
            comp_amp = component['Amplitude']
            t_off, em_vel = get_temp(emission, comp_vel*1000)
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

            spectrum = spectra_map[get_spectra_key(component['Day'], component['Field'], component['Source'])]
            gas.rating = spectrum['Rating']

            loc = SkyCoord(gas.longitude, gas.latitude, frame='galactic', unit="deg")
            gas.ra = loc.icrs.ra.degree
            gas.dec = loc.icrs.dec.degree

            all_gas.append(gas)

            if t_off > 0 and comp_amp < 0.98:
                # Calculate spin temperature and column density
                t_s = t_off / comp_amp

                # Record
                gas.t_off = t_off
                gas.t_s = t_s
                #component['t_s '] = t_s
                print ("src %s at velocity %.4f has t_s %.3f (%.3f/%.3f)" % (component['Field'], comp_vel, t_s, t_off, comp_amp))
    return all_gas


def output_gas_catalogue(all_gas):
    num_gas = len(all_gas)
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
    temps_off = ma.array(np.zeros(num_gas))
    temps_spin = ma.array(np.zeros(num_gas))
    tau = np.zeros(num_gas)
    filenames = np.empty(num_gas, dtype=object)
    local_paths = np.empty(num_gas, dtype=object)
    local_emission_paths = np.empty(num_gas, dtype=object)
    local_spectra_paths = np.empty(num_gas, dtype=object)

    base_path = os.path.realpath('.')

    for i in range(len(all_gas)):
        gas = all_gas[i]
        days.append(gas.day)
        field_names.append(gas.field)
        sources.append(gas.src)
        longitudes.append(gas.longitude)
        latitudes.append(gas.latitude)
        ras.append(gas.ra)
        decs.append(gas.dec)
        velocities[i] = gas.comp_vel
        em_velocities[i] = gas.em_vel/1000
        optical_depths[i] = gas.comp_amp
        comp_widths[i] = gas.comp_width
        if gas.t_off is None:
            temps_off[i] = ma.masked
        else:
            temps_off[i] = gas.t_off
        if gas.t_s is None:
            temps_spin[i] = ma.masked
        else:
            temps_spin[i] = gas.t_s
        tau[i] = gas.tau
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
    vel_diff = np.abs(velocities-em_velocities)
    equiv_width = np.abs((1-optical_depths) * comp_widths)

    temp_table = Table(
        [days, field_names, sources, velocities, em_velocities, optical_depths, temps_off, temps_spin, longitudes,
         latitudes, ras, decs, comp_widths, vel_diff, equiv_width, tau,
         filenames, local_paths, local_emission_paths, local_spectra_paths],
        names=['Day', 'Field', 'Source', 'Velocity', 'em_velocity', 'Optical_Depth', 'temp_off', 'temp_spin',
               'longitude', 'latitude', 'ra', 'dec', 'fwhm', 'vel_diff', 'equiv_width', 'tau',
               'Filename', 'Local_Path', 'Local_Emission_Path', 'Local_Spectrum_Path'],
        meta={'ID': 'magmo_gas',
              'name': 'MAGMO Gas ' + str(datetime.date.today())})
    votable = from_table(temp_table)
    table = votable.get_first_table()
    table.get_field_by_id('ra').ucd = 'pos.eq.ra;meta.main'
    table.get_field_by_id('dec').ucd = 'pos.eq.dec;meta.main'
    filename = "magmo-gas.vot"
    writeto(votable, filename)


def main():
    # Parse command line options
    args = parseargs()

    start = time.time()
    print("#### Started examining MAGMO components at %s ####" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    # Load component catalogue
    components = read_components("magmo-components.vot")
    spectra_map = read_spectra("magmo-spectra.vot")

    all_gas = analyse_components(components, spectra_map)

    # Output a catalogue
    output_gas_catalogue(all_gas)

    # Report
    end = time.time()
    print('#### Examination completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Processed %d components and output %d gas stats in %.02f s' %
          (len(components), len(all_gas), end - start))
    return 0


if __name__ == '__main__':
    exit(main())

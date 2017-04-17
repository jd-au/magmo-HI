#!/usr/bin/env python -u

# Take the Gaussian components and use them to examine the gas characteristics at each position.
#

# Author James Dempsey
# Date 26 Mar 2017


from __future__ import print_function, division

from astropy.io.votable import parse, from_table, writeto
from astropy.table import Table, Column
from string import Template

import argparse
import datetime
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


def analyse_components(components):
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

            if t_off > 0 and comp_amp < 0.98:
                # Calculate spin temperature and column density
                t_s = t_off / comp_amp

                # Record
                gas = Gas(component['Day'], component['Field'], component['Source'])
                gas.comp_vel = comp_vel
                gas.comp_amp = comp_amp
                gas.comp_width = comp_width
                gas.optical_depth = optical_depth
                gas.t_off = t_off
                gas.t_s = t_s
                gas.em_vel = em_vel
                gas.longitude = component['Longitude']
                gas.latitude = component['Latitude']
                #component['t_s '] = t_s
                print ("src %s at velocity %.4f has t_s %.3f (%.3f/%.3f)" % (component['Field'], comp_vel, t_s, t_off, comp_amp))
                all_gas.append(gas)
    return all_gas


def output_gas_catalogue(all_gas):
    days = []
    field_names = []
    sources = []
    longitudes = []
    latitudes = []
    velocities = []
    em_velocities = []
    optical_depths = []
    comp_widths = []
    temps_off = []
    temps_spin = []

    for i in range(len(all_gas)):
        gas = all_gas[i]
        days.append(gas.day)
        field_names.append(gas.field)
        sources.append(gas.src)
        longitudes.append(gas.longitude)
        latitudes.append(gas.longitude)
        velocities.append(gas.comp_vel)
        em_velocities.append(gas.em_vel/1000)
        optical_depths.append(gas.comp_amp)
        comp_widths.append(gas.comp_width)
        temps_off.append(gas.t_off)
        temps_spin.append(gas.t_s)

    temp_table = Table(
        [days, field_names, sources, velocities, em_velocities, optical_depths, temps_off, temps_spin, longitudes,
         latitudes, comp_widths],
        names=['Day', 'Field', 'Source', 'Velocity', 'em_velocity', 'Optical_Depth', 'temp_off', 'temp_spin',
               'longitude', 'latitude', 'fwhm'],
        meta={'ID': 'magmo_gas',
              'name': 'MAGMO Gas ' + str(datetime.date.today())})
    votable = from_table(temp_table)
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

    all_gas = analyse_components(components)

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

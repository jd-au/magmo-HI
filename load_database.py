#!/usr/bin/env python -u

# Load the MAGMO data into a database
#
# This program reads the previously generated fieldand spectra summaries and writes
# them out to a Postgres database.

# Author James Dempsey
# Date 9 May 2017

from __future__ import print_function, division

from astropy.io.votable import parse, from_table, writeto
from astropy.table import Table, Column

import psycopg2
import numpy as np
import argparse
import time


def select_all(table_name):
  # Establish the connection
  conn = psycopg2.connect(dbname='db', user='grok')
  cursor = conn.cursor()

  # Execute an SQL query and receive the output
  cursor.execute('Select * from ' + table_name)
  records = cursor.fetchall()
  return records


def parseargs():
    """
    Parse the command line arguments
    :return: An args map with the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Load the MAGMO data into a database")

    args = parser.parse_args()
    return args


def truncate_database():
    table_names = ['gas', 'spectrum', 'observation', 'field']

    # Establish the connection
    conn = psycopg2.connect(dbname='jamesdempsey', user='postgres')
    cursor = conn.cursor()

    cmd = 'truncate ' + ",".join(table_names)
    cursor.execute(cmd)
    conn.commit()


def read_votable_results(filename):
    votable = parse(filename, pedantic=False)
    results = next(resource for resource in votable.resources if
                   resource.type == "results")
    results_array = results.tables[0].array
    return results_array


def read_components(filename):
    return read_votable_results(filename)


def load_data_to_table(csv_file, table_name):
    # Establish the connection
    conn = psycopg2.connect(dbname='jamesdempsey', user='postgres')
    cursor = conn.cursor()

    cursor.copy_from(csv_file, table_name)
    print(cursor.statusmessage)
    cursor.execute("select count(*) from " + table_name)
    print(cursor.fetchall())
    conn.commit()


def load_observations():
    field_obs = read_votable_results("magmo-fields.vot")

    # extract list of fields
    field_names = np.unique(field_obs['Field'])
    # save as a csv
    filename = "database/fields.tsv"
    np.savetxt(filename, field_names, delimiter="\t", fmt="%s")
    # load
    with open(filename, 'r') as the_file:
        load_data_to_table(the_file, 'field')

    # Save observation data
    filename = "database/observations.tsv"
    np.savetxt(filename, field_obs, delimiter="\t", fmt="%s")
    with open(filename, 'r') as the_file:
        load_data_to_table(the_file, 'observation')


def load_spectra():
    field_obs = read_votable_results("magmo-spectra.vot")

    # Save observation data
    filename = "database/spectra.tsv"
    np.savetxt(filename, field_obs, delimiter="\t", fmt="%s")
    with open(filename, 'r') as the_file:
        load_data_to_table(the_file, 'spectrum')


def load_gas():
    field_obs = read_votable_results("magmo-gas.vot")

    # Save observation data
    filename = "database/gas.tsv"
    np.savetxt(filename, field_obs, delimiter="\t", fmt="%s")
    with open(filename, 'r') as the_file:
        load_data_to_table(the_file, 'gas')


def main():
    # Parse command line options
    args = parseargs()

    start = time.time()
    print("#### Started loading MAGMO database at %s ####" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    # Delete the old data
    truncate_database()

    # Load fields and observations
    load_observations()

    # Load spectra
    load_spectra()
    load_gas()


    # Report
    end = time.time()
    print('#### Load completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    #print('Processed %d components and output %d gas stats in %.02f s' %
    #      (len(components), len(all_gas), end - start))
    return 0


if __name__ == '__main__':
    exit(main())

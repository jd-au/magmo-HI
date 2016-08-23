# Set of general use functions for processing MAGMO HI data.
#

# Author James Dempsey
# Date 29 Jul 2016

import csv
import sys
import os
import subprocess


class CommandFailedError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


# functions to be made common
def get_day_file_data(day):
    """
    Read the magmo-days-full.csv file and find the row for the requested day.
    :param day: The day to be found
    :return: The day's row, or None if the day is not defined
    """

    with open('magmo-days-full.csv', 'rb') as magmodays:
        reader = csv.reader(magmodays)
        for row in reader:
            if row[0] == day:
                return row
    return None


def get_day_obs_data(day):
    """
    Read the magmo-obs.csv file and find the rows for the requested day.
    :param day: The day to be found
    :return: The day's rows, or None if the day is not defined
    """

    sources = []
    with open('magmo-obs.csv', 'rb') as magmo_obs:
        reader = csv.reader(magmo_obs)
        for row in reader:
            if row[0] == day:
                src = dict()
                src['source'] = row[4]
                src['phase_cal'] = row[10]
                src['gal_l'] = row[2]
                sources.append(src)
    return sources


def run_os_cmd(cmd, failOnErr=True):
    """
    Run an operating system command ensuring that it finishes successfully.
    If the comand fails, the program will exit.
    :param cmd: The command to be run
    :return: None
    """
    print ">", cmd
    sys.stdout.flush()
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0:
            message = "Command '"+cmd+"' failed with code " + str(retcode)
            print >>sys.stderr, message
            if failOnErr:
                raise CommandFailedError(message)
    except OSError as e:
        message = "Command '" + cmd + "' failed " + e
        print >> sys.stderr, message
        if failOnErr:
            raise CommandFailedError(message)
    return None


def ensure_dir_exists(dirname):
    """
    Check if a folder exists, and if it doesn't, create it. Fail the
    program if the folder could not be created.
    :param dirname: The name of the folder to be created.
    :return: None
    """
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    if not os.path.isdir(dirname):
        print "Directory %s could not be created." % dirname
        exit(1)
    return None

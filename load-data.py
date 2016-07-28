# Uses atlod to extract the HI spectra and the 2GHz continuum for a day
#
# This program extracts the HI spectra (1.42 GHz) and continuum (2 GHz) data
# from the day's RPFITS files into miriad format. It then backs up the metadata
# to allow any processing to be simply rolled back.

# Author James Dempsey
# Date 7 Jul 2016

import csv
import sys
import os
import shutil
import subprocess


# functions
def get_day_data(day):
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


def run_os_cmd(cmd):
    """
    Run an operating system command ensuring that it finishes successfully.
    If the comand fails, the program will exit.
    :param cmd: The command to be run
    :return: None
    """
    print ">", cmd
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0:
            print >>sys.stderr, "Command '"+cmd+"' failed with code", retcode
            exit(1)
    except OSError as e:
        print >>sys.stderr, "Command '"+cmd+"' failed:", e
        exit(1)


def load_rpfits(freq, ifsel, restFreq, dayNum, inFilePatterns):
    """
    Run atlod to extract data files for a particular frequency from a set of
    RPFITS into Miriad format.

    :param freq: The short title of the frequency
    :param ifsel: IF number(s) to select
    :param restFreq:  The rest frequency, in GHz, for line observations.
    :param dayNum: The MAGMO observation day number
    :param inFilePatterns: The patterns to match the RPFITS files.
    :return: None
    """
    outFile = "day" + dayNum + "/MAGMO_day" + dayNum + "_" + freq + ".uv"
    atlodCmd = "atlod in='" + inFilePatterns + "' out='" + outFile + "' " \
               + "ifsel=" + str(ifsel) + " restfreq=" + restFreq \
               + " options=birdie,xycorr,noauto,rfiflag"
    run_os_cmd(atlodCmd)


def ensure_dir_exists(dirname):
    """
    Check if a folder does not exists, and if it doesn't, create it. Fail the
    program if the folder could not be created.
    :param dirname: The name of the folder to be created.
    :return: None
    """
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    if not os.path.isdir(dirname):
        print "Directory %s could not be created." % dirname
        exit(1)


# ### Script starts here ###

# Read day parameter
if len(sys.argv) != 2:
    print("Incorrect number of parameters.")
    print("Usage: python load-data.py day")
    exit(1)
day = sys.argv[1]

# Read metadata for the day (file pattern fragments etc)
dayRow = get_day_data(day)
if dayRow is None:
    print "Day %s is not defined." % (day)
    exit(1)
print dayRow

# Make day directory
dayDirName = "day" + day
ensure_dir_exists(dayDirName)

# Load files
freqList = ["1421", "1720"]
ifList = [7, 4]
restFreqList = ["1.420405752", "1.72053"]
inFilePatterns = ""
startDatePatterns = 2
for i in range(startDatePatterns, len(dayRow)):
    inFilePatterns += "rawdata/" + dayRow[i] + "* "

for i in range(0, len(freqList)):
    freq = freqList[i]
    ifsel = ifList[i]
    restFreq = restFreqList[i]
    load_rpfits(freq, ifsel, restFreq, dayRow[0], inFilePatterns)

# Split files
os.chdir("day" + dayRow[0])
for freq in freqList:
    uvFile = "MAGMO_day" + dayRow[0] + "_" + freq + ".uv"
    uvsplitCmd = "uvsplit vis=" + uvFile
    run_os_cmd(uvsplitCmd)
os.chdir("..")

# Backup
backupDirName = dayDirName+"/backup"
ensure_dir_exists(backupDirName)
for uvDir in os.listdir(dayDirName):
    uvDirName = dayDirName + "/" + uvDir
    if (not uvDirName.endswith("backup")) and os.path.isdir(uvDirName):
        uvBackupDir = backupDirName + "/" + uvDir
        ensure_dir_exists(uvBackupDir)
        for name in ['flags', 'header', 'history']:
            shutil.copy2(uvDirName + "/" + name, uvBackupDir)

# Report

# Uses atlod to extract the HI spectra and the 2GHz continuum for a day
#
# This program extracts the HI spectra (1.42 GHz) and continuum (2 GHz) data
# from the day's RPFITS files into miriad format. It then backs up the metadata
# to allow any processing to be simply rolled back.

# Author James Dempsey
# Date 7 Jul 2016

import csv
import glob
import sys
import magmo
import os
import shutil
import time


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
    magmo.run_os_cmd(atlodCmd)


# ### Script starts here ###

# Read day parameter
if len(sys.argv) != 2:
    print("Incorrect number of parameters.")
    print("Usage: python load-data.py day")
    exit(1)
day = sys.argv[1]
start = time.time()

# Read metadata for the day (file pattern fragments etc)
dayRow = get_day_data(day)
if dayRow is None:
    print "Day %s is not defined." % (day)
    exit(1)
print "#### Started loading MAGMO day %s at %s ####" % \
      (day, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))
print dayRow

# Make day directory
dayDirName = "day" + day
magmo.ensure_dir_exists(dayDirName)

# Load files
freqList = ["1421", "1757"]
ifList = [7, 1]
restFreqList = ["1.420405752", "1.757"]
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
    magmo.run_os_cmd(uvsplitCmd)
# Rename any 1420.5 folders to remove the decimal
uv_dirs = glob.glob('*.[0-9][0-9][0-9][0-9].[0-9]')
for uvdir in uv_dirs:
    os.rename(uvdir, uvdir[:-2])
os.chdir("..")


# Backup
backupDirName = dayDirName+"/backup"
magmo.ensure_dir_exists(backupDirName)
for uvDir in os.listdir(dayDirName):
    uvDirName = dayDirName + "/" + uvDir
    if (not uvDirName.endswith("backup")) and (not uvDirName.endswith(".uv")) and os.path.isdir(uvDirName):
        uvBackupDir = backupDirName + "/" + uvDir
        magmo.ensure_dir_exists(uvBackupDir)
        for name in ['flags', 'header', 'history']:
            shutil.copy2(uvDirName + "/" + name, uvBackupDir)

# Cleanup
#for freq in freqList:
#    cmd = 'rm -rf day' + dayRow[0] + '/MAGMO_day' + dayRow[0] + '_' + freq + '.uv'
#    magmo.run_os_cmd(cmd)

# Report
end = time.time()
print '#### Loading completed at %s ####' \
      % (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
print 'Processed in %.02f s' % (end - start)

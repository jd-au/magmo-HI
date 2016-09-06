# Delete the raw data files for a day
#
# Author James Dempsey
# Date 6 Sep 2016

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


def main():
    # Read day parameter
    if len(sys.argv) != 2:
        print("Incorrect number of parameters.")
        print("Usage: python delete_raw.py day")
        exit(1)
    day = sys.argv[1]
    start = time.time()

    # Read metadata for the day (file pattern fragments etc)
    dayRow = get_day_data(day)
    if dayRow is None:
        print "Day %s is not defined." % (day)
        exit(1)
    print "#### Started deleting raw data for MAGMO day %s at %s ####" % \
          (day, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    print dayRow

    # Delete files
    startDatePatterns = 2
    for i in range(startDatePatterns, len(dayRow)):
        os.system("rm -r rawdata/" + dayRow[i] + "* ")

    # Report
    end = time.time()
    print '#### Deleting completed at %s ####' \
          % (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print 'Processed in %.02f s' % (end - start)

    return 0


# Run the script if it is called from the command line
if __name__ == "__main__":
    exit(main())
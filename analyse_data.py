
import magmo
import os
import sys
import time
import csv

from astropy.io import ascii


sn_min = 1.3
num_chan = 1053


def get_high_signal_fields(day_dir_name):
    """
    Retrieve a list of fields observed in a particular day that have sufficient signal
    to noise to search for background sources.
    :param day_dir_name: The name of the day's directory.
    :return: A list of high signal fields.
    """
    field_list = []
    print "Fields of interest:"
    with open(day_dir_name + '/stats.csv', 'rb') as stats:
        reader = csv.reader(stats)
        first = True
        for row in reader:
            if first:
                first = False
            else:
                if float(row[3]) > sn_min:
                    print row
                    field_list.append(row[0])

    return field_list


def find_sources(day_dir_name, field_name):
    """
    Search a contoinuum file for sources using the Aegean source finder. A
    VOTable file containing the list of discovered sources will be written out
    for the field. This function will use the Aegean source finder
    ( https://github.com/PaulHancock/Aegean ) to identify the sources.

    :param day_dir_name: The name of the day's directory.
    :param field_name:  The name fo the field to be searched for sources.
    :return: A list of error messages, if any
    """
    error_list = []
    cont_file = day_dir_name + "/1757/magmo-" + field_name + "_1757_restor.fits"
    table_file = day_dir_name + "/" + field_name + '_src.vot'
    try:
        print "##--## Searching continuum image " + cont_file + " ##--##"
        magmo.run_os_cmd('bane ' + cont_file)
        aegean_cmd = 'aegean ' + cont_file + ' --autoload --telescope=ATCA --cores=1 ' \
                     '--table=' + table_file
        magmo.run_os_cmd(aegean_cmd)
    except magmo.CommandFailedError as e:
        error_list.append(str(e))
    return error_list



def main():
    """
    Main script for src_find
    :return: The exit code
    """
    # Read day parameter
    if len(sys.argv) != 2:
        print("Incorrect number of parameters.")
        print("Usage: python analyse_data.py day")
        return 1
    day = sys.argv[1]
    start = time.time()

    # Check metadata against file system
    day_dir_name = "day" + day
    if not os.path.isdir(day_dir_name):
        print "Directory %s could not be found." % day_dir_name
        return 1

    print "#### Started source finding on MAGMO day %s at %s ####" % \
          (day, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))
    error_list = []

    # Read list of sources, filter for ones to be processed
    field_list = get_high_signal_fields(day_dir_name)

    # For each file, find the sources
    for field in field_list:
        error_list.extend(find_sources(day_dir_name, field))

    # Report
    end = time.time()
    print '#### Processing completed at %s ####' \
          % (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print 'Searched %d images in %.02f s' % (len(field_list),
                                               end - start)
    if len(error_list) == 0:
        print "Hooray! No errors found."
    else:
        print "%d errors were encountered:" % (len(error_list))
        for err in error_list:
            print err
    return 0


# Run the script if it is called from the command line
if __name__ == "__main__":
    exit(main())
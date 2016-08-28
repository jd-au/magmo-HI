# Remove the intermediate miriad files so as to conserve disk space.
# The fits, png, dat and log files used by later stages are retained.
#

# Author James Dempsey
# Date 28 Jul 2016

import sys
import os

# ### Script starts here ###

# Read day parameter
if len(sys.argv) != 2:
    print("Incorrect number of parameters.")
    print("Usage: python compress_data.py day")
    exit(1)
day = sys.argv[1]
dayDirName = "day" + day

# Delete the generated files
print "Removing intermediate Miriad files"
os.chdir(dayDirName)
for band in ['1420', '1757']:
    if os.path.exists(band):
        os.chdir(band)
        os.system("rm -r *beam")
        os.system("rm -r *dirty")
        os.system("rm -r *clean")
        os.system("rm -r *restor")
        os.system("rm -r *_ave")
        os.chdir('..')

os.chdir('..')

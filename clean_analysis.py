# Remove all analysis output files.
#
# Author James Dempsey
# Date 6 Aug 2016

import sys
import os
import shutil

# ### Script starts here ###

# Read day parameter
if len(sys.argv) != 2:
    print("Incorrect number of parameters.")
    print("Usage: python clean_analysis.py day")
    exit(1)
day = sys.argv[1]
dayDirName = "day" + day

# Delete the generated files
print "Removing analysis files"
os.chdir(dayDirName)
os.system("rm *plot.png")
os.system("rm *.vot")
os.system("rm *.xml")
os.system("rm *.pdf")

os.chdir('..')

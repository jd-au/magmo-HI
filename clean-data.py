# Remove all produced images and cubes, reset header, history and flags back
# to backed up files.
#

# Author James Dempsey
# Date 9 Jul 2016

import sys
import os
import shutil

# ### Script starts here ###

# Read day parameter
if len(sys.argv) != 2:
    print("Incorrect number of parameters.")
    print("Usage: python clean-data.py day")
    exit(1)
day = sys.argv[1]
dayDirName = "day" + day

# Delete the generated files
print "Removing generated files"
os.chdir(dayDirName)
for band in ['1420', '1757']:
    if os.path.exists(band):
        os.chdir(band)
        os.system("rm -r *beam")
        os.system("rm -r *dirty")
        os.system("rm -r *clean")
        os.system("rm -r *restor")
        os.system("rm -r *_ave")
        os.system("rm -r *.fits")
        os.system("rm -r *.png")
        os.system("rm -r *.png_2")
        os.system("rm -r *.pdf")
        os.chdir('..')

# Delete the analysis files
os.system("rm *.vot")
os.system("rm *.xml")
os.system("rm *.png")
os.system("rm *.png_2")
os.system("rm *.pdf")
os.chdir('..')

# Copy from the backups back to the orignal files
os.system("rm " + dayDirName + '/*/bandpass')
os.system("rm " + dayDirName + '/*/gains')
os.system("rm " + dayDirName + '/*/leakage')
backupDirName = dayDirName+"/backup"
for uvDir in os.listdir(backupDirName):
    print "Restoring backups for ", uvDir
    uvDirName = dayDirName + "/" + uvDir
    uvBackupDir = backupDirName + "/" + uvDir
    for bkFile in os.listdir(uvBackupDir):
        shutil.copy2(uvBackupDir+"/"+bkFile, uvDirName)

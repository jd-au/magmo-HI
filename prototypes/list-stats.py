# List stats for all images in a folder
# Assumes that clean images have the suffix _clean

# Author James Dempsey
# Date 29 Jun 2016

patternClean='.*_clean$'
patternRestored='.*_restor$'

import os
import re
import subprocess
import sys
import math

class SourceStats:
    source = ""
    rms = []
    max_sig = []

    def __init__(self):
        self.rms = [0,0]
        self.max_sig = [0,0]

def getSignalNoiseRatio(f):
        result = subprocess.check_output(["imstat", "in="+f, "options=guaranteespaces"])
        totalNext = False
        for line in result.splitlines():
            # print line
            if totalNext:
                stats = line.strip().split()
                #print stats
                #mean = float(stats[1])
                rms = float(stats[2])
                max = float(stats[3])
                return rms, max
            if re.match('^[ \t]*Total', line):
                totalNext = True;
        print "Unable to find totals in:"
        print result

def getSrcStats(allstats, srcName):
    srcStats = allstats.get(srcName, None)
    if (srcStats is None):
        srcStats = SourceStats()
        srcStats.source = srcName
        allstats[srcName] = srcStats
    return srcStats


# Scirpt starts here
# Read day parameter
if len(sys.argv) != 2:
    print("Incorrect number of parameters.")
    print("Usage: python list-stats.py day")
    exit(1)
day = sys.argv[1]
dayDirName = "day" + day

allstats = dict()
dir_name = dayDirName + "/1757/"
for f in os.listdir(dir_name):
    match = re.search('[0-9]+\.[0-9]+[+-][0-9.]+', f)
    if match is not None:
        srcName = match.group()
    if re.match(patternClean, f):
        #print "Found file %s for src %s" % (f, srcName)
        rms, max = getSignalNoiseRatio(dir_name+f)
        srcStats = getSrcStats(allstats, srcName)
        #print rms
        srcStats.rms[0] = rms
        srcStats.max_sig[0] = max
    if re.match(patternRestored, f):
        rms, max = getSignalNoiseRatio(dir_name+f)
        srcStats = getSrcStats(allstats, srcName)
        srcStats.rms[1] = rms
        srcStats.max_sig[1] = max

channels = 793
keys = allstats.keys()
keys.sort()
for srcName in keys:
    stats = allstats.get(srcName)
    sn_ratio = [0,0]
    sn_ratio[0] = stats.max_sig[0]/(1 if stats.rms[0] == 0 else stats.rms[0])
    sn_ratio[1] = stats.max_sig[1]/(1 if stats.rms[1] == 0 else stats.rms[1])
    print '%s rms: %.3e max: %.3e s/n: %.2e, s/n/chan: %.2e Rs/n: %.2e, Rs/n/chan: %.2e' % \
          (stats.source, stats.rms[0], stats.max_sig[0], sn_ratio[0], sn_ratio[0]/math.sqrt(channels), sn_ratio[1], sn_ratio[1]/math.sqrt(channels))

print '\n\n"File","RMS","MAX","Clean S/N","Clean S/N/Chan","Restored S/N","Restored S/N/Chan"'
for srcName in keys:
    stats = allstats.get(srcName)
    sn_ratio = [0,0]
    sn_ratio[0] = stats.max_sig[0]/(1 if stats.rms[0] == 0 else stats.rms[0])
    sn_ratio[1] = stats.max_sig[1]/(1 if stats.rms[1] == 0 else stats.rms[1])
    sn_chan = sn_ratio[1]/math.sqrt(channels)
    use = sn_chan > 3
    print '"%s",%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%s' % \
          (stats.source, stats.rms[0], stats.max_sig[0], sn_ratio[0], sn_ratio[0]/math.sqrt(channels), sn_ratio[1],
           sn_chan,("Use" if use else ""))

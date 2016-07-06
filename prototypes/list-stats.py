# List stats for all images in a folder
# Assumes that clean images have the suffix _clean

# Author James Dempsey
# Date 29 Jun 2016

patternClean='.*_clean$'
patternRestored='.*_restor$'

import os
import re
import subprocess
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


allstats = dict()
for f in os.listdir('.'):
    match = re.search('[0-9]+\.[0-9]+[+-][0-9.]+', f)
    if match is not None:
        srcName = match.group()
    if re.match(patternClean, f):
        #print "Found file %s for src %s" % (f, srcName)
        rms, max = getSignalNoiseRatio(f)
        srcStats = getSrcStats(allstats, srcName)
        #print rms
        srcStats.rms[0] = rms
        srcStats.max_sig[0] = max
    if re.match(patternRestored, f):
        rms, max = getSignalNoiseRatio(f)
        srcStats = getSrcStats(allstats, srcName)
        srcStats.rms[1] = rms
        srcStats.max_sig[1] = max

channels = 793
keys = allstats.keys()
keys.sort()
for srcName in keys:
    stats = allstats.get(srcName)
    sn_ratio = [0,0]
    sn_ratio[0] = stats.max_sig[0]/stats.rms[0]
    sn_ratio[1] = stats.max_sig[1]/stats.rms[1]
    print '%s rms: %.3e max: %.3e s/n: %.2e, s/n/chan: %.2e Rs/n: %.2e, Rs/n/chan: %.2e' % \
          (stats.source, stats.rms[0], stats.max_sig[0], sn_ratio[0], sn_ratio[0]/math.sqrt(channels), sn_ratio[1], sn_ratio[1]/math.sqrt(channels))

print '\n\n"File","RMS","MAX","Clean S/N","Clean S/N/Chan","Restored S/N","Restored S/N/Chan"'
for srcName in keys:
    stats = allstats.get(srcName)
    sn_ratio = [0,0]
    sn_ratio[0] = stats.max_sig[0]/stats.rms[0]
    sn_ratio[1] = stats.max_sig[1]/stats.rms[1]
    print '"%s",%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,' % \
          (stats.source, stats.rms[0], stats.max_sig[0], sn_ratio[0], sn_ratio[0]/math.sqrt(channels), sn_ratio[1], sn_ratio[1]/math.sqrt(channels))

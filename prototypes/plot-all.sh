#!/bin/csh

# Produce standard plots for all sources in a folder

# Author James Dempsey
# Date 28 Jun 2016

# Cleanup prev run
rm *-amp.ps
rm *-uv.ps

# Loop through all identified sources
set suffix = ".1420"
foreach srcname (232.620+0.996 254.880+0.451 259.939-0.041 263.250+0.514 264.140+2.018 264.289+1.469 269.153-1.128 269.456-1.467 269.658-1.270 270.255+0.835 281.710-1.104 284.352-0.419 284.694-0.361 285.337-0.002 286.383-1.834 287.371+0.644 290.374+1.661 290.411-2.915 291.270-0.719)
#foreach srcnam (347.817+0.018 347.902+0.052 348.195+0.768)
  echo "##--## Processing visibility ${srcname} ##--##"

  uvplt device=${srcname}-amp.ps/vcps options=nobase vis=${srcname}${suffix}
  uvplt device=${srcname}-uv.ps/vcps options=nobase axis=uc,vc vis=${srcname}${suffix}
end

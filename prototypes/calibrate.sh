#!/bin/csh

# Author James Dempsey
# Date 27 Jun 2016

# Calibrate the data using
# 1934-638 as the flux and bandpass cal and
# 0727-115, 0823-500, 1049-53 as the phase cals.

# Process the flux and bandpass cal
mfcal vis=1934-638.1420 options=interpolate
gpcal vis=1934-638.1420 options=xyvary

# Process each phase cal
foreach calnam (0727-115 0823-500 1049-53)
#foreach srcnam (347.817+0.018 347.902+0.052 348.195+0.768)
  echo "##--## Processing phase cal ${calnam} ##--##"

  gpcopy vis=1934-638.1420 out=${calnam}.1420
  gpcal vis=${calnam}.1420 options=xyvary,qusolv

  gpboot vis=${calnam}.1420 cal=1934-638.1420
  mfboot vis=${calnam}.1420,1934-638.1420 "select=source(1934-638)"

  echo "#### Validation ####"
  uvflux vis=1934-638.1420 stokes=i,q,u,v
  uvflux vis=${calnam}.1420 stokes=i,q,u,v
  uvplt vis=${calnam}.1420 stokes=i,q,u,v axis=real,imag options=equal,nobase device=/xs
  echo "gpplt vis=${calnam}.1420 device=/xs yaxis=phase options=xygains"
end

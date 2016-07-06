#!/bin/csh

# Produce spectral cubes for particular sources
# Author James Dempsey
# Date 27 Jun 2016

foreach srcnam (285.337-0.002 281.710-1.104)
#foreach srcnam (347.817+0.018 347.902+0.052 348.195+0.768)
  echo "##--## Processing ${srcnam} ##--##"
  rm -r ${srcnam}.1420.ave
  rm -r magmo-${srcnam}_1420_sl_dirty magmo-${srcnam}_1420_sl_beam
  rm -r magmo-${srcnam}_1420_sl_clean
  rm -r magmo-${srcnam}_1420_sl_restor magmo-${srcnam}_1420_sl*.fits

  uvaver vis=${srcnam}.1420 line=felocity,1053,-250.0,0.4,0.4 out=${srcnam}.1420.ave
  invert vis=${srcnam}.1420.ave map=magmo-${srcnam}_1420_sl_dirty beam=magmo-${srcnam}_1420_sl_beam cell=5 robust=0.5 options=systemp,nopol,mosaic,double stokes=i slop=0.5 line=felocity,1053,-250.0,0.4,0.4
  clean map=magmo-${srcnam}_1420_sl_dirty beam=magmo-${srcnam}_1420_sl_beam  out=magmo-${srcnam}_1420_sl_clean niters=2000 speed=+1
  restor model=magmo-${srcnam}_1420_sl_clean beam=magmo-${srcnam}_1420_sl_beam map=magmo-${srcnam}_1420_sl_dirty out=magmo-${srcnam}_1420_sl_restor
  fits in=magmo-${srcnam}_1420_sl_restor op=xyout out=magmo-${srcnam}_1420_sl_restor.fits
end

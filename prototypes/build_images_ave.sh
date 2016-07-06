#!/bin/csh

# Produce averaged continuum images for each named source
# Author James Dempsey
# Date 29 Jun 2016

foreach srcnam (281.710-1.104 285.337-0.002)
  echo "##--## Processing ${srcnam} ##--##"
  rm -r magmo-${srcnam}_1420_ave_dirty magmo-${srcnam}_1420_ave_beam
  rm -r magmo-${srcnam}_1420_ave_clean
  rm -r magmo-${srcnam}_1420_ave_restor magmo-${srcnam}_1420_ave_*.fits
  switch ($srcnam)
    case "232.620+0.996":
      gpcopy vis=0727-115.1420.ave out=${srcnam}.1420.ave
      breaksw
    case "254.880+0.451":
    case "259.939-0.041":
    case "263.250+0.514":
    case "264.140+2.018":
    case "264.289+1.469":
    case "269.153-1.128":
    case "269.456-1.467":
    case "269.658-1.270":
    case "270.255+0.835":
      gpcopy vis=0823-500.1420.ave out=${srcnam}.1420.ave
      breaksw
    default:
      gpcopy vis=1049-53.1420.ave out=${srcnam}.1420.ave
      breaksw
  endsw

  invert robust=0.5 options=systemp,mfs,double stokes=ii vis=${srcnam}.1420.ave map=magmo-${srcnam}_1420_ave_dirty beam=magmo-${srcnam}_1420_ave_beam
  clean map=magmo-${srcnam}_1420_ave_dirty beam=magmo-${srcnam}_1420_ave_beam out=magmo-${srcnam}_1420_ave_clean niters=2000 speed=+1
  restor model=magmo-${srcnam}_1420_ave_clean beam=magmo-${srcnam}_1420_ave_beam map=magmo-${srcnam}_1420_ave_dirty options=mfs out=magmo-${srcnam}_1420_ave_restor
  fits in=magmo-${srcnam}_1420_ave_restor op=xyout out=magmo-${srcnam}_1420_ave_restor.fits
end

#!/bin/csh

# Produce continuum images for each named source
# Author James Dempsey
# Date 27 Jun 2016

foreach srcnam (232.620+0.996 254.880+0.451 259.939-0.041 263.250+0.514 264.140+2.018 264.289+1.469 269.153-1.128 269.456-1.467 269.658-1.270 270.255+0.835 281.710-1.104 284.352-0.419 284.694-0.361 285.337-0.002 286.383-1.834 287.371+0.644 290.374+1.661 290.411-2.915 291.270-0.719)
#foreach srcnam (347.817+0.018 347.902+0.052 348.195+0.768)
  echo "##--## Processing ${srcnam} ##--##"
  rm -r magmo-${srcnam}_1420_dirty magmo-${srcnam}_1420_beam
  rm -r magmo-${srcnam}_1420_clean
  rm -r magmo-${srcnam}_1420_restor magmo-${srcnam}_1420_*.fits
  switch ($srcnam)
    case "232.620+0.996":
      gpcopy vis=0727-115.1420 out=${srcnam}.1420
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
      gpcopy vis=0823-500.1420 out=${srcnam}.1420
      breaksw
    default:
      gpcopy vis=1049-53.1420 out=${srcnam}.1420
      breaksw
  endsw

  invert robust=0.5 options=systemp,mfs,double stokes=ii vis=${srcnam}.1420 map=magmo-${srcnam}_1420_dirty beam=magmo-${srcnam}_1420_beam
  clean map=magmo-${srcnam}_1420_dirty beam=magmo-${srcnam}_1420_beam out=magmo-${srcnam}_1420_clean niters=2000 speed=+1
  restor model=magmo-${srcnam}_1420_clean beam=magmo-${srcnam}_1420_beam map=magmo-${srcnam}_1420_dirty options=mfs out=magmo-${srcnam}_1420_restor
  fits in=magmo-${srcnam}_1420_restor op=xyout out=magmo-${srcnam}_1420_restor.fits
end

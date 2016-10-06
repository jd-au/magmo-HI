#!/bin/csh

# Produce continuum images at 1420 and 1757 for a named source
# Author James Dempsey
# Date 7 Oct 2016

foreach srcnam (1740-517)
  foreach freq in (1420, 1757)
      echo "##--## Processing ${srcnam} at ${freq} ##--##"
      rm -r magmo-${srcnam}_{freq}_dirty magmo-${srcnam}_{freq}_beam
      rm -r magmo-${srcnam}_{freq}_clean
      rm -r magmo-${srcnam}_{freq}_restor magmo-${srcnam}_{freq}_*.fits

      # Presume gpcopy already done

      invert robust=0.5 options=systemp,mfs,double stokes=ii vis=${srcnam}.{freq} map=magmo-${srcnam}_{freq}_dirty beam=magmo-${srcnam}_{freq}_beam
      clean map=magmo-${srcnam}_{freq}_dirty beam=magmo-${srcnam}_{freq}_beam out=magmo-${srcnam}_{freq}_clean niters=2000 speed=+1
      restor model=magmo-${srcnam}_{freq}_clean beam=magmo-${srcnam}_{freq}_beam map=magmo-${srcnam}_{freq}_dirty options=mfs out=magmo-${srcnam}_{freq}_restor
      fits in=magmo-${srcnam}_{freq}_restor op=xyout out=magmo-${srcnam}_{freq}_restor.fits
  end
end

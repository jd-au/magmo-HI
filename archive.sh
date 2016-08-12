#!/usr/bin/env bash

daydir='day21'
today=`date +%Y%m%d`
archivedir="prev_runs/${daydir}-${today}"
# Copy all
mkdir ${archivedir}
mkdir  ${archivedir}/1420
mkdir  ${archivedir}/1757
cp day21/*.html ${archivedir}/
cp day21/*.log ${archivedir}/
cp day21/*.png* ${archivedir}/
cp day21/1420/*.png* ${archivedir}/1420/
cp day21/1757/*.png* ${archivedir}/1757/

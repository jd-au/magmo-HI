#!/usr/bin/env bash

day=$1
daydir="day${day}"
today=`date +%Y%m%d`
archivedir="prev_runs/${daydir}-${today}"
# Copy all
mkdir ${archivedir}
mkdir  ${archivedir}/1420
mkdir  ${archivedir}/1757
cp -a ${daydir}/*.html ${archivedir}/
cp -a ${daydir}/*.log ${archivedir}/
cp -a ${daydir}/*.png* ${archivedir}/
cp -a ${daydir}/1420/*.png* ${archivedir}/1420/
cp -a ${daydir}/1757/*.png* ${archivedir}/1757/

#!/usr/bin/env bash

day=$1
daydir="day${day}"
today=${2-`date +%Y%m%d`}
archivedir="prev_runs/${today}/${daydir}"
echo "Archiving ${daydir} to ${archivedir}"
# Copy all
mkdir "prev_runs/${today}"
mkdir ${archivedir}
mkdir  ${archivedir}/1420
mkdir  ${archivedir}/1757
cp -a ${daydir}/*.html ${archivedir}/
cp -a ${daydir}/*.log ${archivedir}/
cp -a ${daydir}/*.png* ${archivedir}/
cp -a ${daydir}/*.xml ${archivedir}/
cp -a ${daydir}/*.vot ${archivedir}/
cp -a ${daydir}/*.csv ${archivedir}/
cp -a ${daydir}/1420/*.png* ${archivedir}/1420/
cp -a ${daydir}/1757/*.png* ${archivedir}/1757/

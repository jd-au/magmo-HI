#!/usr/bin/env bash

# Run the processing steps for a day

if [[ -z "$1" ]] ; then
    echo "Usage: reprocess.sh day [archivefoldername]"
    exit 1
fi

day=$1
daydir="day${day}"

date
echo "Reprocessing day ${day}"
python clean-data.py ${day}

echo "Processing ... "
python process_data.py ${day} >& ${daydir}/process.log

echo "Imaging ... "
python image-1420.py ${day} >& ${daydir}/image.log

echo "Analysing ... "
python analyse_data.py ${day} >& ${daydir}/analyse.log

echo "Archiving ..."
./archive.sh ${day} ${2-`date +%Y%m%d`}

echo "Done"
date
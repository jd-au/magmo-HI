#!/usr/bin/env bash

# Script to download all files from the MAGMO survey and extract the IFs of interest

#for i in `seq 6 43`; do
for i in `seq 8 20`; do
    echo "Downloading day $i"
    daydir = day${i}
    cd rawdata
    wget -c --tries=4 -i ../filelist/day${i}.txt
    cd ..
    python load-data.py ${i} >& logs/load-d${i}.log &
done

echo "Done"

#!/usr/bin/env bash

# Script to download all files from the MAGMO survey and extract the IFs of interest

#for i in `seq 6 43`; do
for i in `seq 6 20`; do
    echo "Downloading day $i"
    python find_data.py ${i} 'james.dempsey@csiro.au' .credentials
    python load-data.py ${i} >& load-d${i}.log
done

echo "Done"

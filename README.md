# MAGMO HI

Scripts for processing the HI observations taken with ATCA during the MAGMO project 

## Directory layout

- readme.md
- magmo-obs.csv (table of days and date prefixes plus any other needed metadata)
- load_day.py
- process-day.py
- analyse-day.py
- clean-day.py
- analyse-spectra.py
- rawdata
    - RPFITS files
- day99
    - 1420 (or 1421)
    - 1757

## Python programs

- load_day.py day
    - Uses atlod to extract the HI spectra and the 2GHz continuum
    - Backup the flag, header and history files for each source
- process-day.py day
    - flag, split, calibrate, produce 2 GHz continuum image, produce HI image cube
- analyse-day.py day
    - Source find on continuum images
    - For each point of interest extract spectra, convert to opacity and output numerical spectrum
- clean-day.py day
    - Remove all produced images and cubes, reset header, history and flags back to backed up files.
- analyse-spectra.py
    - Read in all numerical spectra and produce reports
    - Longitude-velocity diagram
    - Histograms
    - Catalogue of HI regions
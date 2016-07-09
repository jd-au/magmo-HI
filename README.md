# MAGMO HI

Scripts for processing the HI observations taken with ATCA during the MAGMO project 

## This Project

This code is part of an honours project being run out of the Australian 
National University's Research School of Astronomy and Astrophysics (RSAA).

In this research project I intend to analyse the composition and 
structure of the cold and warm HI gas at the locations of the MAGMO 
project observations (Green et al. , 2012). To enable this analysis I 
will produce HI absorption and emission spectra from the MAGMO observations.

-- | --
| Project Title | The Cold Neutral Medium As Seen By The MAGMO Project |
| Author | James Dempsey |
| Supervisor | Naomi McClure-Griffiths |


## MAGMO

The MAGMO project was an Australia Telescope Compact Array (ATCA) observing 
program run from 2010 to 2012 which had the aim of studying “Magnetic 
fields of the Milky Way through OH masers” (Green et al. , 2012). It was 
primarily concerned with measuring Zeeman splitting of hydroxyl (OH) 
ground state spectral lines. Many of the OH masers observed by the MAGMO 
project were associated with 6.7 GHz methanol masers. These emissions are 
linked to high mass star formation (HMSF) regions (Minier et al. , 2003). 
As a secondary objective, neutral hydrogen observations were taken towards 
most of the same sources (Green et al. , 2010). As a result these HI 
observations will be largely of HMSF regions.


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
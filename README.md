# MAGMO HI

Scripts for processing the HI observations taken with ATCA during the MAGMO project 

## This Project

This code is part of an honours project that was run out of the Australian 
National University's Research School of Astronomy and Astrophysics (RSAA).

In this research project I analysed the composition and 
structure of the cold and warm HI gas at the locations of the MAGMO 
project observations (Green et al. , 2012). To enable this analysis I 
produced HI absorption and emission spectra from the MAGMO and SGPS observations.


| -- | -- |
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
- magmo-flagging.csv (List of data blocks to be flagged, e.g. due to RFI)
- *.py (Python programs as described below)
- *.sh (Bash scripts to run the pipeline across all data)
- prototypes (Various proof of concept and exploratory programs)
- rawdata
    - RPFITS files
- day99
    - 1420 (or 1421)
    - 1757

## Python programs

### MAGMO Pipeline

The parts of the pipeline for processing magmo data, in processing order is as follows:

- prep-days.py
    - Prepares an extended list of days with the patterns to match the RPFITS files.
- load_data.py day
    - Python program to control the Miriad data reduction process for a specific day
    - Uses the atlod to extract the HI spectral data and the 2GHz continuum data
    - Split out the uv data into sources (uvsplit)
    - Backup the flag, header and history files for each source
- process_day.py day
    - Flag, calibrate, produce 2 GHz continuum image, produce HI image cube
- analyse_data.py day
    - Source find on continuum images
    - For each point of interest extract spectra, convert to opacity and output numerical spectrum
    - Extract emission spectrum for each point of interest from the SGPS data
- analyse_spectra.py
    - Read in all numerical spectra and produce reports
    - Longitude-velocity diagram
    - Histograms
    - Catalogue of HI regions
- decompose.py
    - Decompose each of the MAGMO spectra into Gaussian components using gausspy
- examine_gas.py
    - Calculate details such as spin temperature, column density, cold gas fraction for each as component
    - Produce population histograms for each gas statistic
    - Produce longitude velocity plots of spin temperature and FWHM
    - Compare population stats with previous studies
    - Compare specific spectra with matching spectra from previous SGPS studies 

The following scripts were used to control the pipeline

- archive.sh day date
    - Copy the pipeline outputs for a specific day to a backup folder. 
    - Backup folders are named by the current date prev_run/YYMMDD or can be overridden by the data parameter
- reprocess.sh day [archivefoldername]
    - Run the processing steps for a day
    - Archive the results in either a default directory based on the current date, or a named directory
- find-all.sh
    - Find the data files in ATOA for all days of MAGMO observing
- download-all.sh
    - Download the data files from ATOA for all days of MAGMO observing


### Common Libraries

- magmo.py
    - Utility functions for processing magmo HI data
- sgps.py
    - Shared library of useful functions for working with SGPS data


### Supporting Code

In addition, there are a number of utility programs
 
- clean_analysis.py day
    - Remove all analysis output files for a day.
- clean-data.py day
    - Remove all produced images and cubes, reset header, history and flags back to backed up files.
- compress_data.py
    - Remove the intermediate Miriad files so as to conserve disk space.
    - The fits, png, dat and log files used by later stages are retained.
- delete_raw.py
    - Delete the raw rpfits files for a day
- find_data.py
    - Query ATOA using the Table Access Protocol interface to find all MAGMO-HI data
    - Create input files for wget to retrieve the actual rpfits files fro ATOA
- image-1420.py
    - Produce images of the HI 1420 data
    - Provide a side by side comparison with the 1757 MHz continuum data
- record_cubes.py
    - Print stats on which fields had cubes produced or had spectra used.

    
### Experimental Database Code

During the initial MAGMO-HI project I experimented with creating a postgres database 
of the observations, spectra and gas characteristics. 
The scripts are not up to date with the current output structures from the pipeline however.  

- create_database.sql
  - SQL create table script to create the database structure
- load_database.py
  - Script to read in the VO tables from the pipeline and populated the database with the data
  - Covers observations, spectra and gas components

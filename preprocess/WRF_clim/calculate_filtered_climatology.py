"""
Filename:    calculate_filtered_climatology.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: preprocess chosen variables from 36 years of wrf output and save as yearly nc files
"""

## Imports
# Standard Python modules
import os, sys

# import personal modules

# Path to modules
sys.path.append('../modules')

# Import my modules
from wrf_funcs_preprocess import wrf_3hr_to_daily
from harmonics import calc_filtered_climatology

## variable/data to process
domain = 'd01'
varname = 'geopotential'

wrf_3hr_to_daily(varname, domain)
calc_filtered_climatology(domain, varname)
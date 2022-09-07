"""
Filename:    wrf_funcs_domain.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Functions used to get domain information for WRF
"""

## Imports

import os, sys
import xarray as xr
import wrf
import pandas as pd
import netCDF4 as nc
import numpy as np

# domain calculations
def check_multiple_of_parent_grid_ratio(parent_grid_ratio, n):
    '''check if e_we for inner domain is multiple of 3 plus one'''
    # For nest 2, (e_we-s_we+1) must be one greater than an integer multiple of the parent_grid_ratio of 3.
    base = parent_grid_ratio
    nearest_multiple = base * round(n/base)
    
    return int(nearest_multiple+1)

def calc_gridpoints_domain(domains, resolution, parent_grid_ratio):
    '''calculate domain gridpoints in east-west and north-south directions''' 
    lon2, lon1, lat2, lat1 = domains[0]
    ref_lat = (lat1 + lat2) / 2. # center lat of the parent domain
    ref_lon = (lon1 + lon2) / 2. # center lon of the parent domain
    
    e_we = []
    e_sn = []
    for i, (d, res) in enumerate(zip(domains, resolution)):
        lon2, lon1, lat2, lat1 = d
        ratio = parent_grid_ratio[i]
        
        if i == 0:
            e_we.append(round((lon1 - lon2) * 111.2 / res))
            e_sn.append(round((lat1 - lat2) * 111.2 / res)) 
        if i > 0:    
            e_we_grd = round((lon1 - lon2) * 111.2 / res)
            e_we.append(check_multiple_of_parent_grid_ratio(ratio, e_we_grd))
            e_sn_grd = round((lat1 - lat2) * 111.2 / res)
            e_sn.append(check_multiple_of_parent_grid_ratio(ratio, e_sn_grd))
    
    return e_we, e_sn, ref_lat, ref_lon

def calc_i_j_parent_start(domains, resolution):
    '''Function to calculate i and j parent start - potentially need to fix for more than 2 domains?'''
    i_parent_start = [1, ]
    j_parent_start = [1, ]
    
    lonmin_d01 = domains[0][0] # lonmin of parent domain
    latmin_d01 = domains[0][2] # latmin of parent domain
    res = resolution[0] # parent domain resolution
    
    for i, d in enumerate(domains[1:]):
        lonmin, lonmax, latmin, latmax = d
    
        i_parent_start.append(round((lonmin - lonmin_d01) * 111.2 / res))
        j_parent_start.append(round((latmin - latmin_d01) * 111.2 / res))
    
    return i_parent_start, j_parent_start

def calc_number_nodes_pes(e_we, e_sn, cores=32):
    '''
    This script finds the largest number of processors and nodes you can use, 
    based on the number of grid points in the i/j directions on your domain. 
    Note: The largest number may not decompose the best way.
    
    Parameters
    ----------
    e_we : int
        integer of the number of grid points in the i direction on your domain
    e_sn : int
        integer of the number of grid points in the j direction on your domain
    cores : int
        number of cores per node on your machine (default=32 - the max on one BW node)
  
    Returns
    -------
    weighted_arrays : list of arrays
        weights equal to the sqrt of the cosine of latitude
        
    ''' 
    #     # Basic estimate of range 
    #     max_proc_est = (small_domain_ewe/25.) * (small_domain_esn/25.)
    #     min_proc_est = (large_domain_ewe/100.) * (large_domain_esn/100.)
    #     print('min processesors is ', min_proc_est)
    #     print('max processesors is ', max_proc_est)

    #     print('min nodes is ', min_proc_est/cores)
    #     print('max nodes is ', max_proc_est/cores)

    # The value for 'cores' gets incremented later, so we want a static variable for the original value 
    cores_orig = cores

    # set upper limit of nodes - the max you want to loop through
    node_max = 200 

    # This is the least number of grid points allowed for each processor. 
    # Dont' change this value.
    smallest_size = 10

    x = 1
    while x <= node_max:

    # finds the factor pairs for the total number of cores
        def f(cores):
            factors = []
            for i in range(1, int(cores**0.5)+1):
                if cores % i == 0:
                    factors.append((i, cores/i ))
            return factors

        factors = f(cores)

    # Of the factor pairs, this finds the closest values (pair) in that array
        closest_factors = factors[-1]

    # Of the set of closest values, assign the i and j values
        i_array_value = closest_factors[0]
        j_array_value = closest_factors[-1]

    # Calculate how the domain will be decomposed
        e_we_decomp = int(e_we / i_array_value )
        e_sn_decomp = int(e_sn / j_array_value )

    # Once the decomposition becomes smaller than the least number of grid points
    # allowed for each processor, the loop will quit and display the max 
    # number of processors and nodes you can use for your domain.
        if ((e_sn_decomp < smallest_size) or (e_we_decomp < smallest_size)):

    # test to see if the max number of processors allowed is within the number for a single node 
            initial_factor_pair = factors[0]
            initial_factor = initial_factor_pair[-1]
            if initial_factor == cores_orig:

    # start with value of cores_orig and decrease by 1 for each iteration
    # until the value is allowed
               y = cores_orig
               while y >= 1:
                    processors = y

    # finds the factor pairs for the total number of processors
    # still testing processor values for a single node
                    def f(processors):
                        factors = []
                        for i in range(1, int(processors**0.5)+1):
                            if processors % i == 0:
                                factors.append((i, processors/i ))
                        return factors

                    factors = f(processors)

    # Of the factor pairs, this finds the closest values (pair) in that array
    # still testing processor values for a single node
                    closest_factors = factors[-1]

    # Of the set of closest values, assign the i and j values
    # still testing processor values for a single node
                    i_array_value = closest_factors[0]
                    j_array_value = closest_factors[-1]

    # Calculate how the domain will be decomposed
    # still testing processor values for a single node
                    e_we_decomp = int(e_we / i_array_value )
                    e_sn_decomp = int(e_sn / j_array_value )

    # Once the decomposition becomes larger or equal to the least number of grid points
    # allowed for each processor, the loop will quit and display the max 
    # number of processors and nodes you can use for your domain.
                    if ((e_sn_decomp >= smallest_size) and (e_we_decomp >= smallest_size)): 
                        max_procs = (i_array_value * j_array_value)
                        max_nodes = (max_procs / cores_orig)
#                         print("max # of processors that can be used is: ", max_procs)
#                         print("max # of nodes that can be used is 1 ")
                        break

    # if you haven't reached your limit, the loop continues
    # still testing processor values for a single node   
                    else:
                        y -= 1

    # if the size of the domain allows multiple nodes
            else:
                max_procs = (i_array_value * j_array_value) - cores_orig
                max_nodes = (max_procs / cores_orig)
#                 print("max # of processors that can be used is: ", max_procs)
#                 print("max # of nodes that can be used is: ", max_nodes)
            break

    # If you haven't reached your limit, the loop continues    
        x += 1
        cores = (cores+cores_orig)
    
    return int(max_nodes)

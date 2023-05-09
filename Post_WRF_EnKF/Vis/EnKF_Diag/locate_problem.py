#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import pickle
import netCDF4 as nc
import math
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time

diag_var = ['obs_type','lat','lon','height','i_grid','j_grid','k_grid','Hroi','Vroi','obs','obs_error_variance','prior_mean','posterior_mean','prior_spread','posterior_spread']

# ------------------------------------------------------------------------------------------------------
#           Operation: Read & process & return diagnostics in fort.10000
# ------------------------------------------------------------------------------------------------------
def read_diag( filename, v_interest ):

    # Initiate list containers for variables of interest
    list_all = [[] for i in range(len(v_interest))]

    # Read in values of each variable
    with open(filename, 'r') as f:
        all_lines = f.readlines()
    for line in all_lines:
        split_line = line[10:].split()
        for var in v_interest:
            idx_v_interest = v_interest.index( var )
            idx_v_indiag = diag_var.index( var )
            if var == 'obs_type':
                # read in the obs type value and get rid of the trailing whitespaces   
                list_all[idx_v_interest].append( line[0:10].strip() )
            else:
                # read in other values
                list_all[idx_v_interest].append( float(split_line[idx_v_indiag-1]) )

    # Assemble a dictionary 
    d_diag = {}
    for var in v_interest:
        v_idx = v_interest.index( var )
        d_diag[var] = list_all[v_idx]

    # Print a record to make sure reading process is correct
    #print('Sanity check of reading process...Printing 888th record...') 
    #for key in d_diag:
    #    print(key, ':', d_diag[key][888])

    return d_diag

def find_nearTarget( filename, v_interest, target_loc):

    # location of the target obs
    target_lat = target_loc[0]
    target_lon = target_loc[1]

    # get fort.10000 content
    d_diag = read_diag( filename, v_interest )

    # sort out the content by lat then lon
    merged_tuple = [(d_diag['obs_type'][i], d_diag['lat'][i],d_diag['lon'][i]) for i in range(0, len(d_diag['obs_type']))]
    sorted_tuple = sorted( merged_tuple, key=lambda x: (x[1],x[2]) )

    list_near = []
    # find the possible candidates
    itv = 3
    for it in sorted_tuple:
        if ( target_lat-itv <= it[1] <= target_lat+itv and target_lon-itv <= it[2] <= target_lon+itv):  
            list_near.append( it )
    print(list_near)


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'

    # configuration
    Storm = 'MARIA'
    filename = '/scratch/06191/tg854905/Pro2_PSU_MW/MARIA/J_DA+J_WRF+J_init-SP-intel17-THO-24hr-hroi900/run/201709161700/enkf/d01/fort.10000'
    v_interest = ['obs_type','lat','lon']
    target = [37.80455,-69.7223] # the location (lat,lon) where PSFC value is unphysical


    find_nearTarget( filename, v_interest, target )

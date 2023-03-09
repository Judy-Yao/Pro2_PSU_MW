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

# if a number ends with 5, round down as in Fortran
def round_half_down(n, decimals=2):
    if n > 0:
        n_int = n * 10 ** 3
        rem = n_int%10
        if rem == 0:
            return n
        else:
            multiplier = 10 ** decimals
            return (int(n * multiplier)) / multiplier
    if n < 0:
        n_int = n * 10 ** 3
        rem = n_int%10
        if rem == 0:
            return n
        else:
            multiplier = 10 ** decimals
            return (int(n * multiplier)-1) / multiplier
    
    #    multiplier = 10 ** decimals
    #    return (n_int / 10 - 0.5) / multiplier
    #elif rem > 5:
    #    return math.ceil(n_int)

def round_half_up(n, decimals=2):
    if n > 0:
        n_int = n * 10 ** 3
        rem = n_int%10
        if rem == 0:
            return n
        else:
            multiplier = 10 ** decimals
            return (int(n * multiplier)+1) / multiplier
    if n < 0:
        n_int = n * 10 ** 3
        rem = n_int%10
        if rem == 0:
            return n
        else:
            multiplier = 10 ** decimals
            return (int(n * multiplier)) / multiplier


    #multiplier = 10 ** decimals
    #return (n * multiplier + 0.5) / multiplier

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
                list_all[idx_v_interest].append( split_line[idx_v_indiag-1] )

    # Assemble a dictionary 
    d_diag = {}
    for var in v_interest:
        v_idx = v_interest.index( var )
        d_diag[var] = list_all[v_idx]

    # Print a record to make sure reading process is correct
    print('Sanity check of reading process...Printing 888th record...') 
    for key in d_diag:
        print(key, ':', d_diag[key][888])

    return d_diag

def return_obs_type( filename ):
    
    obs_type_all = []

    with open(filename, 'r') as f:
        all_lines = f.readlines()
    for line in all_lines:
        # read in the obs type value and get rid of the trailing whitespaces
        obs_type_all.append( line[0:10].strip() )

    # Find the unique obs type values
    obs_type_uni_set = set(obs_type_all)
    obs_type_uni = sorted(list(obs_type_uni_set))
    
    print('Unique obs types are: ')
    print(obs_type_uni)

# ---------------------------------------------------------------------------------------------------------
#           Operation: Assign microwave records from diagnostics with channel number from microwave_so file
# ---------------------------------------------------------------------------------------------------------
def read_mw_so( MW_SO ):

    # Read the content inside the microwave obs to a list
    with open(MW_SO) as f:
        all_lines = f.readlines()
    
    # Declare empty lists
    sensor_all = []
    chnum_all = []
    lat_all_max = []
    lon_all_max = []
    lat_all_min = []
    lon_all_min = []
    obs_all = []

    # Read in values
    for line in all_lines:
        split_line = line.split()
        sensor_all.append( split_line[1] )
        chnum_all.append( split_line[2] )
        # round digits after point to match diagnostics
        lat_all_min.append( "{0:.2f}".format(round_half_down(float(split_line[3])))) 
        lon_all_min.append( "{0:.2f}".format(round_half_down(float(split_line[4])))) 
        lat_all_max.append( "{0:.2f}".format(round_half_up(float(split_line[3])))) 
        lon_all_max.append( "{0:.2f}".format(round_half_up(float(split_line[4]))))
        obs_all.append( "{0:.2f}".format(float(split_line[5])))

    # Create a list tuples and sort
    mwso_tup = []
    for i in range(len(sensor_all)):
        mwso_tup.append( (sensor_all[i],chnum_all[i],lat_all_min[i],lat_all_max[i],lon_all_min[i],lon_all_max[i],obs_all[i]) )
    mwso_tup.sort(key=lambda a: (a[2], a[4]))
    return mwso_tup


def label_mw_obs( file_Diag, MW_SO, v_interest ):

    # Read the diagnostics from fort.10000 (obs assmilated by enkf)
    d_diag = read_diag( file_Diag, v_interest )
    # Get Microwave part 
    idx_mw = []
    irc = 0
    for iobs in d_diag['obs_type']:
        if d_diag['obs_type'][irc] == 'Microwave':
            idx_mw.append( irc )
        irc = irc + 1
    # Create a list of tuples
    Diag_tup = []    
    for idx in idx_mw: 
        Diag_tup.append( (d_diag['lat'][idx], d_diag['lon'][idx], d_diag['obs'][idx], d_diag['prior_mean'][idx], d_diag['posterior_mean'][idx]) )
    # Sort the tuples as lat, lon to speed up the search and match at next step
    Diag_tup.sort(key=lambda a: (a[0], a[1]))
    Diag_list = [ list( idia ) for idia in Diag_tup ]  

    # Read recrods from microwave_so file
    mwso_tup = read_mw_so( MW_SO )
    mwso_list = [list( iso ) for iso in mwso_tup]

    # For each record in Diag_list, find the match in mwso_list
    idx_match = []
    for idi in Diag_list:
        diag_lat = idi[0]
        diag_lon = idi[1]
        diag_obs = idi[2]
        #print('----------')       
        #print(diag_lat)
        #print(diag_lon)
        #print(diag_obs)
        irec = 0
        Get_match = False
        for iso in mwso_list:
            #print( iso[2],iso[3])
            #print( iso[4],iso[5])
            #print( iso[6] )
            if iso[6] != diag_obs:
                irec = irec + 1
                continue
            else:
                lat_log = diag_lat == iso[2] or diag_lat == iso[3]
                lon_log = diag_lon == iso[4] or diag_lon == iso[5]
                if lat_log and lon_log:
                    #print('idx: ', irec)
                    idx_match.append( irec )
                    Get_match = True
                    break
                else:
                    irec = irec + 1
        if Get_match == False:
            print('Not able to find the match!')

    # Sanity check
    if len(idx_mw) != len(idx_match):
        raise ValueError('Error matching records in Diag_bundle and so_bundle!')
    print('Check the 500th record in Diagnostics and the according record in MW SO...')
    print('In the Diagnostics:', Diag_list[500])
    print('In the MW so:', mwso_list[idx_match[500]])

    # Label the Diag_list with 
    for irec in range(len(Diag_list)):
        Diag_list[irec].insert(0, mwso_list[idx_match[irec]][1] )
        Diag_list[irec].insert(0, mwso_list[idx_match[irec]][0] )
    print('Check the 500th record in Diagnostics now...')
    print(Diag_list[500]) 

    # Reconstruct the lists
    ss = []
    chnum= []
    lat = []
    lon = []
    obs = []
    prior = []
    posterior = []
    for it in Diag_list:
        ss.append( it[0] )
        chnum.append( it[1] )
        lat.append( float(it[2]) )
        lon.append( float(it[3]) )
        obs.append( float(it[4]) )
        prior.append( float(it[5]) )
        posterior.append( float(it[6]) )
    
    Diag_list_reshape = [ss, chnum, lat, lon, obs, prior, posterior]

    return Diag_list_reshape

# ---------------------------------------------------------------------------------------------------------
#           Operation: Find Radiance records from diagnostics 
# ---------------------------------------------------------------------------------------------------------

def Find_IR( file_Diag, v_interest ):
    # Read the diagnostics from fort.10000 (obs assmilated by enkf)
    d_diag = read_diag( file_Diag, v_interest )
    # Get the indices of IR part 
    idx_ir = []
    irc = 0
    for iobs in d_diag['obs_type']:
        if d_diag['obs_type'][irc] == 'Radiance':
            idx_ir.append( irc )
        irc = irc + 1
    # Get the IR part
    IR_list = [[] for i in range(len(v_interest))]
    for var in  v_interest: 
        idx_var = v_interest.index(var)
        for it in idx_ir:
            if var == 'obs_type':
                IR_list[idx_var].append( d_diag[var][it])
            else:
                IR_list[idx_var].append( float(d_diag[var][it]))

    Diag_IR = {}
    print('Sanity check of reading process...Printing 222nd record...')
    for key in v_interest:
        idx_var = v_interest.index( key )
        Diag_IR[key] = IR_list[idx_var]
        print(key, ':', Diag_IR[key][222])




    return Diag_IR


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'

    # configuration
    Storm = 'HARVEY'
    #Exper_name = ['IR+MW-J_DA+J_WRF+J_init-SP-intel19',]
    filename = '/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/JerryRun/MW_THO/run/201708221200/enkf/d03/fort.10000'
    MW_SO = '/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/JerryRun/MW_THO/run/201708221200/enkf/d03/microwave_201708221200_so'
    IR_SO = '/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/JerryRun/MW_THO/run/201708221200/enkf/d03/radiance_201708221200_so'
    v_interest = ['obs_type','lat','lon','height','i_grid','j_grid','k_grid','Hroi','Vroi','obs','obs_error_variance','prior_mean','posterior_mean','prior_spread','posterior_spread']
    #label_mw_obs( filename, MW_SO, v_interest )
    Find_IR( filename, v_interest )







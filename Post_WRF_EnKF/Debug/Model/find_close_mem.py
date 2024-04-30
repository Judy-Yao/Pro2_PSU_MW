#!/usr/bin/env python3

# The purpose of this script is to, based on the minimum surface pressure, find the most similar member to the mean.

import os
import glob
import numpy as np
import scipy as sp
import scipy.ndimage
import netCDF4 as nc
from wrf import getvar
import math

def min_psfc(Hxa_file):
    ncdir = nc.Dataset( Hxa_file,'r')
    
    slp = getvar(ncdir, 'slp')
    min_slp = np.amin( slp )
    slp_smooth = sp.ndimage.filters.gaussian_filter(slp, [11, 11] )
    idx = np.nanargmin( slp_smooth )
    lat_x = ncdir.variables['XLAT'][:].flatten()[idx]
    lon_x = ncdir.variables['XLONG'][:].flatten()[idx]
    
    return lon_x,lat_x,int(min_slp)

def compare_mems(wrf_dir, DAtime):
    
    wrf_mean = wrf_dir + DAtime + '/wrf_enkf_output_d03_mean'
    #wrf_mean = wrf_dir + DAtime + '/enkf/d03/fort.90071'

    mean_lon, mean_lat, mean_psfc = min_psfc( wrf_mean )
    print('posterior mean:')
    print(mean_lon, mean_lat, mean_psfc)
    # Loop all members
    num_mem = 60
    mem = np.zeros([num_mem,3])
    #for idm in range(80011,80071,1):
    for imem in range(0,60,1):
        #Hxa_file = wrf_dir + DAtime + '/enkf/d03/fort.' + str(idm)
        #print( Hxa_file )
        #imem = idm - 80011
        Hxa_file = wrf_dir + DAtime + '/wrf_enkf_output_d03_' +'{:03}'.format(imem+1)
        mem[imem,:] = min_psfc( Hxa_file )      
        print(mem[imem,2])
    # Find the closest member
    print(mem[:,2] - mean_psfc)
    diff = np.amin(np.abs(mem[:,2] - mean_psfc))
    idx_closest = np.where(np.abs(mem[:,2] - mean_psfc) == diff)
    #print('Member closest to the mean is: ')
    #print('Member:', str(idx_closest[0] + 1 + 80010))
    print('Member:', str(idx_closest[0] + 1))
    print(mem[idx_closest, 0], mem[idx_closest, 1], mem[idx_closest, 2])


if __name__ == '__main__':

    Storm = 'MARIA'
    Exper_name = 'IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900/'
    wrf_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/fc/' 
    #wrf_dir = '/scratch/06191/tg854905/Pro3_MPcsst_CRTM/DORIAN/analyses/gts-tdr/'

    DAtimes = ['201709160600','201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000'] 

    for DAtime in DAtimes:
        print('At' + DAtime)
        compare_mems(wrf_dir, DAtime)

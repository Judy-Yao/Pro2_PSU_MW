#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import sys
import os # functions for interacting with the operating system
import numpy as np
import glob
import netCDF4 as nc
import math
import scipy as sp
import scipy.ndimage
import time
from wrf import getvar
# It might be possible that you are not able to conda install wrf-var with a pretty new python version
# Solution:
# 1. conda create -n $PYTHON34_ENV_NAME python=3.4 anaconda 
# 2. conda activate python=3.4 (use wrf-python in this python environment)

def onFC(wrfout):

    ncdir = nc.Dataset( wrfout )
    
    # domain
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    lat_min = np.amin(lat)
    lon_min = np.amin(lon)
    lat_max = np.amax(lat)
    lon_max = np.amax(lon)
    # Wind at 10 meters
    u10 = ncdir.variables['U10'][0,:,:]
    v10 = ncdir.variables['V10'][0,:,:]
    windspeed = (u10 ** 2 + v10 ** 2) ** 0.5
    # sea level pressure
    slp = getvar(ncdir, 'slp')
    min_slp = np.amin( slp )
    max_slp = np.amax( slp )
    slp_smooth = sp.ndimage.gaussian_filter(slp, [11,11])
    idx = np.nanargmin( slp_smooth )
    lat_minslp = ncdir.variables['XLAT'][:].flatten()[idx]
    lon_minslp = ncdir.variables['XLONG'][:].flatten()[idx]

def onFC_bdy(wrfbdy):

     ncdir = nc.Dataset( wrfbdy )
     u_bxe = ncdir.variables['U_BXE'][0,:,:,:]
     v_bye = ncdir.variables['V_BYE'][0,:,:,:]


if __name__ == '__main__':

    # Start timing
    start_time=time.process_time()

    Storm = sys.argv[1]
    Big_dir = sys.argv[2]
    Exper = sys.argv[3]
    print( '--------- Freshing directory: ', Big_dir, Storm, '/', Exper, '--------------' )

    dir_times = sorted(glob.glob( Big_dir + Storm + '/' + Exper +'/fc/*' ))
    for dir_time in dir_times:
        DAtime = os.path.split( dir_time )[1]
        print(DAtime, ' ...' )

        # Handle wrf_enkf files
        wrfouts = sorted(glob.glob( Big_dir + Storm + '/' + Exper + '/fc/' + DAtime + '/wrf_enkf*'  ))
        for wrfout in wrfouts:
            print(wrfout)
            onFC(wrfout)

        # Handle wrfinput files if there is any
        wrfouts = sorted(glob.glob( Big_dir + Storm + '/' + Exper + '/fc/' + DAtime + '/wrfinput*'  ))
        if wrfouts is not None:
            for wrfout in wrfouts:
                print(wrfout)
                onFC(wrfout)

        # Handle wrfbdy file
        wrfbdys = sorted(glob.glob( Big_dir + Storm + '/' + Exper + '/fc/' + DAtime + '/wrfbdy*'  ))
        if wrfbdys is not None:
            for wrfbdy in wrfbdys:
                print(wrfbdy)
                onFC_bdy(wrfbdy) 

    # End timing
    end_time = time.process_time() 
    print ('time needed: ', end_time-start_time, ' seconds')


























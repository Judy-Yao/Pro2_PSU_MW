
from numba import njit, prange
import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
from wrf import getvar
from scipy import interpolate
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from matplotlib import colors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import pickle
import random 

import Util_Vis
import Util_data as UD
import point2point_var_cov_corr_IR_x as p2p
#import Read_Obspace_IR as ROIR
import Diagnostics as Diag
#import matlab.engine

## Calculate the covariance*(N-1) between a Tb and a model variable over a 3D area (nlevel*lon*lat)
@njit(parallel=True)
def hor_obs_model( xb_ens,hxb_ens ):
    assert xb_ens.shape[1] == hxb_ens.shape[1]
    res = np.zeros( (xb_ens.shape[0],xb_ens.shape[2]),  )# levels, nobs
    res[:] = np.nan
    # for each level of variable: calculate corr at each grid point
    for n in prange( xb_ens.shape[2] ): # loop thru samples
        for i in range( xb_ens.shape[0] ): # loop thru model levels
            for m in range( xb_ens.shape[1] ): # loop thru ens
                res[i,n] += xb_ens[i,m,n] * hxb_ens[m,n]
    return res


## For each obs loc, find the nearest model grid point (good for mercator)
@njit(parallel=True)
def nearest_axis( obs,model ):

    res = np.zeros( (obs.shape[0]), )
    res[:] = np.nan
    for io in prange( obs.shape[0] ):
        for i in range( model.shape[0]-1 ):
            if model[i+1] < obs[io]:
                continue
            elif model[i] > obs[io]:
                continue
            else:
                if model[i] <= obs[io] and model[i+1] >= obs[io]:
                    if abs(model[i]-obs[io]) < abs(model[i+1]-obs[io]):
                        res[io] = i
                    else:
                        res[io] = i+1
    # Make sure every obs has a nearest model grid
    assert res.any() != np.nan
    return res


def Find_nearest_grid( lon_obs,lat_obs ): #lon_obs,lat_obs: lists

    print('------- Search for the nearest model grid for the obs ------')
    
    # Read model lon and lat
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    lon_x1d = ncdir.variables['XLONG'][0,0,:]
    lon_x1d = lon_x1d.filled(np.nan)
    lat_x1d = ncdir.variables['XLAT'][0,:,0]
    lat_x1d = lat_x1d.filled(np.nan)
    lon_x = ncdir.variables['XLONG'][0,:,:].flatten()
    lat_x = ncdir.variables['XLAT'][0,:,:].flatten()

    # Loop thru all obs and search each's left- and bottom-nearest model grid along x and y direction
    # returned is i,j in model domain
    print(lon_obs)
    Idx_i = nearest_axis( lon_obs,lon_x1d )
    Idx_j = nearest_axis( lat_obs,lat_x1d )

    # Transform to a 2d meshgrid and find the corresponding idx
    Idx_nearest = []
    for io in range(len(lon_obs)):
        Idx_nearest.append( int(Idx_j[io]*len(lon_x1d)+Idx_i[io]) )

    # check
    check_idx = random.randint(0, len(lon_obs))-1
    print('Checking the '+str(check_idx)+'th obs... lon: '+str(lon_obs[check_idx])+' lat: '+str(lat_obs[check_idx]))
    print('Its nearest model grid -- lon: '+str( lon_x[Idx_nearest[check_idx]] )+' lat: '+str( lat_x[Idx_nearest[check_idx]] ))

    if len(np.unique(Idx_nearest)) != len(lon_obs):
        warnings.warn('The nearest model grids might be repetitive!')
        print('Number of obs is '+str(len(lon_obs))+' while the number of the unique model location is '+str(len(np.unique(Idx_nearest))))

    return Idx_nearest

# Calculate the horizontal correlation centered at obs locations
def cal_hor_corr( DAtime, var_name, obs_type, d_obs, dim=None ):

    if dim == '3D':
        nLevel = 42
    elif dim == '2D':
        nLevel = 1

    # Find the flattened grid index near the obs
    lon_obs = np.array( d_obs[DAtime]['lon'] )
    lat_obs = np.array( d_obs[DAtime]['lat'] )
    idx_xb = Find_nearest_grid( lon_obs,lat_obs )

    # --- Find the model related statistics
    # Read ensemble perturbations of xb
    des_path = wrf_dir+ "xb_d03_3D_ensPert_" + DAtime + '_' + var_name +  '.pickle'
    with open( des_path,'rb' ) as f:
        xb_ens = pickle.load( f )
    xb_pert = xb_ens[:,:num_ens,:]
    print('Shape of xb_pert: '+ str(np.shape(xb_pert)))
    # Read ensemble standard deviation of xb
    des_path = wrf_dir+ "xb_d03_3D_ensStddev_" + DAtime + '_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_xb = pickle.load( f )
    print('Shape of stddev_xb: '+ str(np.shape(stddev_xb)))

    # --- Read the ensemble perturbation/spread of obs 
    if obs_type == 'slp':
        Hx_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
        #!!!!!! only works in the open ocean: approximate slp with PSFC!!!!!
        des_path = wrf_dir+ "xb_d03_2D_ensPert_" + DAtime + '_PSFC.pickle'
        with open( des_path,'rb' ) as f:
            tmp_hxb_pert = pickle.load( f )
        # Read ensemble standard deviation of xb
        des_path = wrf_dir+ "xb_d03_2D_ensStddev_" + DAtime + '_PSFC.pickle'
        with open( des_path,'rb' ) as f:
            tmp_stddev_hxb = pickle.load( f )

    # --- Calculate covariance
    cov_hor = np.zeros( (len(idx_xb),nLevel,xmax*ymax) )
    cov_hor[:] = np.nan

    for iobs in range(len(idx_xb)): 

        # Find the ensemble perturbation/spread of obs at obs locations
        if obs_type == 'slp':
            hxb_pert = tmp_hxb_pert[0,:num_ens,idx_xb[iobs]]
            print('Shape of hxb_pert: '+ str(np.shape(hxb_pert)))
            stddev_hxb = tmp_stddev_hxb[:,idx_xb[iobs]]
            print('Shape of stddev_hxb: '+ str(np.shape(stddev_hxb)))
        else:
            pass
        # Calculate covariance
        hro_cov[iobs,:,:] = hor_obs_model( xb_pert,hxb_pert )

    cov_hor = hor_obs_model( xb_pert,hxb_pert )





if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR'
    MP = 'THO'
    fort_v = ['obs_type','lat','lon','obs']
    sensor = 'abi_gr'

    # observation type 
    obs_assimilated = True
    if obs_assimilated:
        obs_type = 'slp' # Radiance
    
    # model variable
    model_v = [ 'QSNOW']
    
    # time
    start_time_str = '201709030000'
    end_time_str = '201709030000'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    #to_obs_res = False
    #ens_Interp_to_obs = False

    If_cal_pert_stddev = False
    If_cal_hor_corr = True
    If_save = True
    If_plot = False
    # -------------------------------------------------------
    Exper_name = UD.generate_one_name( Storm,DA,MP )

    # Identify DA times in the period of interest
    if not Consecutive_times:
        DAtimes = ['201709140000',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]


    # Read info of obs point of interest
    d_obs = {}
    if obs_assimilated:
        print('------------ Read obs info from enkf diagnostics fort.10000 --------------')
        for DAtime in DAtimes:
            # Read assimilated obs
            file_Diag = big_dir+Storm+'/'+Exper_name+'/run/'+DAtime+'/enkf/d03/fort.10000'
            if obs_type == 'slp':
                d_obs[DAtime] = Diag.Find_min_slp( file_Diag, fort_v )


    # Calculate ensemble perturbations and variances
    if If_cal_pert_stddev:
        print('------------ Calculate the ensemble perturbations --------------')
        for DAtime in DAtimes:
            wrf_dir =  big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            # Hxb
            if obs_type == 'slp':
                Hx_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
                #!!!!!! only works in the open ocean: approximate slp with PSFC!!!!!
                output_dir = wrf_dir+ "xb_d03_3D_ensPert_" + DAtime + '_PSFC.pickle'
                output_exists = os.path.exists( output_dir )
                if output_exists == False:
                    p2p.cal_pert_stddev_xb( DAtime, wrf_dir, 'PSFC', If_save, '2D')
            elif obs_type == 'Radiance':
                Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
                output_dir = Hx_dir+ "Hxb_ensStddev_obsRes_" + DAtime + '_' +  sensor + '.pickle'
                output_exists = os.path.exists( output_dir )
                if output_exists == False:
                    p2p.cal_pert_stddev_obsRes_Hxb( DAtime,sensor,Hx_dir,If_save,fort_v,wrf_dir)

            # Xb
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in v_interest:
                output_dir = wrf_dir+ "xb_d03_3D_ensPert_" + DAtime + '_' + var_name + '.pickle'
                output_exists = os.path.exists( output_dir )
                if output_exists == False:
                    p2p.cal_pert_stddev_xb( DAtime, wrf_dir, var_name, If_save, '3D')

    # Calculate horizontal correlations between the obs at the obs location and  model field across the specified domain
    if If_cal_hor_corr:
        print('------------ Calculate the horizonal correlation centered at the obs location --------------')
        for DAtime in DAtimes:
            #Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in model_v:

                print('Calculate '+var_name+'...')
                cal_hor_corr( DAtime, var_name, obs_type, d_obs, '3D')










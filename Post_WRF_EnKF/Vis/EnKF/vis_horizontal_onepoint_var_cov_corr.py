
from numba import njit, prange
import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
from wrf import getvar,interplevel
#from scipy import interpolate
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
import calculate_pert_stddev_x_IR as stat
import Read_Obspace_IR as ROIR
import Diagnostics as Diag
#import matlab.engine

def def_vardim( var_name ):
    if var_name == 'PSFC':
        return '2D'
    elif 'Q' in var_name:
        return '3D'

## Calculate the ensemble covariance*(N-1) between a obs and a model variable over a 3D area (nLevel*nLon*nLat)
# xb_ens: nLevel, num_ens, nLon*nLat
# hxb_ens: nLevel,num_ens,
@njit(parallel=True)
def hor_cov_obs_x( xb_ens,hxb_ens ):

    assert xb_ens.shape[1] == hxb_ens.shape[0]
    res = np.zeros( (xb_ens.shape[0],xb_ens.shape[2]), dtype=np.float64) # levels, nobs
    # for each level of variable: calculate corr at each grid point
    for n in prange( xb_ens.shape[2] ): # loop thru samples
        for i in range( xb_ens.shape[0] ): # loop thru model levels
            for m in range( xb_ens.shape[1] ): # loop thru ens
                res[i,n] = res[i,n] + xb_ens[i,m,n] * hxb_ens[m]
    return res

## Calculate the correlation between an obs and a 3D variable (nLevel*nLon*nLat)
# cov_xb_hxb: nLevel,nLon*nLat
# stddev_xb: nLevel, nLon*nLat
# stddev_hxb: 1
@njit(parallel=True)
def hor_corr_obs_x( cov_xb_hxb,stddev_xb,stddev_hxb ):

    #  use a threshold to avoid dividing by a number that is too small
    epsilon=1e-10
 
    assert cov_xb_hxb.shape[0] == stddev_xb.shape[0]
    assert cov_xb_hxb.shape[1] == stddev_xb.shape[1]
    res = np.zeros( (cov_xb_hxb.shape[0],cov_xb_hxb.shape[1]),dtype=np.float64)
    res[:] = np.nan
    for n in prange(cov_xb_hxb.shape[1]):
        for i in range(cov_xb_hxb.shape[0]):
            #print('cov_xb_hxb[i,n]',cov_xb_hxb[i,n])
            #print('stddev_xb[i,n]',stddev_xb[i,n])
            #print('stddev_hxb',stddev_hxb)

            if (abs(stddev_xb[i,n]) < epsilon):
                continue
            else:
                res[i,n] = cov_xb_hxb[i,n]/stddev_xb[i,n]

            #print('res[i,n]',res[i,n])

            if  (abs(stddev_hxb) < epsilon):
                continue
            else:
                res[i,n] = res[i,n]/stddev_hxb

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

# ---------------------------------------------------------------------------------------------------------------
#           Object: ensemble horizontal correlations of one obs and 2D/3D model variable
# ------------------------------------------------------------------------------------------------------------------

# Calculate the horizontal correlation centered at obs locations
def cal_hor_corr( DAtime, var_name, obs_type, d_obs, xdim=None ):

    if xdim == '3D':
        nLevel = 42
    elif xdim == '2D':
        nLevel = 1

    # Find the flattened grid index near the obs
    idx_obs_inX = []
    lon_obs = np.array( d_obs[DAtime]['lon'] )
    lat_obs = np.array( d_obs[DAtime]['lat'] )
    idx_obs_inX.append( Find_nearest_grid( lon_obs,lat_obs ) )

    # --- Find the model related statistics
    # Read ensemble perturbations of xb
    des_path = wrf_dir+ "xb_d03_3D_ensPert_" + DAtime + '_' + var_name +  '.pickle'
    with open( des_path,'rb' ) as f:
        xb_ens = pickle.load( f )
    xb_pert = xb_ens[:,:num_ens,:]
    print('Shape of xb_pert: '+ str(np.shape(xb_pert)))
    assert not np.isnan(xb_pert).any()
    # Read ensemble standard deviation of xb
    des_path = wrf_dir+ "xb_d03_3D_ensStddev_" + DAtime + '_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_xb = pickle.load( f )
    print('Shape of stddev_xb: '+ str(np.shape(stddev_xb)))
    assert not np.isnan(stddev_xb).any()

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
    print('Calculating the horizontal covariance between tho obs and ' + var_name + '.....' )
    cov_hor = np.zeros( (len(idx_obs_inX),nLevel,xmax*ymax) )
    #cov_hor[:] = np.nan

    start = time.perf_counter()
    for iobs in range(len(idx_obs_inX)): 

        # Find the ensemble perturbation/spread of obs at obs locations
        if obs_type == 'slp':
            hxb_pert = tmp_hxb_pert[0,:num_ens,idx_obs_inX[iobs]]
            print('Shape of hxb_pert: '+ str(np.shape(hxb_pert)))
            stddev_hxb = tmp_stddev_hxb[:,idx_obs_inX[iobs]]
            print('Shape of stddev_hxb: '+ str(np.shape(stddev_hxb)))
        else:
            pass
        # Calculate covariance
        assert not np.isnan(hxb_pert).any()
        assert not np.isnan(stddev_hxb).any()
        tmp_cov_hor = hor_cov_obs_x( xb_pert,hxb_pert[0,:] )
        cov_hor[iobs,:,:] = tmp_cov_hor / ( num_ens-1 )

    end = time.perf_counter()
    print("Elapsed (after compilation) of covariance calculation = {}s".format((end - start)))

    # Check if there are any NaN values using assert
    #assert not np.isnan(cov_hor).any()

    if If_save:
        des_path = wrf_dir+DAtime+"_d03_hroi_cov_Obs_" + obs_type +'_model_' +  var_name + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( cov_hor, f )
        f.close()
        print('Save '+des_path)

    # --- Calculate correlation
    print('Calculating the horizontal correlation between tho obs and ' + var_name + '.....' )
    hori_corr = np.zeros( (len(idx_obs_inX),nLevel,xmax*ymax) )
    #hori_corr[:] = np.nan

    start = time.perf_counter()
    for iobs in range(len(idx_obs_inX)):

        # Find the ensemble perturbation/spread of obs at obs locations
        if obs_type == 'slp':
            hxb_pert = tmp_hxb_pert[0,:num_ens,idx_obs_inX[iobs]]
            print('Shape of hxb_pert: '+ str(np.shape(hxb_pert)))
            stddev_hxb = tmp_stddev_hxb[:,idx_obs_inX[iobs]].item()
            print('Shape of stddev_hxb: '+ str(np.shape(stddev_hxb)))
        else:
            pass
        # Calculate covariance
        hori_corr[iobs,:,:] = hor_corr_obs_x( cov_hor[iobs,:,:],stddev_xb,stddev_hxb )
    end = time.perf_counter()
    print('Min of correlation: '+str(np.nanmin( hori_corr )))
    print('Max of correlation: '+str(np.nanmax( hori_corr )))

    print("Elapsed (after compilation) of correlation calculation = {}s".format((end - start)))    

    # sanity check
    assert  0 <= abs(hori_corr).any() and abs(hori_corr).any() <= 1

    # May save the correlation
    if If_save:
        des_path = wrf_dir+DAtime+"_d03_hroi_corr_Obs_" + obs_type +'_model_' +  var_name + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( hori_corr, f )
        f.close()
        print('Save '+des_path)

    return None


# ---------------------------------------------------------------------------------------------------------------
#           Object: ensemble horizontal correlations of one obs and 2D/3D model variable;Operation: Plot the snapshot
# ------------------------------------------------------------------------------------------------------------------
def HroiCorr_snapshot( DAtime,var_name,xdim=None,ver_coor=None):

    if xdim == '3D':
        nLevel = 42
    elif xdim == '2D':
        nLevel = 1

    # Read correlations between a Tb and a column of model var
    des_path = wrf_dir+DAtime+"_d03_hroi_corr_Obs_" + obs_type +'_model_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        hori_corr = pickle.load( f )
    print('Shape of hori_corr: '+ str(np.shape(hori_corr)))

    # Find the flattened grid index near the obs
    idx_obs_inX = []
    lon_obs = np.array( d_obs[DAtime]['lon'] )
    lat_obs = np.array( d_obs[DAtime]['lat'] )
    idx_obs_inX.append( Find_nearest_grid( lon_obs,lat_obs ) )

    # Plot
    # read model attributes
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    # read lon and lat
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()

    if xdim == '2D':
        if If_plot_corr_snapshot:
            for iobs in range(len(idx_obs_inX)):
                plot_2Dcorr_snapshot( xlat,xlon,hori_corr[iobs,:,:],lon_obs,lat_obs)
    else:
        if interp_H and not interp_P:
            # interpolate
            if If_plot_corr_snapshot:
                for iobs in range(len(idx_obs_inX)):
                    Interp_hori_corr = stat.vertical_interp( ncdir,hori_corr[iobs,:,:],ver_coor)
                    if If_plot_corr_snapshot:
                        #plot_3Dcorr_snapshot( xlat,xlon,hori_corr[iobs,25:31,:],H_of_interest )
                        plot_3Dcorr_snapshot( xlat,xlon,Interp_hori_corr,ver_coor,lon_obs,lat_obs)
        elif interp_P and not interp_H:

            if If_plot_corr_snapshot:
                plot_3Dcorr_snapshot( xlat,xlon,Interp_hori_corr,P_of_interest )
        else:
            pass

    return None

# Plot 2Dcorr per snapshot
def plot_2Dcorr_snapshot( lat,lon,corr,lon_obs,lat_obs ):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5,figsize=(6.5,6), dpi=300)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    min_corr = -0.6
    max_corr = 0.6
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)
    #cs = ax.flat[isub].scatter(lon,lat,5,Interp_corr[isub,:],cmap='RdBu_r',edgecolors='none',transform=ccrs.PlateCarree(),)
    cs = ax.scatter(lon,lat,5,corr,cmap='RdBu_r',edgecolors='none',vmin=min_corr,vmax=max_corr,transform=ccrs.PlateCarree(),)
    # Mark the observed location
    ax.scatter(lon_obs, lat_obs, c='red', s=5, marker='s', edgecolors='red', transform=ccrs.PlateCarree())
    #if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
    #    ax.flat[isub].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    color_ticks = np.linspace(min_corr, max_corr, 5, endpoint=True)
    cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, extend='both')
    cbar.set_ticks( color_ticks )
    cbar.ax.tick_params(labelsize=12)

    #subplot title
    font = {'size':12,}
    ax.set_title( 'Horizontal Corr: obs ' + obs_type +' & model ' +  var_name, font, fontweight='bold')

    #title for all
    title_name = Storm+': '+Exper_name
    fig.suptitle(title_name, fontsize=10, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    gl.top_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.right_labels = False

    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 12}


    # Save the figure
    save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/hori_Corr/'+DAtime+'_HroiCorr_obs_' + obs_type +'_model_' +  var_name + '.png'
    #save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/hori_Corr/Interp_H_corr_ms_'+DAtime+'_'+var_name+'_'+sensor+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


# Plot 3Dcorr per snapshot
def plot_3Dcorr_snapshot( lat,lon,Interp_corr,ver_coor,lon_obs,lat_obs):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # Read location from TCvitals
    #if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
    #    tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir,Storm, DAtime)

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(4, 5, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0.05, 'hspace':0.05}, linewidth=0.5, sharex='all', sharey='all',  figsize=(15,12), dpi=200)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    min_corr = -0.6
    max_corr = 0.6
    for isub in range(20):
        ax.flat[isub].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.flat[isub].coastlines(resolution='10m', color='black',linewidth=0.5)
        #cs = ax.flat[isub].scatter(lon,lat,5,Interp_corr[isub,:],cmap='RdBu_r',edgecolors='none',transform=ccrs.PlateCarree(),)
        cs = ax.flat[isub].scatter(lon,lat,5,Interp_corr[isub,:],cmap='RdBu_r',edgecolors='none',vmin=min_corr,vmax=max_corr,transform=ccrs.PlateCarree(),)
        # Mark the observed location
        ax.flat[isub].scatter(lon_obs, lat_obs, c='red', s=10, marker='s', edgecolors='red', transform=ccrs.PlateCarree())
        #if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        #    ax.flat[isub].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    color_ticks = np.linspace(min_corr, max_corr, 5, endpoint=True)
    cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, extend='both')
    cbar.set_ticks( color_ticks )
    cbar.ax.tick_params(labelsize=18)

    #subplot title
    font = {'size':15,}
    for isub in range(20):
        ax.flat[isub].set_title( str(ver_coor[isub])+' KM', font, fontweight='bold')

    #title for all
    title_name = Storm+': '+Exper_name+'\nHorizontal Corr: obs ' + obs_type +' & model ' +  var_name
    fig.suptitle(title_name, fontsize=20, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for i in range(4):
        for j in range(5):
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

            gl.top_labels = False
            gl.right_labels = False
            if j == 0:
                gl.left_labels = True
            else:
                gl.left_labels = False
            if i == 3:
                gl.bottom_labels = True
            else:
                gl.bottom_labels = False

            gl.ylocator = mticker.FixedLocator(lat_ticks)
            gl.xlocator = mticker.FixedLocator(lon_ticks)
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 10}
            gl.ylabel_style = {'size': 15}

    # Save the figure
    save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/hori_Corr/'+DAtime+'_Interp_H_HroiCorr_obs_' + obs_type +'_model_' +  var_name + '.png'
    #save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/hori_Corr/Interp_H_corr_ms_'+DAtime+'_'+var_name+'_'+sensor+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()

# ---------------------------------------------------------------------------------------------------------------
#           Object: ensemble horizontal covariance of one obs and 2D/3D model variable;Operation: Plot the snapshot
# ------------------------------------------------------------------------------------------------------------------

def HroiCov_snapshot( DAtime,var_name,xdim ):

    if xdim == '3D':
        nLevel = 42
    elif xdim == '2D':
        nLevel = 1

    # Read covelations between a Tb and a column of model var
    des_path = wrf_dir+DAtime+"_d03_hroi_cov_Obs_" + obs_type +'_model_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        hori_cov = pickle.load( f )
    print('Shape of hori_cov: '+ str(np.shape(hori_cov)))
    print( np.shape(hori_cov ))

    # Find the flattened grid index near the obs
    idx_obs_inX = []
    lon_obs = np.array( d_obs[DAtime]['lon'] )
    lat_obs = np.array( d_obs[DAtime]['lat'] )
    idx_obs_inX.append( Find_nearest_grid( lon_obs,lat_obs ) )

    # Plot
    # read model attributes
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    # read lon and lat
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()

    if xdim == '2D':
        if If_plot_cov_snapshot:
            for iobs in range(len(idx_obs_inX)):
                plot_2Dcov_snapshot( xlat,xlon,hori_cov[iobs,:,:],)
    else:
        if interp_H and not interp_P:
            # interpolate
            H_of_interest = [8,9,10,11,12,13]
            if If_plot_cov_snapshot:
                for iobs in range(len(idx_obs_inX)):
                    #Interp_hori_cov = #vertical_interp( ncdir,hori_cov[iobs,:,:],H_of_interest ) 
                    if If_plot_cov_snapshot:
                        # Count the number of NaN values
                        nan_count = np.isnan(hori_cov[iobs,0:6,:]).sum()
                        print("Number of NaN elements:", nan_count)
                        #print(np.nanmin(hori_cov[iobs,0:6,:]))
                        #print(np.nanmax(hori_cov[iobs,0:6,:]))
                        plot_3Dcov_snapshot( xlat,xlon,hori_cov[iobs,0:6,:],H_of_interest ) 
                        #plot_3Dcov_snapshot( xlat,xlon,Interp_hori_cov,H_of_interest )
        elif interp_P and not interp_H:

            if If_plot_cov_snapshot:
                plot_3Dcov_snapshot( xlat,xlon,Interp_hori_cov,P_of_interest ) 
        else:
            pass

    return None


# Plot 3D covariance per snapshot
def plot_3Dcov_snapshot( lat,lon,Interp_cov,ver_coor ):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # Read location from TCvitals
    #if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
    #    tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir,Storm, DAtime)

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(2, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(9.75,6.5), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    min_cov = -1
    max_cov = 1
    for isub in range(6):
        ax.flat[isub].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.flat[isub].coastlines(resolution='10m', color='black',linewidth=0.5)

        #cs = ax.flat[isub].scatter(lon,lat,5,Interp_cov[isub,:],cmap='RdBu_r',edgecolors='none',transform=ccrs.PlateCarree(),)
        cs = ax.flat[isub].scatter(lon,lat,5,Interp_cov[isub,:],cmap='RdBu_r',edgecolors='none',vmin=min_cov,vmax=max_cov,transform=ccrs.PlateCarree(),)
        # Mark the observed location
        ax.flat[isub].scatter(lon,lat,5,Interp_cov[isub,:],cmap='RdBu_r',edgecolors='none',vmin=min_cov,vmax=max_cov,transform=ccrs.PlateCarree(),)
        #    edgecolors='none', cmap='RdBu_r',transform=ccrs.PlateCarree())
        #if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        #    ax.flat[isub].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    color_ticks = np.linspace(min_cov, max_cov, 5, endpoint=True)
    cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, )
    cbar.set_ticks( color_ticks )
    cbar.ax.tick_params(labelsize=11)
    
    #subplot title
    font = {'size':15,}
    for isub in range(6):
        ax.flat[isub].set_title( str(ver_coor[isub])+' KM', font, fontweight='bold')

    #title for all
    title_name = Storm+': '+Exper_name+'(ver_cov of '+var_name+'&IR)'
    fig.suptitle(title_name, fontsize=10, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for j in range(6):
        gl = ax.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

        gl.top_labels = False
        gl.bottom_labels = True
        if j==0 or j==3:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False

        if j==3 or j==4 or j==5:
            gl.bottom_labels = True
            gl.top_labels = False
        else:
            gl.bottom_labels = False
            gl.top_labels = False

        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10}
        gl.ylabel_style = {'size': 12}

    # Save the figure
    save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/hori_Corr/Interp_H_HroiCov_obs_' + obs_type +'_model_' +  var_name + '.png'
    #save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/hori_Cov/Interp_H_cov_ms_'+DAtime+'_'+var_name+'_'+sensor+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/' #'/expanse/lustre/scratch/zuy121/temp_project/Pro2_PSU_MW/' #'/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/' #'/expanse/lustre/projects/pen116/zuy121/Pro2_PSU_MW/'  #'/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'CONV'
    MP = 'THO'
    fort_v = ['obs_type','lat','lon','obs']
    sensor = 'abi_gr'

    # observation type 
    obs_assimilated = True
    if obs_assimilated:
        obs_type = 'slp' # Radiance
    
    # model variable
    model_v = [ 'QSNOW',]#'QSNOW','QCLOUD','QRAIN','QICE','QGRAUP']
    
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

    # vertical interpolation if needed
    interp_P = False
    P_range = np.arange( 995,49,-20 )
    interp_H = True
    H_range = list(np.arange(1,21,1))

    If_cal_pert_stddev = False
    If_cal_hor_corr = False
    If_save = True

    If_plot_corr_snapshot = True
    If_plot_cov_snapshot = False

    # ensemble spread in vis_point2point_stat_IR_toColumnModel.py 

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
                output_dir = wrf_dir+ "xb_d03_2D_ensPert_" + DAtime + '_PSFC.pickle'
                output_exists = os.path.exists( output_dir )
                if output_exists == False:
                    stat.cal_pert_stddev_xb( DAtime, wrf_dir, 'PSFC', If_save, '2D')
            elif obs_type == 'Radiance':
                Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
                output_dir = Hx_dir+ "Hxb_ensStddev_obsRes_" + DAtime + '_' +  sensor + '.pickle'
                output_exists = os.path.exists( output_dir )
                if output_exists == False:
                    stat.cal_pert_stddev_obsRes_Hxb( DAtime,sensor,Hx_dir,If_save,fort_v,wrf_dir)

            # Xb
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in model_v:
                var_dim = def_vardim( var_name )
                output_dir = wrf_dir+ 'xb_d03_'+var_dim+'_ensPert_' + DAtime + '_' + var_name + '.pickle'
                output_exists = os.path.exists( output_dir )
                if output_exists == False:
                    stat.cal_pert_stddev_xb( DAtime, wrf_dir, var_name, If_save, '3D')

    # Calculate horizontal correlations between the obs at the obs location and  model field across the specified domain
    if If_cal_hor_corr:
        print('------------ Calculate the horizonal correlation centered at the obs location --------------')
        for DAtime in DAtimes:
            #Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in model_v:

                print('Calculate '+var_name+'...')
                var_dim = def_vardim( var_name )
                cal_hor_corr( DAtime, var_name, obs_type, d_obs, var_dim )

    # Plot the ensemble spread of xb per snapshot
    if If_plot_stddev_xb_snapshot:
        print('------------ Plot the ensemble spread of xb --------------')
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Ens_stddev_xb/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

        for DAtime in DAtimes:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            #Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            print('At '+DAtime)
            for var_name in model_v:
                print('Plot '+var_name+'...')
                var_dim = def_vardim( var_name )
                stddev_xb_snapshot( wrf_dir,DAtime,var_name,var_dim )







    # Plot the horizontal covariance per snapshot
    if If_plot_cov_snapshot:
        print('------------ Plot the horizontal correlation --------------')
        for DAtime in DAtimes:
            #Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            print('At '+DAtime)
            for var_name in model_v:

                print('Plot horizontal covariance: '+var_name+'...')
                HroiCov_snapshot( DAtime,var_name,'3D' )

    # Plot the horizontal correlations per snapshot
    if If_plot_corr_snapshot:
        print('------------ Plot the horizontal correlation --------------')
        for DAtime in DAtimes:
            #Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            print('At '+DAtime)
            for var_name in model_v:

                print('Plot horizontal correlation: '+var_name+'...')
                if interp_H and not interp_P:
                    HroiCorr_snapshot( DAtime,var_name,'3D',H_range)







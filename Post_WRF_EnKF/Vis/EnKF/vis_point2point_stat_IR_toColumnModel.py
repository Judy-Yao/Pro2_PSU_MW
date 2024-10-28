
from numba import njit, prange
import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
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
import warnings
from wrf import getvar, interplevel
import random

import Util_Vis
import Util_data as UD
import Read_Obspace_IR as ROIR
import Diagnostics as Diag

# ------------------------------------------------------------------------------------------------------
#           Operation: Perform Operations
# ------------------------------------------------------------------------------------------------------

## Calculate the covariance*(N-1) between a Tb and a column of model variale over a 2D area (lon*lat)
@njit(parallel=True)
def cross_Tb_toCol( xb_ens,hxb_ens ):
    assert xb_ens.shape[1] == hxb_ens.shape[0]
    res = np.zeros( (xb_ens.shape[0],xb_ens.shape[2]),  )# levels, nobs
    # for each level of variable: calculate corr at each grid point
    for n in prange( xb_ens.shape[2] ): # loop thru samples
        for i in range( xb_ens.shape[0] ): # loop thru model levels
            for m in range( xb_ens.shape[1] ): # loop thru ens
                res[i,n] += xb_ens[i,m,n] * hxb_ens[m,n]
    return res

## Calculate the correlation*(N-1) between a Tb and a column of model variale over a 2D area (lon*lat)
@njit(parallel=True)
def corr_Tb_toCol( cov_xb_hxb,stddev_xb,stddev_hxb ):
    assert cov_xb_hxb.shape[0] == stddev_xb.shape[0]
    assert cov_xb_hxb.shape[1] == stddev_hxb.shape[0]
    res = np.zeros( (cov_xb_hxb.shape[0],cov_xb_hxb.shape[1]), )
    for n in prange(cov_xb_hxb.shape[1]):
        for i in range(cov_xb_hxb.shape[0]):
                res[i,n] = cov_xb_hxb[i,n]/stddev_xb[i,n]
                res[i,n] = res[i,n]/stddev_hxb[n]
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


# vertical interpolation
def vertical_interp( ncdir,array,levels ):

    if interp_H and not interp_P:
        #Interp_corr_xb_hxb = np.zeros( [len(H_of_interest),len(idx_xb)] )
        z = getvar(ncdir, 'z', units='km')
        Interp_corr_xb_hxb = np.zeros( (len(levels),z.shape[1],z.shape[2]) )
        array =  array.reshape( (z.shape) )
        start_time=time.process_time()
        for ih in levels:
            Interp_corr_xb_hxb[levels.index(ih),:,:] = interplevel(array, z, ih)
        end_time = time.process_time()
        print ('time needed for the interpolation: ', end_time-start_time, ' seconds')
        print('Min of correlation: '+str(np.amin( Interp_corr_xb_hxb )))
        print('Max of correlation: '+str(np.amax( Interp_corr_xb_hxb )))
        return Interp_corr_xb_hxb
    elif interp_P and not interp_H:
        # pressure levels
        PB = ncdir.variables['PB'][0,:,:,:]
        P = ncdir.variables['P'][0,:,:,:]
        P_hpa_all = (PB + P)/100 # 0 dimension: bottom to top
        P_hpa_all = P_hpa_all.reshape(nLevel,xmax*ymax)
        P_hpa = P_hpa_all[:,idx_xb]
        # Calculate the corr at specified levels
        start_time=time.process_time()
        array_P = np.zeros( (len(P_of_interest),array.shape[1]),  )
        for im in range( array.shape[1] ):
            f_interp = interpolate.interp1d( P_hpa[:,im], array[:,im])
            Interp_corr_xb_hxb[:,im] = f_interp( P_of_interest )
        #corr_colxb_hxb_cloud_P = corr_colxb_hxb_cloud[:3,:] # test....
        end_time = time.process_time()
        print ('time needed for the interpolation: ', end_time-start_time, ' seconds')
        print('Min of correlation: '+str(np.amin( Interp_corr_xb_hxb )))
        print('Max of correlation: '+str(np.amax( Interp_corr_xb_hxb )))

    else:
        pass



# ------------------------------------------------------------------------------------------------------
#           Object: ensemble correlations of columns of Xb and Hxb in 2D; Operation: Calculation
# ------------------------------------------------------------------------------------------------------
def Find_nearest_col( Hx_dir, wrf_dir, DAtime, sensor):

    start_time=time.process_time()
    print('------- Search for the nearest model grid for each obs ------')
    # Read the obs locations
    lat_obs = []
    lon_obs = []
    mean_hxb = Hx_dir + "mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    with open(mean_hxb) as f:
        next(f)
        all_lines = f.readlines()
    for line in all_lines:
        split_line = line.split()
        lat_obs.append( float(split_line[0]) )
        lon_obs.append( float(split_line[1]) ) 
    lat_obs = np.array( lat_obs )
    lon_obs = np.array( lon_obs )

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
    Idx_i = nearest_axis( lon_obs,lon_x1d )
    Idx_j = nearest_axis( lat_obs,lat_x1d )

    # Transform to a 2d meshgrid and find the corresponding idx
    Idx_nearest = []
    for io in range(len(lon_obs)):
        Idx_nearest.append( int(Idx_j[io]*len(lon_x1d)+Idx_i[io]) )

    # check 
    check_idx = random.randint(0, len(lon_obs))
    print('Checking the '+str(check_idx)+'th obs... lon: '+str(lon_obs[check_idx])+' lat: '+str(lat_obs[check_idx]))
    print('Its nearest model grid -- lon: '+str( lon_x[Idx_nearest[check_idx]] )+' lat: '+str( lat_x[Idx_nearest[check_idx]] ))

    if len(np.unique(Idx_nearest)) != len(lon_obs):
        warnings.warn('The nearest model grids might be repetitive!')
        #raise ValueError('The nearest model grids might be repetitive!')
        print('Number of obs is '+str(len(lon_obs))+' while the number of the unique model location is '+str(len(np.unique(Idx_nearest))))

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    return Idx_nearest    


def cal_2Dcorr_IR_ColVar( DAtime, var_name):

    # Read ensemble perturbations of xb
    des_path = wrf_dir+ "xb_d03_3D_ensPert_" + DAtime + '_' + var_name +  '.pickle'
    with open( des_path,'rb' ) as f:
        xb_ens = pickle.load( f )
    print('Shape of xb_ens: '+ str(np.shape(xb_ens)))
    # Read ensemble standard deviation of xb
    des_path = wrf_dir+ "xb_d03_3D_ensStddev_" + DAtime + '_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_xb = pickle.load( f )
    print('Shape of stddev_xb: '+ str(np.shape(stddev_xb)))
    # Read ensemble perturbations of Hxb 
    if to_obs_res:
        des_path = Hx_dir+ "Hxb_ensPert_obsRes_" + DAtime + '_' +  sensor + '.pickle'
    else:
        des_path = Hx_dir+ "Hxb_ensPert_modelRes_" + DAtime + '_' +  sensor + '.pickle'
    with open( des_path,'rb' ) as f:
        hxb_ens = pickle.load( f )
    print('Shape of hxb_ens: '+ str(np.shape(hxb_ens)))
    # Read ensemble stand deviation of Hxb
    if to_obs_res:
        des_path = Hx_dir+ "Hxb_ensStddev_obsRes_" + DAtime + '_' +  sensor + '.pickle'
    else:
        des_path = Hx_dir+ "Hxb_ensStddev_modelRes_" + DAtime + '_' +  sensor + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_hxb = pickle.load( f )
    print('Shape of stddev_hxb: '+ str(np.shape(stddev_hxb)))

    # Find the location of model grid of interest
    if to_obs_res: # nearest for each obs
        idx_xb = Find_nearest_col( Hx_dir, wrf_dir, DAtime, sensor)
    else: # every model grid point
        idx_xb = np.arange(xmax*ymax) 

    # Calculate the covariance between Xb and Hxb (a column of var and a Tb)
    print('Calculating the covariance between Tbs and ' + var_name + '(a column of var and a Tb)......' )
    start = time.perf_counter()
    cov_xb_hxb = cross_Tb_toCol( xb_ens[:,:num_ens,idx_xb],hxb_ens[:num_ens,:] )
    end = time.perf_counter()
    print("Elapsed (after compilation) of covariance calculation = {}s".format((end - start)))
    cov_xb_hxb = cov_xb_hxb / ( num_ens-1 )

    # Check if there are any NaN values using assert
    assert not np.isnan(cov_xb_hxb).any()

    # Calculate the correlation between Xb and Hxb
    print('Calculating the correlation between Tbs and ' + var_name + '......' )
    start = time.perf_counter()
    corr_xb_hxb = corr_Tb_toCol( cov_xb_hxb,stddev_xb[:,idx_xb],stddev_hxb )
    end = time.perf_counter()
    print("Elapsed (after compilation) = {}s".format((end - start)))

    # sanity check
    assert  0 <= abs(corr_xb_hxb).any() and abs(corr_xb_hxb).any() <= 1

    # May save the correlation
    if If_save:
        if to_obs_res:
            des_path = wrf_dir+ "d03_2Dcorr_obsRes_Hxb_" + DAtime + '_Column_' +  var_name + '.pickle'
        else:
            des_path = wrf_dir+ "d03_2Dcorr_modelRes_Hxb_" + DAtime + '_Column_' +  var_name + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( corr_xb_hxb, f )
        f.close()
        print('Save '+des_path)

    return None



# ------------------------------------------------------------------------------------------------------
#           Object: ensemble correlations of columns of Xb and Hxb in 2D; Operation: Plot the snapshot
# ------------------------------------------------------------------------------------------------------

# Plot 3Dcorr per snapshot
def plot_3Dcorr_snapshot( lat,lon,Interp_corr,ver_coor ):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # Read Tbs of Hxb
    Tb_file = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/' + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    d_all = ROIR.read_Tb_obsRes(Tb_file, sensor )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir,Storm, DAtime)

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(2, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(9.75,6.5), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    min_corr = -1
    max_corr = 1
    for isub in range(6):
        ax.flat[isub].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.flat[isub].coastlines(resolution='10m', color='black',linewidth=0.5)
        #cs = ax.flat[isub].scatter(lon,lat,5,Interp_corr[isub,:],cmap='RdBu_r',edgecolors='none',transform=ccrs.PlateCarree(),)
        cs = ax.flat[isub].scatter(lon,lat,5,Interp_corr[isub,:],cmap='RdBu_r',edgecolors='none',vmin=min_corr,vmax=max_corr,transform=ccrs.PlateCarree(),)
        #    edgecolors='none', cmap='RdBu_r',transform=ccrs.PlateCarree())
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
            ax.flat[isub].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    color_ticks = np.linspace(min_corr, max_corr, 5, endpoint=True)
    cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, )
    cbar.set_ticks( color_ticks )
    cbar.ax.tick_params(labelsize=11)

    #subplot title
    font = {'size':15,}
    for isub in range(6):
        ax.flat[isub].set_title( str(ver_coor[isub])+' KM', font, fontweight='bold')

    #title for all
    title_name = Storm+': '+Exper_name+'(ver_corr of '+var_name+'&IR)'
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
    if to_obs_res and interp_H:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/Interp_H_corr_os_'+DAtime+'_'+var_name+'_'+sensor+'.png'
    else:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/Interp_H_corr_ms_'+DAtime+'_'+var_name+'_'+sensor+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()
    return None

# Plot corr per snapshot
def plot_2Dcorr_snapshot( lat,lon,corr_2d ):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # Read Tbs of Hxb
    Tb_file = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/' + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    d_all = ROIR.read_Tb_obsRes(Tb_file, sensor )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir,Storm, DAtime)

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5,figsize=(6.5,6), dpi=300)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    min_corr = -0.5
    max_corr = 0.5
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax.scatter(lon,lat,15,corr_2d,cmap='RdBu_r',edgecolors='none',vmin=min_corr,vmax=max_corr,transform=ccrs.PlateCarree(),)
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax.scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    color_ticks = np.linspace(min_corr, max_corr, 5, endpoint=True)
    cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, )
    cbar.set_ticks( color_ticks )
    cbar.ax.tick_params(labelsize=11)

    #subplot title
    font = {'size':15,}
    ax.set_title( 'Point-to-point Corr: '+var_name+'& IR Tbs', font, fontweight='bold')

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
    if to_obs_res:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/corr_os_'+DAtime+'_'+var_name+'_'+sensor+'.png'
    else:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/corr_ms_'+DAtime+'_'+var_name+'_'+sensor+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()
    return None


def corr_snapshot( DAtime,var_name,var_dim ):

    if var_dim == '3D':
        nLevel = 42
    elif var_dim == '2D':
        nLevel = 1

    # Read correlations between a Tb and a column of model var
    if to_obs_res:
        des_path = wrf_dir+ "d03_2Dcorr_obsRes_Hxb_" + DAtime + '_Column_' +  var_name + '.pickle'
    else:
        des_path = wrf_dir+ "d03_2Dcorr_modelRes_Hxb_" + DAtime + '_Column_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        corr_colxb_hxb = pickle.load( f )
    print('Shape of corr_colxb_hxb: '+ str(np.shape(corr_colxb_hxb)))

    # Find the location of model grid of interest
    if to_obs_res: # nearest for each obs
        idx_xb = Find_nearest_col( Hx_dir, wrf_dir, DAtime, sensor)
    else: # every model grid point
        idx_xb = np.arange(xmax*ymax)    
    print('Shape of idx_xb: '+ str(np.shape(idx_xb)))

    # Plot
    # read model attributes
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    # read lon and lat
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()

    if var_dim == '2D':
        if If_plot_corr_snapshot:
            plot_2Dcorr_snapshot( xlat[idx_xb],xlon[idx_xb],corr_colxb_hxb,)
    else:
        if interp_H and not interp_P:
            # interpolate
            H_of_interest = [1,2,3,4,5,6]
            Interp_corr_xb_hxb = vertical_interp( ncdir,corr_colxb_hxb )
            if If_plot_corr_snapshot:
                plot_3Dcorr_snapshot( xlat[idx_xb],xlon[idx_xb],Interp_corr_xb_hxb,H_of_interest )
        elif interp_P and not interp_H:
             
            if If_plot_corr_snapshot: 
                plot_3Dcorr_snapshot( xlat[idx_xb],xlon[idx_xb],Interp_cov_xb_hxb,P_of_interest )
        else:
            pass
        
    return None

# ------------------------------------------------------------------------------------------------------
#           Object: ensemble correlations of columns of Xb and Hxb in 2D; Operation: Plot the evolution
# ------------------------------------------------------------------------------------------------------

# Plot
def plot_var_incre_timeseries( ave_corr_overT,geoHkm=None ):

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)
    # Set up coordinates
    if not interp_P:
        xv = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]
        #y_bottom = np.arange(0,10,1)
        y_bottom = np.arange(0,5,0.5)
        #y_middle = np.arange(10,15,0.5)
        y_top = np.arange(5,31,1)
        #y_top = np.arange(15,31,1)
        y_range = np.concatenate( (y_bottom,y_top),axis=0 )
        #y_range = np.concatenate( (y_bottom,y_middle,y_top),axis=0 ) # control the scale
        y_axis_rg = range(len(y_range))
        f_interp = interpolate.interp1d( y_range, y_axis_rg)
        yv = f_interp( geoHkm )
        xcoor, ycoor = np.meshgrid( xv, yv )
    else:
        # Create a coordinate matrix from coordinate vectors
        xv = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]#range( np.shape(ave_norm_overT)[0] )
        yv = range( np.shape(ave_corr_overT)[1])
        xcoor, ycoor = np.meshgrid( xv, yv )

    max_corr = 0.3
    min_corr = -0.3
    bounds = np.linspace(min_corr, max_corr, 13)
    cmap = 'bwr'
    incre_contourf = ax.contourf( xcoor, ycoor, np.transpose( ave_corr_overT ), cmap=cmap, vmin=min_corr, vmax=max_corr,levels=bounds,extend='both')
    # Add color bar below the plot
    color_bar = fig.colorbar(incre_contourf,orientation = 'horizontal',pad=0.15)
    color_bar.ax.tick_params(labelsize=12)
    color_bar.ax.set_xlabel('Domain-mean Vertical Corr',fontsize=15)
    
    # set X label
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45)
    # set Y label
    if not interp_P:
        ylabel_like = [0.0,1.0,2.0,3.0,4.0,5.0,10.0,15.0,20.0]
        #ylabel_like = [0.0,1.0,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.0]
        #ylabel_like = [0.0,5.0,10.0,11.0,12.0,13.0,14.0,15.0,20.0,25.0,30.0]
        yticks = []
        list_y_range = list(y_range)
        for it in ylabel_like:
            yticks.append( list_y_range.index(it) )
        ax.set_yticks( yticks )
        ax.set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
        ax.set_ylabel('Height (KM)',fontsize=15)
        ax.set_ylim(ymin=0,ymax=25) # cut off data above 25km
    else:
        ax.set_yticks( yv[::10] )
        ax.set_yticklabels( [str(it) for it in P_interest[::10]],fontsize=15 )
        ax.set_ylabel('Pressure (hPa)',fontsize=15)

    # Set title
    if 'Q' in var_name:
        title_name = 'Vertical Correlation: '+var_name+' - IR Ch'+ch_list[0]
    ax.set_title( title_name,fontweight="bold",fontsize='12' )
    fig.suptitle(Storm+': '+Exper_name, fontsize=10, fontweight='bold')   

 
    # Save the figure
    if not interp_P:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/IR/ML_vCorr_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    else:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/IR/Interp_Vcorr_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()

# Plot the time-averaged area-mean correlation
def plot_Tmean_VD_corr( ave_corr_overT, ver_coor ):

    # Calculate the time-averaged area-averaged corr
    if limit:
        ave_corr = np.mean( ave_corr_overT, axis=1 )
    else:
        ave_corr = np.mean( ave_corr_overT, axis=0 )   

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)

    ## x axis: correlation value
    x_range = np.arange(-1,1.05,0.05)
    x_axis_rg = range(len(x_range))
    fx_interp = interpolate.interp1d( x_range, x_axis_rg)
    ## y axis: vertical coordinate
    if interp_H and not interp_P:
        y_range = H_range
        #y_bottom = np.arange(0,10,1)
        #y_middle = np.arange(10,15,0.5)
        #y_top = np.arange(15,31,1)
        #y_range = np.concatenate( (y_bottom,y_middle,y_top),axis=0 ) # control the scale
    else:
        y_range = P_range
    y_axis_rg = range(len(y_range))
    fy_interp = interpolate.interp1d( y_range, y_axis_rg)
    loc_iny = fy_interp( ver_coor)

    # Plot
    if limit:
        loc_inx = fx_interp( ave_corr[0,:] )
        ax.plot( loc_inx,loc_iny,linewidth=5,color='#A9A9A9',label='all-sky')
        loc_inx = fx_interp( ave_corr[1,:] )
        ax.plot( loc_inx,loc_iny,linewidth=5,color='#097969',label='cloud_goodX') 
        loc_inx = fx_interp( ave_corr[2,:] )
        ax.plot( loc_inx,loc_iny,linewidth=5,color='#BF565A',label='cloud_badX')
    else:
        loc_inx = fx_interp( ave_corr )
        ax.plot( loc_inx,loc_iny,linewidth=5,color='#9A96CC')
    nocorr_loc = fx_interp(0)
    ax.axvline(x=nocorr_loc,color='black',linestyle='-',linewidth=3)

    # Legend
    ax.legend(loc='upper right',fontsize='16')

    # set X label
    if 'Q' in var_name:
        xlabel_like = [-0.8,-0.4,0.0,0.2,0.4]
    else:
        xlabel_like = [-0.4,-0.2,0.0,0.2,0.4]
    xticks = []
    for it in xlabel_like:
        xticks.append( fx_interp( it ) )
    ax.set_xticks( xticks )
    ax.set_xticklabels( [str(it) for it in xlabel_like],fontsize=20 )
    # small ticks
    if 'Q' in var_name:
        xlabel_like = [-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4]
    else:
        xlabel_like = [-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4]
    xticks_minor = []
    for it in xlabel_like:
        xticks_minor.append( fx_interp( it ) )
    ax.xaxis.set_minor_locator(mticker.FixedLocator(xticks_minor))
    ax.set_xlabel('Correlation',fontsize=20)
    #ax.set_xlim(xmin=-0.5,xmax=0.5) 
    # grids
    ax.grid(which='both',axis='x',linewidth=1,alpha=0.5)
    # set Y label
    if interp_H and not interp_P:
        ylabel_like = [0.0,5.0,10.0,15.0,20.0,25.0]
        yticks = []
        for it in ylabel_like:
            yticks.append( fy_interp( it ) ) 
        ax.set_yticks( yticks )
        ax.set_yticklabels( [str(it) for it in ylabel_like],fontsize=20 )
        # small ticks
        ylabel_like = np.arange(0,25.5,1)
        yticks_minor = []
        for it in ylabel_like:
            yticks_minor.append( fy_interp( it ) )
        ax.yaxis.set_minor_locator(mticker.FixedLocator(yticks_minor)) 
        ax.set_ylabel('Height (KM)',fontsize=20)
        ax.set_ylim(ymin=0,ymax=25) # cut off data above 25km
    else:
        pass
    # grids
    ax.grid(which='both',linewidth=1,alpha=0.7,axis='y')

    # Set title
    title_name = 'Tmean_Dmean_Vertical Corr: '+var_name+'&IR Ch'+ch_list[0] +'\n Cloud:obs<220K;goodX:abs(H(Xb)-obs)<=3K;badX:H(Xb)-obs>=3K'
    ax.set_title( title_name,fontweight="bold",fontsize='12' )
    fig.suptitle(Storm+': '+Exper_name, fontsize=10, fontweight='bold')

    # Save the figure
    if not interp_P:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/IR/ML_Tmean_vCorr_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    else:        
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/IR/Interp_Tmean_vCorr_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


# Calculate domain-averaged correlation
def Domain_ave_corr( var_name ):

    nLevel = 42

    # Construct array
    vercoor_overT = np.zeros( [len(DAtimes),nLevel] )
    if limit:
        ave_corr_overT = np.zeros( [num_limit,len(DAtimes),nLevel] ) 
    else:
        ave_corr_overT = np.zeros( [len(DAtimes),nLevel] )
    ave_corr_overT[:] = np.nan

    # Load correlation over time and interpolate if necessary
    for DAtime in DAtimes:

        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
        t_idx = DAtimes.index( DAtime )

        # Read correlations between a Tb and a column of model var
        if to_obs_res:
            des_path = wrf_dir+ "d03_2Dcorr_obsRes_Hxb_" + DAtime + '_Column_' +  var_name + '.pickle'
        else:
            des_path = wrf_dir+ "d03_2Dcorr_modelRes_Hxb_" + DAtime + '_Column_' +  var_name + '.pickle' 
        with open( des_path,'rb' ) as f:
            corr_colxb_hxb = pickle.load( f )
        print('Shape of corr_colxb_hxb: '+ str(np.shape(corr_colxb_hxb)))

        # Read domain-mean vertical coordinate 
        mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
        ncdir = nc.Dataset( mean_xb, 'r')
        if interp_H:
            PHB = ncdir.variables['PHB'][0,:,:,:]
            PH = ncdir.variables['PH'][0,:,:,:]
            geoHkm = (PHB+PH)/9.8/1000 # in km
            geoHkm = geoHkm.reshape( geoHkm.shape[0],-1)
            geoHkm_Dmean = np.mean( geoHkm, axis=1 )
            geoHkm_half_eta = (geoHkm_Dmean[:-1]+geoHkm_Dmean[1:])/2
            geoHkm_half_eta = np.ma.getdata(geoHkm_half_eta)
            ver_coor = geoHkm_half_eta
        else:
            PB = ncdir.variables['PB'][0,:,:,:]
            P = ncdir.variables['P'][0,:,:,:]
            P_hpa = ((PB + P)/100)
            P_hpa = P_hpa.reshape( P_hpa.shape[0],-1) # 0 dimension: bottom to top
            P_hpa = np.mean( P_hpa, axis=1 )
            ver_coor = P_hpa
        vercoor_overT[t_idx,:] = ver_coor

        # Impose limitations
        if limit:
            d_idx = {}
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            Tb_file = Hx_dir + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
            d_all = ROIR.read_Tb_obsRes(Tb_file, sensor )
            condi1 = d_all['Yo_obs'] <= 220
            idx_1 = np.where( condi1 )[0]
            condi2 = abs(d_all['meanYb_obs'] - d_all['Yo_obs']) <= 3 #d_all['Yo_obs'] <= 210
            idx_2 = np.where( condi2 )[0]
            idx_x1 = list(set(idx_1)&set(idx_2))
            d_idx['limit0'] = range(corr_colxb_hxb.shape[1]) # all points
            d_idx['limit1'] = idx_x1
            condi3 = d_all['meanYb_obs'] - d_all['Yo_obs'] > 3
            idx_3 = np.where( condi3 )[0]
            idx_x2 = list(set(idx_1)&set(idx_3))
            d_idx['limit2'] = idx_x2
        else:
            idx_x = range(corr_colxb_hxb.shape[1])      
            d_idx = idx_x
 
        # Calculate the domain-mean or time-mean 
        ave_vercoor_overT = np.mean( vercoor_overT, axis=0 )
        if limit:
            for i in range(num_limit):
                ave_corr_overT[i,t_idx,:] = np.nanmean( corr_colxb_hxb[:,d_idx['limit'+str(i)]], axis=1 )
        else: 
            ave_corr_overT[t_idx,:] = np.nanmean( corr_colxb_hxb[:,idx_x], axis=1 )

    # Plot the time-averaged and time-evolution of vertical distribution of corr
    plot_Tmean_VD_corr( ave_corr_overT, ave_vercoor_overT )
    #plot_var_incre_timeseries( ave_corr_overT,ver_coor )

    return None

# ------------------------------------------------------------------------------------------------------
#    Object: ensemble spread of columns of hxb; Operation: Plot the snapshot
# ------------------------------------------------------------------------------------------------------

def plot_stddev_hxb_snapshot( DAtime ):

    # Read ensemble stand deviation of Hxb
    if to_obs_res:
        des_path = Hx_dir+ "Hxb_ensStddev_obsRes_" + DAtime + '_' +  sensor + '.pickle'
    else:
        des_path = Hx_dir+ "Hxb_ensStddev_modelRes_" + DAtime + '_' +  sensor + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_hxb = pickle.load( f )
    print('Shape of stddev_hxb: '+ str(np.shape(stddev_hxb)))

    # Find the location of model grid of interest
    if to_obs_res: # nearest for each obs
        idx_xb = Find_nearest_col( Hx_dir, wrf_dir, DAtime, sensor)
    else: # every model grid point
        idx_xb = np.arange(xmax*ymax)

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_input_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )
    ncdir = nc.Dataset( wrf_file, 'r')
    # read lon and lat
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir,Storm, DAtime)

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5,figsize=(6.5,6), dpi=300)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    #min_corr = -0.5
    #max_corr = 0.5
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax.scatter(xlon[idx_xb],xlat[idx_xb],15,stddev_hxb,cmap='gist_heat_r',edgecolors='none',transform=ccrs.PlateCarree(),)
    #cs = ax.scatter(lon,lat,15,corr_2d,cmap='RdBu_r',edgecolors='none',vmin=min_corr,vmax=max_corr,transform=ccrs.PlateCarree(),)
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax.scatter(tc_lon, tc_lat, s=5, marker='*', edgecolors='blue', transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    #color_ticks = np.linspace(min_corr, max_corr, 5, endpoint=True)
    cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, )
    #cbar.set_ticks( color_ticks )
    #cbar.ax.tick_params(labelsize=11)

    #subplot title
    font = {'size':15,}
    ax.set_title( 'Ensemble Spread: IR Tbs', font, fontweight='bold')

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
    if to_obs_res:
        save_des = plot_dir+'stddev_hxb_os_'+DAtime+'_IR_'+sensor+'.png'
    else:
        save_des = plot_dir+'stddev_hxb_ms_'+DAtime+'_IR_'+sensor+'.png'
    
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()
    return None

# ------------------------------------------------------------------------------------------------------
#    Object: ensemble spread of columns of Xb; Operation: Plot the snapshot
# ------------------------------------------------------------------------------------------------------
def plot_2D_stddev_xb_snapshot( lat,lon,stddev_xb_2d ):

    if var_name == 'PSFC':
        stddev_xb_2d = stddev_xb_2d/100
    
    # Replace zeros as NaN
    stddev_xb_2d[stddev_xb_2d == 0] = np.nan

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir,Storm, DAtime)

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5,figsize=(6.5,6), dpi=300)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    #min_corr = -0.5
    #max_corr = 0.5
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax.scatter(lon,lat,15,stddev_xb_2d,cmap='copper_r',edgecolors='none',transform=ccrs.PlateCarree(),)
    #cs = ax.scatter(lon,lat,15,corr_2d,cmap='RdBu_r',edgecolors='none',vmin=min_corr,vmax=max_corr,transform=ccrs.PlateCarree(),)
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax.scatter(tc_lon, tc_lat, s=5, marker='*', edgecolors='blue', transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    #color_ticks = np.linspace(min_corr, max_corr, 5, endpoint=True)
    cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, )
    #cbar.set_ticks( color_ticks )
    #cbar.ax.tick_params(labelsize=11)

    #subplot title
    font = {'size':15,}
    ax.set_title( 'Ensemble Spread: '+var_name, font, fontweight='bold')

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
    if to_obs_res:
        save_des = plot_dir+'stddev_xb_os_'+DAtime+'_'+var_name+'.png'
    else:
        save_des = plot_dir+'stddev_xb_ms_'+DAtime+'_'+var_name+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()

# Plot 3Dcorr per snapshot
def plot_3D_stddev_xb_snapshot( lat,lon,Interp_var,ver_coor,):

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

    #min_corr = -0.6
    #max_corr = 0.6
    for isub in range(20):
        ax.flat[isub].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.flat[isub].coastlines(resolution='10m', color='black',linewidth=0.5)
        cs = ax.flat[isub].scatter(lon,lat,5,Interp_var[isub,:],cmap='RdBu_r',edgecolors='none',transform=ccrs.PlateCarree(),)
        #cs = ax.flat[isub].scatter(lon,lat,5,Interp_corr[isub,:],cmap='RdBu_r',edgecolors='none',vmin=min_corr,vmax=max_corr,transform=ccrs.PlateCarree(),)
        #if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        #    ax.flat[isub].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    #color_ticks = np.linspace(min_corr, max_corr, 5, endpoint=True)
    cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04) #, extend='both')
    #cbar.set_ticks( color_ticks )
    #cbar.ax.tick_params(labelsize=18)

    #subplot title
    font = {'size':15,}
    for isub in range(20):
        ax.flat[isub].set_title( str(ver_coor[isub])+' KM', font, fontweight='bold')

    #title for all
    title_name = Storm+': '+Exper_name+'\nEnsemble Spread: '+var_name
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
    if to_obs_res:
        save_des = plot_dir+'stddev_xb_os_'+DAtime+'_'+var_name+'.png'
    else:
        save_des = plot_dir+'stddev_xb_ms_'+DAtime+'_'+var_name+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


def stddev_xb_snapshot( wrf_dir,DAtime,var_name,var_dim,ver_coor=None):

    if var_dim == '3D':
        nLevel = 42
    elif var_dim == '2D':
        nLevel = 1

    # Read ensemble standard deviation of xb
    des_path = wrf_dir+ "xb_d03_"+var_dim+"_ensStddev_" + DAtime + '_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_xb = pickle.load( f )
    print('Shape of stddev_xb: '+ str(np.shape(stddev_xb)))

    # Find the location of model grid of interest
    if to_obs_res: # nearest for each obs
        idx_xb = Find_nearest_col( Hx_dir, wrf_dir, DAtime, sensor)
    else: # every model grid point
        idx_xb = np.arange(xmax*ymax)

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( wrf_file, 'r')
    # read lon and lat
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()

    if var_dim == '2D':
        if If_plot_stddev_xb_snapshot:
            plot_2D_stddev_xb_snapshot( xlat[idx_xb],xlon[idx_xb],stddev_xb[idx_xb],)
    else:
        if interp_H and not interp_P:
            # interpolate
            Interp_cov = vertical_interp( ncdir,stddev_xb,ver_coor )
            Interp_cov = Interp_cov.reshape( (len(ver_coor),xmax*ymax) )
            print( np.shape( Interp_cov ) )
            if If_plot_stddev_xb_snapshot:
                plot_3D_stddev_xb_snapshot( xlat[idx_xb],xlon[idx_xb],Interp_cov[:,idx_xb],ver_coor)
        elif interp_P and not interp_H:
            pass
            #if If_plot_stddev_xb_snapshot:
            #    plot_3Dcorr_snapshot( xlat[idx_xb],xlon[idx_xb],Interp_cov_xb_hxb,P_of_interest )
        else:
            pass

    return None



if __name__ == '__main__':

    big_dir = '/expanse/lustre/scratch/zuy121/temp_project/Pro2_PSU_MW/' #'/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/expanse/lustre/projects/pen116/zuy121/Pro2_PSU_MW/'  #'/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'CONV'
    MP = 'THO'

    v_interest = [ 'QVAPOR','QSNOW','QCLOUD','QRAIN','QICE','QGRAUP'] #[ 'PSFC',]#'rt_vo']
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709030000'
    end_time_str = '201709030000'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297
   
    # field on model/obs locations 
    to_obs_res = False

    # vertical interpolation if needed
    var_dim = '3D'
    if var_dim == '3D':
        interp_P = False
        P_range = np.arange( 995,49,-20 )
        interp_H = True
        H_range = list(np.arange(1,21,1))

    # limitations
    limit = False
    if limit:  
        num_limit = 3

    If_cal_corr = False
    If_save = True

    If_plot_stddev_hxb_snapshot = False
    If_plot_stddev_xb_snapshot = True
    
    If_plot_corr_snapshot = False
    If_plot_meanCorr = False
    # -------------------------------------------------------    
    
    Exper_name = UD.generate_one_name( Storm,DA,MP )

    if not Consecutive_times:
        DAtimes = ['201709160600','201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Calculate correlations between obs and their nearest model columns
    if If_cal_corr:
        print('------------ Calculate the correlation between IR Tbs and columns of model variables --------------')
        for DAtime in DAtimes:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in v_interest:
                print('Calculate '+var_name+'...')
                cal_2Dcorr_IR_ColVar( DAtime, var_name)

    # Plot the ensemble spread of hxb per snapshot
    if If_plot_stddev_hxb_snapshot:
        print('------------ Plot the ensemble spread of hxb --------------')
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Ens_stddev_hxb/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

        for DAtime in DAtimes:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            print('At '+DAtime)
            plot_stddev_hxb_snapshot( DAtime )

    # Plot the ensemble spread of xb per snapshot
    if If_plot_stddev_xb_snapshot:
        print('------------ Plot the ensemble spread of xb --------------')
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Ens_stddev_xb/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

        for DAtime in DAtimes:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            print('At '+DAtime)
            for var_name in v_interest:
                print('Plot '+var_name+'...')
                if interp_H and not interp_P:
                    stddev_xb_snapshot( wrf_dir,DAtime,var_name,var_dim,H_range)

    # Plot the correlations per snapshot
    if If_plot_corr_snapshot:
        print('------------ Plot the correlation --------------')
        for DAtime in DAtimes:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            print('At '+DAtime)
            for var_name in v_interest:
                print('Plot '+var_name+'...')
                corr_snapshot( DAtime,var_name,var_dim)

    # Plot the domain-averaged correlations
    if If_plot_meanCorr:
        for var_name in v_interest:
            print('Plot '+var_name+'...')
            Domain_ave_corr( var_name )






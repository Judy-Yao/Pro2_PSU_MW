#!/usr/bin/env python3

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

import Util_Vis
import Read_Obspace_IR as ROIR
import corr_2DIR_to3Dmodel as Ens_Stat
import Diagnostics as Diag
import matlab.engine

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

# ------------------------------------------------------------------------------------------------------
#           Object: ensemble correlations of columns of Xb and Hxb in 2D; Operation: Calculation
# ------------------------------------------------------------------------------------------------------
def Find_nearest_col( Hx_dir, wrf_dir, DAtime, sensor):

    start_time=time.process_time()
    print('------- Search for the nearest model grid for each obs ------')
    # Read the obs locations
    ens_file = Hx_dir + "/Hxb_ens_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    d_obs = Ens_Stat.read_obspace( ens_file, sensor ) 
    lon_obs = d_obs['lon_obs']
    lat_obs = d_obs['lat_obs']
    # Read model lon and lat
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    lon_x1d = ncdir.variables['XLONG'][0,0,:]
    lat_x1d = ncdir.variables['XLAT'][0,:,0]
    lon_x = ncdir.variables['XLONG'][0,:,:].flatten()
    lat_x = ncdir.variables['XLAT'][0,:,:].flatten()

    # Loop thru all obs and search each's left- and bottom-nearest model grid along x and y direction 
    # returned is i,j in model domain
    Idx_i = nearest_axis( lon_obs,lon_x1d )
    Idx_j = nearest_axis( lat_obs,lat_x1d )

    # Transform to a 2d meshgrid and find the corresponding idx
    Idx_nearest = []
    for io in range(len(lon_obs)):
        Idx_nearest.append( int(Idx_i[io]*len(lon_x1d)+Idx_j[io]) )
    # check 
    print('Checking the 3000th obs... lon: '+str(lon_obs[99])+' lat: '+str(lat_obs[99]))
    print('Its nearest model grid -- lon: '+str( lon_x[Idx_nearest[99]] )+' lat: '+str( lat_x[Idx_nearest[99]] ))

    if len(np.unique(Idx_nearest)) != len(lon_obs):
        warnings.warn('The nearest model grids might be repetitive!')
        print('Number of obs is '+str(len(lon_obs))+' while the number of the unique model location is '+str(len(np.unique(Idx_nearest))))

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    return Idx_nearest    


def cal_2Dcorr_IR_ColVar( DAtime, var_name):

    if 'Q' in var_name:
        nLevel = 42

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
    if corr_obs_res:
        des_path = Hx_dir+ "Hxb_ensPert_obsRes_" + DAtime + '_' +  sensor + '.pickle'
    else:
        des_path = Hx_dir+ "Hxb_ensPert_modelRes_" + DAtime + '_' +  sensor + '.pickle'
    with open( des_path,'rb' ) as f:
        hxb_ens = pickle.load( f )
    print('Shape of hxb_ens: '+ str(np.shape(hxb_ens)))
    # Read ensemble stand deviation of Hxb
    if corr_obs_res:
        des_path = Hx_dir+ "Hxb_ensStddev_obsRes_" + DAtime + '_' +  sensor + '.pickle'
    else:
        des_path = Hx_dir+ "Hxb_ensStddev_modelRes_" + DAtime + '_' +  sensor + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_hxb = pickle.load( f )
    print('Shape of stddev_hxb: '+ str(np.shape(stddev_hxb)))

    # Find the location of model grid of interest
    if corr_obs_res: # nearest for each obs
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

    # Calculate the correlation between Xb and Hxb
    print('Calculating the correlation between Tbs and ' + var_name + '......' )
    start = time.perf_counter()
    corr_xb_hxb = corr_Tb_toCol( cov_xb_hxb,stddev_xb[:,idx_xb],stddev_hxb )
    end = time.perf_counter()
    print("Elapsed (after compilation) = {}s".format((end - start)))

    # sanity check
    assert  0 <= abs(corr_xb_hxb).all() and abs(corr_xb_hxb).all() <= 1

    # May save the correlation
    if If_save:
        if corr_obs_res:
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

# Plot corr per snapshot
def plot_corr_snapshot(big_dir, small_dir, Storm, Exper_name, var_name, DAtime, sensor, P_of_interest, d_model_res ):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # Read Tbs of Hxb
    Tb_file = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/' + "/mean_obs_res_d03" + DAtime + '_' +  sensor + '.txt'
    d_all = ROIR.read_allTb(Tb_file, sensor )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = ROIR.read_TCvitals(small_dir+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(6.5,6.5), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    ### ---Plot Hxb---
    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )
    ax[0,0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax[0,0].coastlines(resolution='10m', color='black',linewidth=0.5)
    xb_Tb = ax[0,0].scatter(d_all['lon_obs'],d_all['lat_obs'],4,c=d_all['meanYb_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax[0,0].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())
    # Colorbar
    #caxes = fig.add_axes([0.12, 0.1, 0.45, 0.02])
    #Tb_bar = fig.colorbar(xb_Tb,ax=ax[0,0],orientation="horizontal", cax=caxes)
    #Tb_bar.ax.tick_params()

    ### ---Plot correlation---
#    Inno = d_all['Yo_obs']-d_all['meanYb_obs']
    min_corr = -1
    max_corr = 1
    for isub in range(1,4):
        ax.flat[isub].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.flat[isub].coastlines(resolution='10m', color='black',linewidth=0.5)
        cs = ax.flat[isub].scatter(d_model_res['lon_model'],d_model_res['lat_model'],4,c=d_model_res['corr'][isub-1,:],\
            edgecolors='none', cmap='RdBu_r', vmin=min_corr, vmax=max_corr, transform=ccrs.PlateCarree())
            #    edgecolors='none', cmap='RdBu_r',transform=ccrs.PlateCarree())
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
            ax.flat[isub].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

    # Colorbar
    caxes = fig.add_axes([0.2, 0.05, 0.6, 0.02])
    #cb_corr_ticks = np.linspace(min_corr, max_corr, 5, endpoint=True)
    corr_bar = fig.colorbar(cs, ax=ax.flat[2:], orientation="horizontal", cax=caxes)
    corr_bar.ax.tick_params()

    #subplot title
    font = {'size':8,}
    ax.flat[0].set_title('H(Xb)', font, fontweight='bold')
    for isub in range(1,4):
        ax.flat[isub].set_title( str(P_of_interest[isub-1])+' hPa', font, fontweight='bold')

    #title for all
    fig.suptitle(Storm+': '+Exper_name+'(c2p_corr of Qvapor&Tb)', fontsize=8, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(4):
        gl = ax.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top = False
        gl.xlabels_bottom = True
        if j==0 or j==2:
            gl.ylabels_left = True
            gl.ylabels_right = False
        else:
            gl.ylabels_left = False
            gl.ylabels_right = False

        if j==2 or j==3:
            gl.xlabels_bottom = True
            gl.xlabels_top = False
        else:
            gl.xlabels_bottom = False
            gl.xlabels_top = False

        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 8}
        gl.ylabel_style = {'size': 8}


    # Save figures
    des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/IR/corr_'+DAtime+'_col'+var_name+'_'+sensor+'.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)
    plt.close()    
    return None

def corr_snapshot( DAtime,var_name ):

    # Read correlations between a Tb and a column of model var
    if corr_obs_res:
        des_path = wrf_dir+ "d03_2Dcorr_obsRes_Hxb_" + DAtime + '_Column_' +  var_name + '.pickle'
    else:
        des_path = wrf_dir+ "d03_2Dcorr_modelRes_Hxb_" + DAtime + '_Column_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        corr_colxb_hxb = pickle.load( f )
    print('Shape of corr_colxb_hxb: '+ str(np.shape(corr_colxb_hxb)))

    # Find the location of model grid of interest
    if corr_obs_res: # nearest for each obs
        idx_xb = Find_nearest_col( Hx_dir, wrf_dir, DAtime, sensor)
    else: # every model grid point
        idx_xb = np.arange(xmax*ymax)    
    print('Shape of idx_xb: '+ str(np.shape(idx_xb)))

    # Read model attributes
    nLevel = corr_colxb_hxb.shape[0]
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    # pressure levels
    PB = ncdir.variables['PB'][0,:,:,:]
    P = ncdir.variables['P'][0,:,:,:]
    P_hpa_all = (PB + P)/100 # 0 dimension: bottom to top
    P_hpa_all = P_hpa_all.reshape(nLevel,xmax*ymax)
    P_hpa = P_hpa_all[:,idx_xb]
    # lon and lat
    lon_all = ncdir.variables['XLONG'][0,:,:].flatten()
    lon = lon_all[idx_xb]
    lat_all = ncdir.variables['XLAT'][0,:,:].flatten()
    lat = lat_all[idx_xb]

    # Specify pressure levels of interest
    P_of_interest = [100,500,850]

    # Calculate the corr at specified levels
    start_time=time.process_time()
    corr_colxb_hxb_P = np.zeros( (len(P_of_interest),corr_colxb_hxb.shape[1]),  )
    for im in range( corr_colxb_hxb.shape[1] ):
        f_interp = interpolate.interp1d( P_hpa[:,im], corr_colxb_hxb[:,im])
        corr_colxb_hxb_P[:,im] = f_interp( P_of_interest )
    #corr_colxb_hxb_cloud_P = corr_colxb_hxb_cloud[:3,:] # test....
    end_time = time.process_time()
    print ('time needed for the interpolation: ', end_time-start_time, ' seconds')
    print('Min of correlation: '+str(np.amin( corr_colxb_hxb_P )))
    print('Max of correlation: '+str(np.amax( corr_colxb_hxb_P )))
    d_model_res = {'lon_model':lon,'lat_model':lat,'corr':corr_colxb_hxb_P}

    # Plot the correlation
    plot_corr_snapshot(big_dir, small_dir, Storm, Exper_name, var_name, DAtime, sensor, P_of_interest, d_model_res )
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

# Plot the time-averaged domain-mean correlation
def plot_Tmean_VD_corr( ave_corr_overT,geoHkm=None ):

    # Calculate the time-averaged domain-averaged corr
    ave_corr = np.mean( ave_corr_overT, axis=0 )   

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)

    # Manually set discrete values on x and y axis and interpolate data to these values
    ## x axis: correlation value
    x_range = np.arange(-1,1.05,0.05)
    x_axis_rg = range(len(x_range))
    f_interp = interpolate.interp1d( x_range, x_axis_rg)
    loc_inx = f_interp( ave_corr )
    ## y axis: model level height
    if not interp_P:
        y_bottom = np.arange(0,10,1)
        y_middle = np.arange(10,15,0.5)
        y_top = np.arange(15,31,1)
        y_range = np.concatenate( (y_bottom,y_middle,y_top),axis=0 ) # control the scale
        y_axis_rg = range(len(y_range))
        f_interp = interpolate.interp1d( y_range, y_axis_rg)
        loc_iny = f_interp( geoHkm )

        ax.plot( loc_inx,loc_iny )
    else:
        pass 
    
    # set X label
    ax.set_xticks( x_axis_rg[::5] )
    ax.set_xticklabels( ['-1.0','-0.75','-0.50','-0.25','0.0','0.25','0.50','0.75','1.0'],fontsize=15 )    
    # set Y label
    if not interp_P:
        ylabel_like = [0.0,5.0,10.0,11.0,12.0,13.0,14.0,15.0,20.0]
        yticks = []
        list_y_range = list(y_range)
        for it in ylabel_like:
            yticks.append( list_y_range.index(it) ) 
        ax.set_yticks( yticks )
        ax.set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
        ax.set_ylabel('Height (KM)',fontsize=15)
        ax.set_ylim(ymin=0,ymax=25) # cut off data above 25km
    else:
        pass

    # Set title
    if 'Q' in var_name:
        title_name = 'Vertical Corr: '+var_name+' - IR Ch'+ch_list[0]
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

    if 'Q' in var_name:
        nLevel = 42

    # Construct a new array 
    if interp_P and corr_obs_res:
        #idx_xb = Find_nearest_col( Hx_dir, wrf_dir, DAtime, sensor)
        corr_overT = np.zeros( [len(DAtimes),len(P_of_interest),xmax*ymax] )
    elif not interp_P and corr_obs_res:
        pass
    elif interp_P and not corr_obs_res: 
        pass
    else: # at model level
        corr_overT = np.zeros( [len(DAtimes),nLevel] )
        
    # Load correlation over time and interpolate if necessary
    for DAtime in DAtimes:

        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
        t_idx = DAtimes.index( DAtime )

        # Read correlations between a Tb and a column of model var
        if corr_obs_res:
            des_path = wrf_dir+ "d03_2Dcorr_obsRes_Hxb_" + DAtime + '_Column_' +  var_name + '.pickle'
        else:
            des_path = wrf_dir+ "d03_2Dcorr_modelRes_Hxb_" + DAtime + '_Column_' +  var_name + '.pickle' 
        with open( des_path,'rb' ) as f:
            corr_colxb_hxb = pickle.load( f )
        print('Shape of corr_colxb_hxb: '+ str(np.shape(corr_colxb_hxb)))

        mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
        ncdir = nc.Dataset( mean_xb, 'r')
        # ---------- at model level (geo height) ------------------
        if not interp_P: 
            PHB = ncdir.variables['PHB'][0,:,:,:]
            PH = ncdir.variables['PH'][0,:,:,:]
            geoHkm = (PHB+PH)/9.8/1000 # in km
            geoHkm = geoHkm.reshape( geoHkm.shape[0],-1)
            geoHkm_Dmean = np.mean( geoHkm, axis=1 )
            geoHkm_half_eta = (geoHkm_Dmean[:-1]+geoHkm_Dmean[1:])/2
            geoHkm_half_eta = np.ma.getdata(geoHkm_half_eta)
            corr_overT[t_idx,:] = np.mean( corr_colxb_hxb, axis=1 )                         
        # ---------- at interpolated levels ----------
        else:
            PB = ncdir.variables['PB'][0,:,:,:]
            P = ncdir.variables['P'][0,:,:,:]
            P_hpa = ((PB + P)/100)
            P_hpa = P_hpa.reshape( P_hpa.shape[0],-1) # 0 dimension: bottom to top
            # Interpolate to P level of interest
            start_time=time.process_time()
            for im in range( corr_colxb_hxb.shape[1] ):
                f_interp = interpolate.interp1d( P_hpa[:,im], corr_colxb_hxb[:,im])
                corr_overT[t_idx,:,im] = f_interp( P_of_interest )
            end_time = time.process_time()
            print ('time needed for the interpolation: ', end_time-start_time, ' seconds') 
    
    # Average the domain
    if not interp_P:
        ave_corr_overT = corr_overT
    else:
        ave_corr_overT = np.mean( corr_overT, axis=2 )  

    # Plot the time-averaged and time-evolution of vertical distribution of corr 
    if not interp_P:
        #plot_Tmean_VD_corr( ave_corr_overT,geoHkm_half_eta )
        plot_var_incre_timeseries( ave_corr_overT,geoHkm_half_eta )
    else:
        #plot_Tmean_VD_corr( ave_corr_overT )
        plot_var_incre_timeseries( ave_corr_overT )

    return None


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'MARIA'
    Exper_name = 'IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900'

    v_interest = [ 'QVAPOR',]
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709160000'
    end_time_str = '201709180000'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297
    
    corr_obs_res = False
    ens_Interp_to_obs = False

    If_cal_pert_stddev = False
    If_cal_corr = False
    interp_P = False
    P_of_interest = list(range( 995,49,-20 ))
    If_save = True

    If_plot_snapshot = False
    If_plot_meanCorr = True
    # -------------------------------------------------------    

    if not Consecutive_times:
        DAtimes = ['201709160600','201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]


    # Interpolate simulated Tb at model resolution to obs resolution
    if ens_Interp_to_obs:
        print('------- For all members, interpolate Hx in model resolution to obs location ------')
        for DAtime in DAtimes:
            # Read assimilated obs 
            file_Diag = big_dir+Storm+'/'+Exper_name+'/run/'+DAtime+'/enkf/d03/fort.10000'
            d_obs = Diag.Find_IR( file_Diag, fort_v )
            # Interpolate
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            Ens_Stat.interp_simu_to_obs_matlab_ens( d_obs, Hx_dir, sensor, ch_list,  DAtime )

    # Calculate ensemble perturbations and variances
    if If_cal_pert_stddev:
        print('------------ Calculate the ensemble perturbations --------------')
        for DAtime in DAtimes:
            # Hxb
            #Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            #if corr_obs_res:
            #    print('At obs space...')
            #    Ens_Stat.cal_pert_stddev_obsRes_Hxb( DAtime, sensor, Hx_dir, If_save)
            #else:
            #    print('At model space...')
            #    Ens_Stat.cal_pert_stddev_modelRes_Hxb( DAtime, sensor, Hx_dir, If_save)
            # Xb
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in v_interest:
                Ens_Stat.cal_pert_stddev_xb( DAtime, wrf_dir, var_name, If_save )


    # Calculate correlations between obs and their nearest model columns
    if If_cal_corr:
        print('------------ Calculate the correlation between IR Tbs and columns of model variables --------------')
        for DAtime in DAtimes:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in v_interest:
                print('Calculate '+var_name+'...')
                cal_2Dcorr_IR_ColVar( DAtime, var_name)


    # Plot the correlations per snapshot
    if If_plot_snapshot:
        print('------------ Plot the correlation --------------')
        for DAtime in DAtimes:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            print('At '+DAtime)
            for var_name in v_interest:
                print('Plot '+var_name+'...')
                corr_snapshot( DAtime,var_name )

    # Plot the domain-averaged correlations
    if If_plot_meanCorr:
        for var_name in v_interest:
            print('Plot '+var_name+'...')
            Domain_ave_corr( var_name )






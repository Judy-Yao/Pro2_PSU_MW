#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
from wrf import getvar, interplevel
# It might be possible that you are not able to conda install wrf-var with a pretty new python version
# Solution:
# 1. conda create -n $PYTHON34_ENV_NAME python=3.4 anaconda 
# 2. conda activate python=3.4 (use wrf-python in this python environment)
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
import subprocess

import Read_Obspace_IR as ROIR
from Util_Vis import HydroIncre
import Diagnostics as Diag


def plot_var_incre_timeseries( ave_var_overT, geoHkm=None ):

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)
    # Set up coordinates
    if not interp_P:
        xv = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]
        #y_bottom = np.arange(0,10,1)
        #y_middle = np.arange(10,15,0.5)
        #y_top = np.arange(15,31,1)
        #y_range = np.concatenate( (y_bottom,y_middle,y_top),axis=0 ) # control the scale
        y_bottom = np.arange(0,5,0.5)
        y_top = np.arange(5,31,1)
        y_range = np.concatenate( (y_bottom,y_top),axis=0 )

        y_axis_rg = range(len(y_range))
        f_interp = interpolate.interp1d( y_range, y_axis_rg)
        yv = f_interp( geoHkm )
        xcoor, ycoor = np.meshgrid( xv, yv )
    else:
        xv = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]#range( np.shape(ave_norm_overT)[0] )
        yv = range( np.shape(ave_var_overT)[1])
        xcoor, ycoor = np.meshgrid( xv, yv )

    
    if_all_posi = all( it[0] > 0 for it in ave_var_overT.tolist() )
    if_all_nega = all( it[0] < 0 for it in ave_var_overT.tolist() )
    if if_all_posi or if_all_nega is True:
        cmap = 'jet'
    else:
        cmap = 'bwr'
    
    bounds = np.linspace(-0.2,0.2,5)
    #bounds = np.linspace( np.floor(ave_var_overT.min()),np.ceil(ave_var_overT.max()),10 )
    #print(bounds)
    bounds_format = [ "{0:.2f}".format( item ) for item in bounds]
    #incre_contourf = ax.contourf( xcoor, ycoor, np.transpose( ave_var_overT ), cmap=cmap, vmin=ave_var_overT.min(), vmax=ave_var_overT.max(),levels=bounds_format,extend='both')
    incre_contourf = ax.contourf( xcoor, ycoor, np.transpose( ave_var_overT ), cmap=cmap,levels=bounds_format,extend='both')
    # Add color bar below the plot
    color_bar = fig.colorbar(incre_contourf,orientation = 'horizontal',pad=0.15,ticks=bounds)
    bounds_str =  [ str(item) for item in bounds_format ]
    color_bar.ax.set_xticklabels( bounds_str, rotation=45)
    color_bar.ax.tick_params(labelsize=12)
    color_bar.ax.set_xlabel('Domain-mean Increment',fontsize=15)

    # Set X/Y labels
    # set X label
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45)
    # set Y label
    if not interp_P:
        ylabel_like = [0.0,1.0,2.0,3.0,4.0,5.0,10.0,15.0,20.0]
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
        ax.set_yticklabels( [str(it) for it in P_of_interest[::10]],fontsize=15 )
        ax.set_ylabel('Pressure (hPa)',fontsize=15)

    # Set title
    if 'rt_vo' in var_name:
        title_name = 'EnKF Increment: relative vorticity ($10^-5$ s-1)'
    else:
        pass
    ax.set_title( title_name,fontweight="bold",fontsize='12' )
    fig.suptitle(Storm+': '+Exper_name, fontsize=10, fontweight='bold')

    # Save the figure
    if not interp_P:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/ML_increment_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    else:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/Interp_increment_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


    return None

def eachVar_plot( ):

    # Dimension
    xmax = 297
    ymax = 297
    nLevel = 42

    if interp_P: # ---------- Interpolate to specified pressure levels ----------
        # Construct a new array (using interpolation)
        ave_var_overT = np.zeros( [len(DAtimes),len(P_of_interest)] )
        #ave_T_profile = np.zeros( [len(DAtimes),len(P_of_interest)] )
    else:
        # Construct a new array at model level
        ave_var_overT = np.zeros( [len(DAtimes),nLevel] )
        #ave_T_profile = np.zeros( [len(DAtimes),nLevel] )

    for DAtime in DAtimes:
        print('At ', DAtime)

        # Get the time index
        t_idx = np.where([DAtime == it for it in DAtimes])[0]

        # Read increment of variable of interest
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
        if var_name == 'rt_vo': # relative vorticity
            # analysis
            xa_file = wrf_dir+'/wrf_enkf_output_d03_mean'
            xa_ncdir = nc.Dataset(xa_file, 'r')
            xa_avo = getvar( xa_ncdir, 'avo') # Absolute vorticity, units: 10-5 s-1 
            xa_coriolis_sin = xa_ncdir.variables['F'][0,:,:]/1e-5 #units: 10-5 s-1
            xa_rtvo = xa_avo - xa_coriolis_sin #relative vorticity, units: 10-5 s-1 
            # background
            xb_file = wrf_dir+'/wrf_enkf_input_d03_mean'
            xb_ncdir = nc.Dataset(xb_file, 'r')
            xb_avo = getvar( xb_ncdir, 'avo') # Absolute vorticity, units: 10-5 s-1 
            xb_coriolis_sin = xb_ncdir.variables['F'][0,:,:]/1e-5 #units: 10-5 s-1
            xb_rtvo = xb_avo - xb_coriolis_sin #relative vorticity, units: 10-5 s-1 
            # increment
            var = xa_rtvo - xb_rtvo
            var = np.array( var )
            var = var.reshape( var.shape[0],-1) # n_level, n_points
        else:
            pass

        # Set up coordinate info
        if interp_P: # Interpolate to P level of interest
            pass
        else:
            # Read height
            mean_xa = wrf_dir+'/wrf_enkf_output_d03_mean'
            ncdir = nc.Dataset( mean_xa, 'r')
            PHB = ncdir.variables['PHB'][0,:,:,:]
            PH = ncdir.variables['PH'][0,:,:,:]
            geoHkm = (PHB+PH)/9.8/1000 # in km
            geoHkm = geoHkm.reshape( geoHkm.shape[0],-1)
            geoHkm_Dmean = np.mean( geoHkm, axis=1 )
            geoHkm_half_eta = (geoHkm_Dmean[:-1]+geoHkm_Dmean[1:])/2
            geoHkm_half_eta = np.ma.getdata(geoHkm_half_eta)
            # Perform domain mean
            var_mean = np.mean( var,axis=1 )
            #idx_zero = np.where( abs(var_mean) <= 1e-8 )[0]
            #var_mean[idx_zero] = 0
            ave_var_overT[t_idx,:] = var_mean

    # Plot the increament in a time series
    if interp_P:
        plot_var_incre_timeseries( ave_var_overT )
    else:
        plot_var_incre_timeseries( ave_var_overT, geoHkm_half_eta )

    return None



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    Exper_name = 'IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'
    v_interest = ['rt_vo',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709030000'
    end_time_str = '201709040000'
    Consecutive_times = True

    interp_P = False
    P_of_interest = list(range( 850,50,-20 ))
    interp_H = True
    interp_to_obs = False

    If_ncdiff = False
    If_plot = True
    # -------------------------------------------------------    

    # Identify DA times in the period of interest
    if not Consecutive_times:
        DAtimes = ['201709140000',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Plot the time evolution of domain-averaged increments of dynamic fields
    if If_plot:
        start_time=time.process_time()
        for var_name in v_interest:
            print('Plot '+var_name+'...')
            eachVar_plot( )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')






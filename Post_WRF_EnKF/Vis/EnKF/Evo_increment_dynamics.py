#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
from wrf import getvar, ll_to_xy
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
import scipy as sp
import scipy.ndimage
import time
import subprocess
import pickle

import metpy.calc as mpcalc
from metpy.units import units

import Read_Obspace_IR as ROIR
import Util_data as UD
import Diagnostics as Diag

def find_minSLP( wrfout ):
    ncdir = nc.Dataset( wrfout )
    slp = getvar(ncdir, 'slp')
    min_slp = np.amin( slp )
    slp_smooth = sp.ndimage.gaussian_filter(slp, [11,11])
    idx = np.nanargmin( slp_smooth )
    lat_minslp = ncdir.variables['XLAT'][:].flatten()[idx]
    lon_minslp = ncdir.variables['XLONG'][:].flatten()[idx]

    min_slp = [lat_minslp,lon_minslp,min_slp.values]
    return min_slp

def plot_2Dvar_incre_timeseries( ave_var_overT):

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)

    # Plot the figure
    xv = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]
    ax.plot( xv,ave_var_overT,marker='o',linewidth=4,color='black' )
    
    ax.axhspan(0, max(ave_var_overT), color='blue', alpha=0.2)
    ax.axhspan(min(ave_var_overT), 0, color='red', alpha=0.2)
    # Set X/Y labels
    ax.set_ylim([min(ave_var_overT),max(ave_var_overT)])
    # set X label
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45, labelsize=10)

    # Set title
    if 'min_slp' in var_name:
        title_name = 'EnKF Increment: value of min SLP'
    else:
        pass
    ax.set_title( title_name,fontweight="bold",fontsize='11' )
    fig.suptitle(Storm+': '+Exper_name, fontsize=10, fontweight='bold')

    # Save the figure
    save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


    return None


def plot_3Dvar_incre_timeseries( ave_var_overT, geoHkm=None ):

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
    
    bounds = np.linspace(-0.2,0.2,11) #np.linspace(-0.3,0.3,7)
    #bounds = np.linspace( np.floor(ave_var_overT.min()),np.ceil(ave_var_overT.max()),10 )
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
    ax.tick_params(axis='x', labelrotation=45, labelsize=10)
    # set Y label
    if not interp_P:
        ylabel_like = [0.0,1.0,2.0,3.0,4.0,5.0,10.0,15.0,20.0]
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
    if specify_area:
        if 'rt_vo' in var_name:
            title_name = 'EnKF Increment: relative vorticity ($10^-5$ s-1) in a circle with R='+str(radius_threshold)+' KM'
        elif 'diver' in var_name:
            title_name = 'EnKF Increment: divergence ($10^-5$ s-1) in a circle with R='+str(radius_threshold)+' KM'
    else:
        pass
    ax.set_title( title_name,fontweight="bold",fontsize='11' )
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

# Area-mean calculation of 2D variable
def each2DVar_timeSeries( ):

    # Dimension
    xmax = 297
    ymax = 297
    nLevel = 1

    # Initiate a new array 
    ave_var_overT = np.zeros( [len(DAtimes),] )

    if specify_area:
        # Read the best track position every 6 hours
        d_btk = UD.read_bestrack(Storm)
        # Read the hourly best track position 
        dict_btk = UD.interpolate_locations( DAtimes, d_btk)

    for DAtime in DAtimes:
        print('At ', DAtime)

        # Get the time index
        t_idx = np.where([DAtime == it for it in DAtimes])[0]

        # Read increment of variable of interest
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime

        # might only use data in a specified area
        if specify_area:
            # Find the best-track position 
            btk_dt = [it_str for it_str in dict_btk['time'] ]#[datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]
            bool_match = [DAtime == it for it in btk_dt]
            if True in bool_match:
                idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
            else:
                idx_btk  = None
            # convert from ll to xy
            wrf_file = wrf_dir+'/wrf_enkf_input_d03_mean'
            ncdir = nc.Dataset( wrf_file )
            tc_ij = ll_to_xy(ncdir, dict_btk['lat'][idx_btk], dict_btk['lon'][idx_btk])
            # What ll_to_xy returns is not the xy coordinate itself but the grid index starting from 0. 
            # (https://forum.mmm.ucar.edu/threads/what-does-wrf-python-function-ll_to_xy-returns.12248/)
            tc_i = tc_ij.values[0]
            tc_j = tc_ij.values[1]
            idx_x = UD.find_circle_area_model_ij( wrf_file, tc_i, tc_j, radius_threshold, 3)
        else:
            idx_x = np.arange(xmax*ymax)

        # variable increment
        if var_name == 'min_slp': # sea level pressure
            # analysis
            xa_file = wrf_dir+'/wrf_enkf_output_d03_mean'
            xa_minSLP = find_minSLP( xa_file )
            # background
            xb_file = wrf_dir+'/wrf_enkf_input_d03_mean'
            xb_minSLP = find_minSLP( xb_file )
            # increment
            var = xa_minSLP[2] - xb_minSLP[2]
        else:
            pass

        # area-mean calculation
        if var_name == 'min_slp':
            ave_var_overT[t_idx] = var
        else:
            var_mean = np.nanmean( var )
            ave_var_overT[t_idx] = var_mean
    
    # May save the data
    if If_save:
        # Metadata
        current_datetime = datetime.now()
        formatted_datetime = current_datetime.strftime('%Y-%m-%d %H:%M:%S')
        metadata = {'created_at':formatted_datetime}
        save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
        # create a dictionary with metadata and data
        meta_and_data = {'metadata':metadata,'ave_var_overT':ave_var_overT}
        # Write the dictionary to a pickle file
        with open(save_des,'wb') as file:
            pickle.dump( meta_and_data, file )
        print( 'Saving the data: ', save_des )


    # Plot the time series of increment     
    if If_plot_series:
        plot_2Dvar_incre_timeseries( ave_var_overT )

    return None


# Area-mean calculation of 3D variable
def each3DVar_timeSeries_cal( ):

    # Dimension
    xmax = 297
    ymax = 297
    nLevel = 42

    if interp_P: # ---------- Interpolate to specified pressure levels ----------
        # Construct a new array (using interpolation)
        ave_var_overT = np.zeros( [len(DAtimes),len(P_of_interest)] )
    else:
        # Construct a new array at model level
        ave_var_overT = np.zeros( [len(DAtimes),nLevel] )

    if specify_area:
        # Read the best track position every 6 hours
        d_btk = UD.read_bestrack(Storm)
        # Read the hourly best track position 
        dict_btk = UD.interpolate_locations( DAtimes, d_btk)

    for DAtime in DAtimes:
        print('At ', DAtime)

        # Get the time index
        t_idx = np.where([DAtime == it for it in DAtimes])[0]

        # Read increment of variable of interest
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime

        # might only use data in a specified area
        if specify_area:
            # Find the best-track position 
            btk_dt = [it_str for it_str in dict_btk['time'] ]#[datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]
            bool_match = [DAtime == it for it in btk_dt]
            if True in bool_match:
                idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
            else:
                idx_btk  = None
            # convert from ll to xy
            wrf_file = wrf_dir+'/wrf_enkf_input_d03_mean'
            ncdir = nc.Dataset( wrf_file )
            tc_ij = ll_to_xy(ncdir, dict_btk['lat'][idx_btk], dict_btk['lon'][idx_btk])
            # What ll_to_xy returns is not the xy coordinate itself but the grid index starting from 0. 
            # (https://forum.mmm.ucar.edu/threads/what-does-wrf-python-function-ll_to_xy-returns.12248/)
            tc_i = tc_ij.values[0]
            tc_j = tc_ij.values[1]
            idx_x = UD.find_circle_area_model_ij( wrf_file, tc_i, tc_j, radius_threshold, 3)     
        else:
            idx_x = np.arange(xmax*ymax)

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
            var = var[:,idx_x]
        elif var_name == 'diver':
            # analysis
            xa_file = wrf_dir+'/wrf_enkf_output_d03_mean'
            xa_ncdir = nc.Dataset(xa_file, 'r')
            Um_xa = getvar( xa_ncdir, 'ua') # U-component of wind on mass points
            Vm_xa = getvar( xa_ncdir, 'va') # V-component of wind on mass points
            Um_xa = Um_xa* units("m/s")
            Vm_xa = Vm_xa* units("m/s")
            dvg_xa = mpcalc.divergence(Um_xa, Vm_xa, dx=3000* units("m"), dy=3000* units("m"))
            # background
            xb_file = wrf_dir+'/wrf_enkf_input_d03_mean'
            xb_ncdir = nc.Dataset(xb_file, 'r')
            Um_xb = getvar( xb_ncdir, 'ua') # U-component of wind on mass points
            Vm_xb = getvar( xb_ncdir, 'va') # V-component of wind on mass points
            Um_xb = Um_xb* units("m/s")
            Vm_xb = Vm_xb* units("m/s")
            dvg_xb = mpcalc.divergence(Um_xb, Vm_xb, dx=3000* units("m"), dy=3000* units("m"))
            # increment
            var = (dvg_xa - dvg_xb)*1e5 # unit: 10-5 s-1
            var = np.array( var )
            var = var.reshape( var.shape[0],-1) # n_level, n_points
            var = var[:,idx_x]
        else:
            pass

        # Set up coordinate info
        if interp_P: # Interpolate to P level of interest
            # Read pressure levels
            mean_xa = wrf_dir+'/wrf_enkf_output_d03_mean'
            ncdir = nc.Dataset( mean_xa, 'r')
            PB = ncdir.variables['PB'][0,:,:,:]
            P = ncdir.variables['P'][0,:,:,:]
            P_hpa = (PB + P)/100
            P_hpa = P_hpa.reshape( P_hpa.shape[0],-1)
            P_hpa = P_hpa[:,idx_x] # only select value in the specified area
            # Quality control: the P_hpa of lowest level is less than 900 mb 
            idx_bad = np.where( P_hpa[0,:] < 900 )[0]
            idx_all = range( P_hpa.shape[1] )
            idx_good = np.delete(idx_all, idx_bad)
            good_P_hpa = P_hpa[:,idx_good]
            good_var = var[:,idx_good]
            # Interpolate 
            T_interp = np.zeros( [len(P_of_interest),len(idx_good)] )
            var_interp = np.zeros( [len(P_of_interest),len(idx_good)] )
            start_time=time.process_time()
            for im in range( len(idx_good) ):
                f_interp = interpolate.interp1d( good_P_hpa[:,im], good_var[:,im] )
                var_interp[:,im] = f_interp( P_of_interest )
            end_time = time.process_time()
            print ('time needed for the interpolation: ', end_time-start_time, ' seconds')
            # Process for mixing ratio
            var_mean = np.nanmean( var_interp,axis=1 )
            # Perform domain mean
            ave_var_overT[t_idx,:] = var_mean
        else:
            # Read height
            mean_xa = wrf_dir+'/wrf_enkf_output_d03_mean'
            ncdir = nc.Dataset( mean_xa, 'r')
            PHB = ncdir.variables['PHB'][0,:,:,:]
            PH = ncdir.variables['PH'][0,:,:,:]
            geoHkm = (PHB+PH)/9.8/1000 # in km
            geoHkm = geoHkm.reshape( geoHkm.shape[0],-1)
            geoHkm = geoHkm[:,idx_x]
            geoHkm_Dmean = np.nanmean( geoHkm, axis=1 )
            geoHkm_half_eta = (geoHkm_Dmean[:-1]+geoHkm_Dmean[1:])/2
            geoHkm_half_eta = np.ma.getdata(geoHkm_half_eta)
            # Perform domain mean
            var_mean = np.nanmean( var,axis=1 )
            ave_var_overT[t_idx,:] = var_mean

    # May save the data
    if If_save:
        # Metadata
        current_datetime = datetime.now()
        formatted_datetime = current_datetime.strftime('%Y-%m-%d %H:%M:%S')
        # create data
        if interp_P:
            metadata = {'created_at':formatted_datetime, 'Interpolated_to': 'Pressure (hPa)','Interpolated_at':P_of_interest,'radius_threshold':radius_threshold}
            save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/Interp_increment_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
            # create a dictionary with metadata and data
            meta_and_data = {'metadata':metadata,'ave_var_overT':ave_var_overT}
        else:
            metadata = {'created_at':formatted_datetime,'radius_threshold':radius_threshold}
            save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/ML_increment_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
            meta_and_data = {'metadata':metadata,'ave_var_overT':ave_var_overT,'geoHkm_half_eta':geoHkm_half_eta}
        # Write the dictionary to a pickle file
        with open(save_des,'wb') as file:
            pickle.dump( meta_and_data, file )
        print( 'Saving the data: ', save_des )

    # Plot the increament in a time series
    if If_plot_series:
        if interp_P:
            plot_3Dvar_incre_timeseries( ave_var_overT )
        else:
            plot_3Dvar_incre_timeseries( ave_var_overT, geoHkm_half_eta )

    return None



if __name__ == '__main__':

    big_dir = '/scratch_S2/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    MP = 'TuneWSM6'
    DA = 'IR'
    v_interest = ['rt_vo',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709030000'
    end_time_str = '201709050000'
    Consecutive_times = True

    interp_P = True
    P_of_interest = list(range( 900,10,-20 ))
    interp_H = False
    specify_area = True
    radius_threshold = 300 #km

    If_ncdiff = False
    If_cal_series = True
    If_save = True
    If_plot_series = True
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

    # Plot the time evolution of domain-averaged increments
    if If_cal_series:
        start_time=time.process_time()
        for var_name in v_interest:
            print('Calculate '+var_name+'...')
            if var_name == 'min_slp':
                each2DVar_timeSeries( )
            else:
                each3DVar_timeSeries_cal( )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')

    # Plot the time evolution of domain-averaged increments
    if If_plot_series:
        for var_name in v_interest:
            print('Plot '+var_name+'...')

            if var_name == 'min_slp':
                save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
                with open(save_des,'rb') as file:
                    meta_and_data = pickle.load( file )
                ave_var_overT = meta_and_data['ave_var_overT']
                plot_2Dvar_incre_timeseries( ave_var_overT )
            else:
                # for 3D var only
                if interp_P:
                    save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/Interp_increment_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
                    with open(save_des,'rb') as file:
                        meta_and_data = pickle.load( file )
                    ave_var_overT = meta_and_data['ave_var_overT']
                    #ave_T_profile = meta_and_data['ave_T_profile']
                    plot_3Dvar_incre_timeseries( ave_var_overT )
                else:
                    save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/ML_increment_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
                    with open(save_des,'rb') as file:
                        meta_and_data = pickle.load( file )
                    ave_var_overT = meta_and_data['ave_var_overT']
                    #ave_T_profile = meta_and_data['ave_T_profile']
                    geoHkm_half_eta = meta_and_data['geoHkm_half_eta']
                    plot_3Dvar_incre_timeseries( ave_var_overT,geoHkm_half_eta )









#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

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

import Read_Obspace_IR as ROIR
from Util_Vis import HydroIncre
import Diagnostics as Diag


# Generate time series
def generate_times( Storms, start_time_str, end_time_str ):

    dict_times = {}
    for istorm in Storms:
        time_diff = datetime.strptime(end_time_str[istorm],"%Y%m%d%H%M") - datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        dict_times[istorm] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times


# ------------------------------------------------------------------------------------------------------
#           Operation: Plot the domain-averaged value of each variable in a time series
# ------------------------------------------------------------------------------------------------------


# Plot increment itself
def plot_var_incre_timeseries( ave_var_overT,N_times,geoHkm=None ):

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)
    # Set up coordinates
    if not interp_P:
        xv = range( np.shape(ave_var_overT)[0] ) 
        y_bottom = np.arange(0,5,0.5)
        y_top = np.arange(5,31,1)
        y_range = np.concatenate( (y_bottom,y_top),axis=0 )

        y_axis_rg = range(len(y_range))
        f_interp = interpolate.interp1d( y_range, y_axis_rg)
        yv = f_interp( geoHkm )
        xcoor, ycoor = np.meshgrid( xv, yv )
    else:
        xv = range( np.shape(ave_var_overT)[0] )
        yv = range( np.shape(ave_var_overT)[1])
        xcoor, ycoor = np.meshgrid( xv, yv )

    if_all_posi = all( it[0] > 0 for it in ave_var_overT.tolist() )
    if_all_nega = all( it[0] < 0 for it in ave_var_overT.tolist() )
    if if_all_posi or if_all_nega is True:
        cmap = 'jet'
    else:
        cmap = 'bwr'
    bounds = np.linspace(-0.2,0.2,5)
    bounds_format = [ "{0:.2f}".format( item ) for item in bounds]
    incre_contourf = ax.contourf( xcoor, ycoor, np.transpose( ave_var_overT ), cmap=cmap,levels=bounds_format,extend='both')
    # Add color bar below the plot
    color_bar = fig.colorbar(incre_contourf,orientation = 'horizontal',pad=0.15,ticks=bounds)
    bounds_str =  [ str(item) for item in bounds_format ]
    color_bar.ax.set_xticklabels( bounds_str, rotation=45)
    color_bar.ax.tick_params(labelsize=12)
    color_bar.ax.set_xlabel('Domain-mean Increment',fontsize=15)

    # Plot T profile
    if not interp_P:
        pass
    else:
        pass
        #T_contour = ax.contour( xcoor[0:-5,:], ycoor[0:-5,:], np.transpose(ave_T_profile[:,0:-5]-273.15),colors='k')
        #ax.clabel(T_contour, inline=True)

    # Set X/Y labels
    # set X label
    range_xtick = np.linspace(0,N_times,6)
    #range_xtick_label = range_xtick+1
    ax.set_xticks( range_xtick )
    ax.set_xticklabels( [str(int(it)) for it in range_xtick],fontsize=15 )
    ax.set_xlabel('WRF-EnKF Cycles',fontsize=15)
    ax.set_xlim(xmin=0,xmax=N_times-1) 
    #ax.tick_params(axis='x', labelrotation=45)
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
    if 'rt_vo' in var_name:
        title_name = 'EnKF Increment: relative vorticity ($10^-5$ s-1) in the circle area (center@best-track) with R='+str(radius_threshold)+' KM'
    else:
        pass
    ax.set_title( title_name,fontweight="bold",fontsize='11' )
    fig.suptitle('Storms: '+MP, fontsize=15, fontweight='bold')

    # Save the figure
    if not interp_P:
        save_des = small_dir+'/SYSTEMS/Vis_analyze/'+MP+'_ML_increment_'+var_name+'_'+str(N_times)+'cycles.png'
    else:
        save_des = small_dir+'/SYSTEMS/Vis_analyze/'+MP+'_Interp_increment_'+var_name+'_'+str(N_times)+'cycles.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()
    return None

# ------------------------------------------------------------------------------------------------------
#           Operation: Read the domain-averaged increment of individual storm and average over storms
# ------------------------------------------------------------------------------------------------------
def eachVar_timeSeries( Exper_names,dict_times,var_name ): 

    if 'Q' in var_name:
        nLevel = 42

    # Check the length of DAtimes is the same for all storms
    for i in range(len(Storms)-1):
        if len(dict_times[Storms[i]]) != len(dict_times[Storms[i+1]]):
            raise ValueError('The length of times is not equal between experiments')
    N_times = len(dict_times[Storms[i]])

    # Array for the mean over all storms
    if interp_P: 
        ave_var = np.zeros( [N_times,len(P_of_interest)] )
    else:
        ave_var = np.zeros( [N_times,nLevel] )
        geoHkm = np.zeros( [N_times,nLevel] )        

    # Loop thru storms
    for ist in Storms:
        print('Loading data for '+ ist)
        if interp_P:
            save_des = small_dir+ist+'/'+Exper_names[ist]+'/Data_analyze/EnKF/Interp_increment_'+var_name+'_'+start_time_str[ist]+'_'+end_time_str[ist]+'.pickle'
            print(save_des)
            # Read data from a pickle file
            with open(save_des,'rb') as file:
                load_data = pickle.load( file )
            ave_var = ave_var + load_data['ave_var_overT'] 
            # Make average over storms
            ave_var = ave_var/len(Storms)
            #ave_T = ave_T/len(Storms)
        else:
            save_des = small_dir+ist+'/'+Exper_names[ist]+'/Data_analyze/EnKF/ML_increment_'+var_name+'_'+start_time_str[ist]+'_'+end_time_str[ist]+'.pickle'
            # Read data from a pickle file
            with open(save_des,'r') as file:
                load_data = pickle.load( file )
            ave_var = ave_var + load_data['ave_var_overT'] 
            geoHkm = geoHkm + load_data['geoHkm_half_eta'] 
            # Make average over storms
            ave_var = ave_var/len(Storms)
            #ave_T = ave_T/len(Storms)
            geoHkm = geoHkm/len(Storms)

    # Plot the increament in a time series
    if If_plot_series:
        if interp_P:
            plot_var_incre_timeseries( ave_var, N_times )
        else:
            plot_var_incre_timeseries( ave_var, N_times,geoHkm )

    return None





if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    MP = 'THO'
    v_interest = ['rt_vo']
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708231200','IRMA':'201709040000','JOSE':'201709060000','MARIA':'201709170000'}
    Consecutive_times = True

    interp_P = True
    P_of_interest = list(range( 900,10,-20 ))
    interp_H = False
    radius_threshold = 250 

    If_plot_series = True
    # -------------------------------------------------------    

    # Create experiment names
    Exper_names = generate_names( MP,Storms )        

    # Identify DA times in the period of interest
    dict_times = generate_times( Storms, start_time_str, end_time_str )

    # Plot time series of domain-average increment for each variable
    if If_plot_series:
        start_time=time.process_time()
        for var_name in v_interest:
            print('Calculate '+var_name+'...')
            eachVar_timeSeries( Exper_names,dict_times,var_name )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')












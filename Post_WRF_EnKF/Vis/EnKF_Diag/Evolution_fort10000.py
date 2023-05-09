#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os,sys,stat # functions for interacting with the operating system
import subprocess
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
import Diagnostics as Diag
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from matplotlib import colors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time


# setting font sizeto 30
plt.rcParams.update({'font.size': 15})


def plot_slp_timeseries( small_dir, Storm, Exper, DAtimes, Evo_slp ):

    # Set up figure
    fig = plt.figure( figsize=(12,9), dpi=300 )
    ax = plt.subplot(1,1,1)
    dates = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]
    # convert the unit from Pa to hPa
    obs_hpa = Evo_slp['obs']/100 
    yb_hpa = Evo_slp['prior_mean']/100 
    ya_hpa = Evo_slp['posterior_mean']/100 
    yb_spead_hpa = Evo_slp['prior_spread']/100 
    ya_spead_hpa = Evo_slp['posterior_spread']/100 
    # Plot obs and mean
    ax.plot_date(dates, obs_hpa, 'black', label='TCvital', linestyle='-')
    ax.plot_date(dates, yb_hpa, 'blue', label='Prior', linestyle='-')
    ax.plot_date(dates, ya_hpa, 'red', label='Posterior', linestyle='-')
    leg = plt.legend(loc='upper right',fontsize=15)
    # Plot the spread area
    ax.fill_between(dates, yb_hpa-yb_spead_hpa, yb_hpa+yb_spead_hpa, color='blue', alpha=0.2)
    ax.fill_between(dates, ya_hpa-ya_spead_hpa, ya_hpa+ya_spead_hpa, color='red', alpha=0.2)
    # Set X/Y labels
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45,labelsize=15)
    ax.set_ylabel('Minimum Sea Level Pressure (hPa)',fontsize=15)
    ax.set_ylim( 1005,1018 )

    title_name = Storm+'('+Exper_name+')'+': mim SLP'
    ax.set_title( title_name,fontweight="bold",fontsize='15' )

    # Save the figure
    save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/minslp_fort10000_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()

# ---------------------------------------------------------------------------------------------------------
#           Operation: Assemble records from files at different times into one 
# ---------------------------------------------------------------------------------------------------------
def Gather_slp( big_dir, Storm, Exper, DAtimes, v_interest ):

    slp_evo = [[] for i in range(len(v_interest))] 
    
    for DAtime in DAtimes:
        file_diag = big_dir+Storm+'/'+Exper+'/run/'+DAtime+'/enkf/d03/fort.10000'
        it = DAtimes.index( DAtime )
        slp_tmp = Diag.Find_min_slp( file_diag, v_interest )        
        for var in v_interest:
            idx_var = v_interest.index(var)
            slp_evo[idx_var].append( slp_tmp[var][0] )

    Evo_slp = {}
    print('Sanity check of reading process...Printing 1st record...')
    for key in v_interest:
        idx_var = v_interest.index( key )
        Evo_slp[key] = np.array( slp_evo[idx_var] )
        print(key, ':', Evo_slp[key][0])

    return Evo_slp


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/' #'/scratch/02191/yuz31/
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # Configuration
    Storm = 'MARIA'
    Exper_name = 'IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr/'
    v_interest = ['obs_type','lat','lon','obs','prior_mean','posterior_mean','prior_spread','posterior_spread']
    start_time_str = '201709160000'
    end_time_str = '201709160900'
    Consecutive_times = True

    Check_slp = True

    # Identify DA times in the period of interest
    if not Consecutive_times:
        DAtimes = ['201709140000',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Investigate the evolution
    if Check_slp:
        # put together slp at all times into a dictionary
        Evo_slp = Gather_slp( big_dir, Storm, Exper_name, DAtimes, v_interest )
        plot_slp_timeseries( small_dir, Storm, Exper_name, DAtimes, Evo_slp )

import os # functions for interacting with the operating system
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
from wrf import getvar, interplevel
import math
import scipy as sp
import scipy.ndimage
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import subprocess
import metpy.calc as mpcalc
from metpy.units import units

import Util_data as UD


# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot the evolution of min slp
# ------------------------------------------------------------------------------------------------------
def find_minSLP( wrfout ):
    ncdir = nc.Dataset( wrfout )
    slp = getvar(ncdir, 'slp')
    min_slp = np.amin( slp )
    slp_smooth = sp.ndimage.gaussian_filter(slp, [11,11])
    idx = np.nanargmin( slp_smooth )
    lat_minslp = ncdir.variables['XLAT'][:].flatten()[idx]
    lon_minslp = ncdir.variables['XLONG'][:].flatten()[idx]

    minSLP = [lat_minslp,lon_minslp,min_slp.values]
    return minSLP

def plot_slp_timeseries( small_dir, Storm, Expers, DAtimes, Evo_slp ):

    # Set up figure
    fig = plt.figure( figsize=(12,9), dpi=300 )
    ax = plt.subplot(1,1,1)
    dates = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]
    # Plot obs and mean
    print(Evo_slp['obs_slp'])
    ax.plot_date(dates, Evo_slp['obs_slp'], 'black', label='TCvital', linestyle='-')
    ax.plot_date(dates, Evo_slp['xa_slp'][0,:], 'red', label='conv_THO', linestyle='--')

    leg = plt.legend(loc='lower left',fontsize=15)
    # Set X/Y labels
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45,labelsize=15)
    ax.set_ylabel('Minimum Sea Level Pressure (hPa)',fontsize=15)
    ax.set_ylim( 980,1020)

    #title_name = Storm+'('+Exper_name+')'+': mim SLP'
    title_name = Storm+': mim SLP'
    ax.set_title( title_name,fontweight="bold",fontsize='15' )
    #fig.suptitle('conv+HPI', fontsize=15, fontweight='bold')


    # Save the figure
    save_des = small_dir+Storm+'/'+Expers[0]+'/Vis_analyze/Model/minslp_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()

def Gather_slp( Storm, Expers, DAtimes, big_dir ):

    obs_minslp_lat = []
    obs_minslp_lon = []
    obs_minslp_value = []

    xa_minslp_lat = [[] for i in range(len(Expers))]
    xa_minslp_lon = [[] for i in range(len(Expers))]
    xa_minslp_value = [[] for i in range(len(Expers))]

    for DAtime in DAtimes:

        # collect min slp found from WRF output
        for Exper in Expers:
            wrf_dir = big_dir+Storm+'/'+Exper+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
            print('Reading the EnKF posterior mean from ', wrf_dir)
            list_wrfout = find_minSLP( wrf_dir )
            idx = Expers.index( Exper )
            xa_minslp_lat[idx].append( list_wrfout[0] )
            xa_minslp_lon[idx].append( list_wrfout[1] )
            xa_minslp_value[idx].append( list_wrfout[2] )

        # collect assimilated min slp obs from TCvital record in fort.10000
        diag_enkf = big_dir+Storm+'/'+Expers[0]+'/run/'+DAtime+'/enkf/d03/fort.10000'
        print('Reading the EnKF diagnostics from ', diag_enkf)
        enkf_minSlp = subprocess.run(["grep","slp",diag_enkf],stdout=subprocess.PIPE,text=True)
        list_enkf_minSlp = enkf_minSlp.stdout.split()
        # condition on the number of assimilated min slp
        if list_enkf_minSlp.count('slp') == 0 :
            raise ValueError('No min slp is assimilated!')
        elif list_enkf_minSlp.count('slp') == 1 :
            obs_minslp_lat.append( float(list_enkf_minSlp[1]) )
            obs_minslp_lon.append( float(list_enkf_minSlp[2]) )
            obs_minslp_value.append( float(list_enkf_minSlp[9])/100 )
        else : # at least two min slp records are assimilated
            print('At least two min slp obs are assimilated!')
            # find the index of the current time
            idx_time = DAtimes.index( DAtime )
            # assemble the diagnosed min slp from an analysis
            xa_ms_lat = xa_minslp_lat[0][idx_time]
            xa_ms_lon = xa_minslp_lon[0][idx_time]
            # ---condition 1: find the nearest TCvital min slp from the analysis
            # find the index/location of 'slp' in fort.10000
            indices = [i for i ,e in enumerate(list_enkf_minSlp) if e == 'slp']
            # assemble a pair of coordinate for each 'slp'
            distances = []
            obs_slp = []
            for it in indices:
                obs_slp.append( float(list_enkf_minSlp[it+9]) )
                lon1 = float(list_enkf_minSlp[it+2])
                lat1 = float(list_enkf_minSlp[it+1])
                distances.append( UD.mercator_distance(lon1, lat1, xa_ms_lon, xa_ms_lat) )
            min_distances = [np.amin(distances) == it for it in distances]
            idx_nearest = min_distances.index(True)
            # ---condition 2: the min slp
            min_obs_slp = [np.amin(obs_slp) == it for it in obs_slp]
            idx_min = min_obs_slp.index(True)
            # ---combine the two conditions
            if idx_min == idx_nearest:
                idx_coor = idx_min
                print('Storm center is choosed with two condtions met!')
            else:
                print('Not sure which obs is the storm center. Go with the min value one!')
                idx_coor = idx_min

            # gather this TCvital min slp
            obs_minslp_lat.append( float(list_enkf_minSlp[indices[idx_coor]+1]) )
            obs_minslp_lon.append( float(list_enkf_minSlp[indices[idx_coor]+2]) )
            obs_minslp_value.append( float(list_enkf_minSlp[indices[idx_coor]+9])/100 )

    dict_minSLP = {'obs_slp':np.array(obs_minslp_value),'obs_lat':np.array(obs_minslp_lat),'obs_lon':np.array(obs_minslp_lon),'xa_slp':np.array(xa_minslp_value),'xa_lat':np.array(xa_minslp_lat), 'xa_lon':np.array(xa_minslp_lon)}
    return dict_minSLP


# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot the evolution of IC water
# ------------------------------------------------------------------------------------------------------




if __name__ == '__main__':

    big_dir = '/scratch_S2/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # -------- Configuration -----------------
    Storm = 'IRMA'
    DA = 'IR'
    MP = 'WSM6'

    Plot_minslp_evo = False
    # -----------------------------------------

    # Time range set up
    start_time_str = '201709030000'
    end_time_str = '201709040000'
    Consecutive_times = True

    if not Consecutive_times:
        DAtimes = ['201709041400','201709041600']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    ## Experiment name
    Expers = UD.generate_one_name( Storm,DA,MP ) 

    # Plot the evolution of minimum sea level pressure
    if Plot_minslp_evo:
        Evo_slp = Gather_slp( Storm, Expers, DAtimes, big_dir )
        plot_slp_timeseries( small_dir, Storm, Expers, DAtimes, Evo_slp )




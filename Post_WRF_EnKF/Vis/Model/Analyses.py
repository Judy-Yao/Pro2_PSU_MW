#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
from wrf import getvar
# It might be possible that you are not able to conda install wrf-var with a pretty new python version
# Solution:
# 1. conda create -n $PYTHON34_ENV_NAME python=3.4 anaconda 
# 2. conda activate python=3.4 (use wrf-python in this python environment)
import math
import matlab.engine
import scipy as sp
import scipy.ndimage
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from Track_intensity import read_bestrack
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time

def d03_domain(wrfout_d03):
    ncdir = nc.Dataset(wrfout_d03, 'r')

    xlat = ncdir.variables['XLAT'][0,:,:]
    xlong = ncdir.variables['XLONG'][0,:,:]

    d03_lat_min = np.min( xlat.flatten() )
    d03_lat_max = np.max( xlat.flatten() )
    d03_lon_min = np.min( xlong.flatten() )
    d03_lon_max = np.max( xlong.flatten() )

    d03_list = [d03_lon_min, d03_lon_max, d03_lat_min, d03_lat_max]
    return d03_list

# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot UV10_slp data
# ------------------------------------------------------------------------------------------------------
def plot_UV10_slp( Storm, DAtime, wrfout, plot_dir ):

    # Read storm center
    dict_btk = read_bestrack(Storm)
    # Find the best-track position
    btk_dt = [it_str for it_str in dict_btk['time'] ]#[datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]
    bool_match = [DAtime == it for it in btk_dt]
    if True in bool_match:
        if_btk_exist = True
        idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
    else:
        if_btk_exist = False

    # ------ Read WRFout -------------------
    ncdir = nc.Dataset( wrfout )

    # domain
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    lat_min = np.amin(lat)
    lon_min = np.amin(lon)
    lat_max = np.amax(lat)
    lon_max = np.amax(lon)
    # sea level pressure
    slp = getvar(ncdir, 'slp')
    min_slp = np.amin( slp )
    max_slp = np.amax( slp )
    slp_smooth = sp.ndimage.gaussian_filter(slp, [11,11])
    idx = np.nanargmin( slp_smooth )
    lat_minslp = ncdir.variables['XLAT'][:].flatten()[idx]
    lon_minslp = ncdir.variables['XLONG'][:].flatten()[idx]
    # Wind at 10 meters
    u10 = ncdir.variables['U10'][0,:,:]
    v10 = ncdir.variables['V10'][0,:,:]
    windspeed = (u10 ** 2 + v10 ** 2) ** 0.5

    # figure
    fig = plt.figure()

    ax = plt.subplot(1,1,1,projection=ccrs.PlateCarree())
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax.coastlines (resolution='10m', color='black', linewidth=1)
    # sea level pressure
    slp_contour = ax.contour(lon,lat,slp_smooth,cmap='Greys_r',vmin=min_slp,vmax=max_slp,transform=ccrs.PlateCarree())
    plt.clabel(slp_contour, inline=1, fontsize=9)
    # Wind at 10 meters
    wind_smooth = sp.ndimage.gaussian_filter(windspeed, [2,2])
    min_wind = 0
    max_wind = 25
    bounds = np.linspace(min_wind, max_wind, 6)
    wind_contourf = ax.contourf(lon,lat,wind_smooth,cmap='hot_r',vmin=min_wind,vmax=max_wind,levels=bounds,extend='both',transform=ccrs.PlateCarree())
    # Adding the colorbar
    cbaxes = fig.add_axes([0.05, 0.1, 0.03, 0.8])
    wind_bar = fig.colorbar(wind_contourf,cax=cbaxes,fraction=0.046, pad=0.04) #Make a colorbar for the ContourSet returned by the contourf call.
    wind_bar.ax.set_ylabel('Wind Speed (m/s)')
    wind_bar.ax.tick_params(labelsize=7)
    ax.barbs(lon.flatten(), lat.flatten(), u10.flatten(), v10.flatten(), length=5, pivot='middle',
         color='royalblue', regrid_shape=20, transform=ccrs.PlateCarree())
    # Mark the best track
    if if_btk_exist:
        ax.scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 8, 'green', marker='*',transform=ccrs.PlateCarree())

    # Title
    wrfout_head_tail = os.path.split( wrfout )
    ax.set_title(wrfout_head_tail[1]+': UV10_slp',  fontweight='bold', fontsize=10)
    #ax.set_title(wrfout_head_tail[1].replace('wrfout_d03_',' '),  fontweight='bold', fontsize=10)

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2))
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.xlabels_bottom = True
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 6}
    gl.ylabel_style = {'size': 6}

    plt.savefig( plot_dir+wrfout_head_tail[1]+'_'+DAtime+'_UV10_slp.png', dpi=300 )
    print('Saving the figure: ', plot_dir+wrfout_head_tail[1]+'_'+DAtime+'_UV10_slp.png')
    plt.close()

def UV10_slp( Storm, Exper_name, DAtimes, big_dir, small_dir ):

    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
        print('Reading WRF analysis: ', wrf_dir)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/UV10_slp/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            plot_UV10_slp( Storm, DAtime, wrf_dir, plot_dir )
        else:
            plot_UV10_slp( Storm, DAtime, wrf_dir, plot_dir )

if __name__ == '__main__':

    Storm = 'MARIA'
    Exper_name = ['IR-J_DA+J_WRF+J_init-SP-intel17',]
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    Plot_UV10_slp = True

    # Time range set up
    start_time_str = '2017091400'
    end_time_str = '2017091700'
    Consecutive_times = True

    if not Consecutive_times:
        DAtimes = ['201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Plot low-level circulation
    if Plot_UV10_slp:
        for iExper in Exper_name:
            UV10_slp( Storm, iExper, DAtimes, big_dir, small_dir )


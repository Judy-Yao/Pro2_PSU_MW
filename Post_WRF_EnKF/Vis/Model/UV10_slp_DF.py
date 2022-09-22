#!/usr/bin/env python3

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import pickle
from netCDF4 import Dataset
from wrf import getvar
# It might be possible that you are not able to conda install wrf-var with a pretty new python version
# Solution:
# 1. conda create -n $PYTHON34_ENV_NAME python=3.4 anaconda 
# 2. conda activate python=3.4 (use wrf-python in this python environment)
import math
import scipy as sp
import scipy.ndimage
import matplotlib
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


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

def plot_UV10_slp( Hx_file ):
    ncdir = nc.Dataset( Hx_file )
    
    # domain
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    min_lat = np.amin(lat)
    min_lon = np.amin(lon)
    max_lat = np.amax(lat)
    max_lon = np.amax(lon)
    # sea level pressure
    slp = getvar(ncdir, 'slp')
    slp_smooth = sp.ndimage.gaussian_filter(slp, [11,11])
    min_slp = np.amin(slp_smooth)
    max_slp = np.amax( slp_smooth )
    # Wind at 10 meters
    u10 = ncdir.variables['U10'][0,:,:]
    v10 = ncdir.variables['V10'][0,:,:]
    windspeed = (u10 ** 2 + v10 ** 2) ** 0.5

    # figure
    fig = plt.figure()

    ax = plt.subplot(1,1,1,projection=ccrs.PlateCarree())
    ax.set_extent([min_lon,max_lon,min_lat,max_lat], crs=ccrs.PlateCarree())
    ax.coastlines (resolution='10m', color='black', linewidth=1)
    # sea level pressure
    contours = ax.contour(lon,lat,slp_smooth,cmap='Greys',vmin=min_slp,vmax=max_slp,transform=ccrs.PlateCarree())
    plt.clabel(contours, inline=1, fontsize=10)
    # Wind at 10 meters
    wind_smooth = sp.ndimage.gaussian_filter(windspeed, [2,2])
    min_wind = np.amin(windspeed)
    max_wind = np.amax(windspeed)
    wind_contourf = ax.contourf(lon,lat,wind_smooth,cmap='jet',vmin=min_wind,vmax=max_wind,transform=ccrs.PlateCarree())
    cbar = fig.colorbar(wind_contourf) #Make a colorbar for the ContourSet returned by the contourf call.
    cbar.ax.set_ylabel('Wind Speed (m/s)')
    cbar.ax.tick_params(labelsize=6)
    ax.barbs(lon.flatten(), lat.flatten(), u10.flatten(), v10.flatten(), pivot='middle',
         color='black', regrid_shape=10, transform=ccrs.PlateCarree())
    # Plot rectangle of domain 3
    d03_name = wrfout.replace('d02','d03')
    [d03_lon_min, d03_lon_max, d03_lat_min, d03_lat_max] = d03_domain( d03_name )
    rec_x = [d03_lon_max, d03_lon_max, d03_lon_min, d03_lon_min, d03_lon_max]
    rec_y = [d03_lat_max, d03_lat_min, d03_lat_min, d03_lat_max, d03_lat_max]
    ax.plot( rec_x,rec_y,color='white',linewidth=1,marker='.',transform=ccrs.PlateCarree())

    ax.set_title(wrfout_d01.replace('wrfout_d02_',' '),  fontweight='bold')

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
    gl.xlabel_style = {'size': 4}
    gl.ylabel_style = {'size': 6}

    plt.savefig( wrfout_d01+'.png', dpi=300 )
    plt.close()

def plot_d02_d03():
    # List all of wrfd02 files
    wrfout_list = glob.glob( Hx_file+'/wrfout_d01_2017*' )
    for wrfout in wrfout_list:
        plot()






if __name__ == '__main__':
    Storm = 'MARIA'
    Exper_name = 'newWRF_MW_THO'
    wrf_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'

    DFtimes = []
    plot_UV10_slp( Storm, Exper_name, wrf_dir )


#!/usr/bin/env python3

import os
import glob
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math


def location_domain(wrfout_d03):
    ncdir = nc.Dataset(wrfout_d03, 'r')

    xlat = ncdir.variables['XLAT_M'][0,:,:]
    xlong = ncdir.variables['XLONG_M'][0,:,:]

    d03_lat_min = np.min( xlat.flatten() )
    d03_lat_max = np.max( xlat.flatten() )
    d03_lon_min = np.min( xlong.flatten() )
    d03_lon_max = np.max( xlong.flatten() ) 

    d03_list = [d03_lon_min, d03_lon_max, d03_lat_min, d03_lat_max]
    return d03_list



wrfout_d01_list = glob.glob('met_em.d01.*')

for wrfout_d01 in wrfout_d01_list:
    
    ncdir = nc.Dataset(wrfout_d01, 'r')
    d01_skintemp = ncdir.variables['SKINTEMP'][:,:]
    d01_lat = ncdir.variables['XLAT_M'][0,:,:]
    d01_lon = ncdir.variables['XLONG_M'][0,:,:]
 
    d01_domain_geo = location_domain(wrfout_d01)
    d01_lon_min = d01_domain_geo[0]
    d01_lon_max = d01_domain_geo[1]
    d01_lat_min = d01_domain_geo[2]
    d01_lat_max = d01_domain_geo[3]

    d02_name = wrfout_d01.replace('d01','d02')
    d02_domain_geo = location_domain(d02_name)
    d02_lon_min = d02_domain_geo[0]
    d02_lon_max = d02_domain_geo[1]
    d02_lat_min = d02_domain_geo[2]
    d02_lat_max = d02_domain_geo[3]

    d03_name = wrfout_d01.replace('d01','d03')
    d03_domain_geo = location_domain(d03_name)
    d03_lon_min = d03_domain_geo[0]
    d03_lon_max = d03_domain_geo[1]
    d03_lat_min = d03_domain_geo[2]
    d03_lat_max = d03_domain_geo[3]

    fig = plt.figure()
    ax = plt.subplot(1,1,1,projection=ccrs.PlateCarree())

    ax.set_extent([d01_lon_min,d01_lon_max,d01_lat_min,d01_lat_max], crs=ccrs.PlateCarree())
    ax.coastlines( resolution='10m', color='black',linewidth=1 )
    cs = ax.scatter(d01_lon.flatten(),d01_lat.flatten(),2,d01_skintemp.flatten(),cmap=plt.cm.get_cmap('jet'),vmin=280, vmax=330,transform=ccrs.PlateCarree())
    cbar = fig.colorbar(cs, orientation="horizontal")
    cbar.ax.tick_params(labelsize=6) 
    # Plot rectangle
    rec_x = [d02_lon_max, d02_lon_max, d02_lon_min, d02_lon_min, d02_lon_max]
    rec_y = [d02_lat_max, d02_lat_min, d02_lat_min, d02_lat_max, d02_lat_max]
    ax.plot( rec_x,rec_y,color='purple',linewidth=1.5,marker='.',transform=ccrs.PlateCarree())
    # Plot rectangle
    rec_x = [d03_lon_max, d03_lon_max, d03_lon_min, d03_lon_min, d03_lon_max]
    rec_y = [d03_lat_max, d03_lat_min, d03_lat_min, d03_lat_max, d03_lat_max]
    ax.plot( rec_x,rec_y,color='white',linewidth=1.5,marker='.',transform=ccrs.PlateCarree())    

    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.xlabels_bottom = True
    gl.ylabels_left = True
    gl.ylabels_right = False


    # Axis labels
    lon_ticks = list(range(math.ceil(d01_lon_min)-2, math.ceil(d01_lon_max)+2,2))
    lat_ticks = list(range(math.ceil(d01_lat_min)-2, math.ceil(d01_lat_max)+2,2))

    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 4}
    gl.ylabel_style = {'size': 6}

    ax.set_title(wrfout_d01.replace('wrfout_d01_',' '),  fontweight='bold')

    plt.savefig( wrfout_d01+'.png', dpi=300 )
    plt.close()

#!/usr/bin/env python3

import os
import glob
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
#from metpy.plots import ctables
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math


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


# Get the 1st domain as the reference
#ncdir = nc.Dataset('wrfout_d02_2017-08-22_18:00:00', 'r')
#d01_lat = ncdir.variables['XLAT'][0,:,:]
#d01_lon = ncdir.variables['XLONG'][0,:,:]
#d01_lat_min = np.min( d01_lat.flatten() )
#d01_lat_max = np.max( d01_lat.flatten() )
#d01_lon_min = np.min( d01_lon.flatten() )
#d01_lon_max = np.max( d01_lon.flatten() )

wrfout_d01_list = glob.glob('wrfout_d01_2017-08*')
#wrfout_d01_list = ['wrfout_d02_2017-08-23_12:00:00','wrfout_d02_2017-08-23_13:00:00']

for wrfout_d01 in wrfout_d01_list:
    ncdir = nc.Dataset(wrfout_d01, 'r')
    d01_reflect = ncdir.variables['REFL_10CM'][0,0,:,:]
    d01_lat = ncdir.variables['XLAT'][0,:,:]
    d01_lon = ncdir.variables['XLONG'][0,:,:]
    d01_lat_min = np.min( d01_lat.flatten() )
    d01_lat_max = np.max( d01_lat.flatten() )
    d01_lon_min = np.min( d01_lon.flatten() )
    d01_lon_max = np.max( d01_lon.flatten() )    


    d03_name = wrfout_d01.replace('d01','d03')
    d03_domain_geo = d03_domain(d03_name)
    d03_lon_min = d03_domain_geo[0]
    d03_lon_max = d03_domain_geo[1]
    d03_lat_min = d03_domain_geo[2]
    d03_lat_max = d03_domain_geo[3]

    fig = plt.figure()
    ax = plt.subplot(1,1,1,projection=ccrs.PlateCarree())

    ax.set_extent([d01_lon_min,d01_lon_max,d01_lat_min,d01_lat_max], crs=ccrs.PlateCarree())
    ax.coastlines( resolution='10m', color='black',linewidth=1 )
    #cmap = ctables.get_colortable('NEWReflectivity')
    cs = ax.scatter(d01_lon.flatten(),d01_lat.flatten(),2,d01_reflect.flatten(),cmap=plt.cm.get_cmap('jet'),vmin=-25, vmax=50,transform=ccrs.PlateCarree())
    cbar = fig.colorbar(cs, orientation="horizontal")
    cbar.ax.tick_params(labelsize=6) 
    # Plot rectangle
    rec_x = [d03_lon_max, d03_lon_max, d03_lon_min, d03_lon_min, d03_lon_max]
    rec_y = [d03_lat_max, d03_lat_min, d03_lat_min, d03_lat_max, d03_lat_max]
    ax.plot( rec_x,rec_y,color='white',linewidth=1.5,marker='.',transform=ccrs.PlateCarree())
    
    ax.set_title(wrfout_d01.replace('wrfout_d01_',' '),  fontweight='bold')

    plt.savefig( wrfout_d01+'.png', dpi=300 )
    plt.close()

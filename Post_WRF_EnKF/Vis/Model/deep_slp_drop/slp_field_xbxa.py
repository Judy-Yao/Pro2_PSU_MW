#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
from wrf import getvar, interplevel
# It might be possible that you are not able to conda install wrf-var with a pretty new python version
# Solution:
# 1. conda create -n $PYTHON34_ENV_NAME python=3.4 anaconda 
# 2. conda activate python=3.4 (use wrf-python in this python environment)
import math
import scipy as sp
import scipy.ndimage
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
import matplotlib.patches as patches
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time

import Util_data as UD


def read_slp( wrf_dir ):

    # set file names for background and analysis
    mean_dir = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
    # read the shared variables: lat/lon
    ncdir = nc.Dataset( mean_dir[0] )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # read UV10 and slp
    slp_original = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    slp_smooth = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    for i in range(len(mean_dir)):
        file_name = mean_dir[i]
        ncdir = nc.Dataset( file_name )
        # sea level pressure
        slp_original[i,:,:] = UD.compute_slp( ncdir )#getvar(ncdir, 'slp')
        slp_smooth[i,:,:] = sp.ndimage.gaussian_filter(slp_original[i,:,:], [11,11])

    d_slp = {'lat':lat,'lon':lon,'slp':slp_original,'slp_smooth':slp_smooth}
    return d_slp


def plot_slp( Storm, Exper_name, DAtime, wrf_dir, plot_dir ):

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(Storm, DAtime)
    # ------ Read WRFout -------------------
    d_field = read_slp( wrf_dir )
    lat = d_field['lat']
    lon = d_field['lon']
    slp = d_field['slp']
    slp_smooth = d_field['slp_smooth']
    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )
    min_slp = np.min( slp )
    max_slp = np.max( slp )

    # ------ Plot Figure -------------------
    fig, ax=plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(12,6), dpi=400)

    # find a box with min slp of enkf output less than 900 hpa
    #idx = np.where( slp[1,:,:].flatten() <= 960)[0]
    #min_lonx = np.amin( lon.flatten()[idx] )-0.05
    #max_lonx = np.amax( lon.flatten()[idx] )+0.05
    #min_latx = np.amin( lat.flatten()[idx] )-0.05
    #max_latx = np.amax( lat.flatten()[idx] )+0.05
    #path = [ [min_lonx,min_latx],[min_lonx,max_latx],[max_lonx,max_latx],[max_lonx,min_latx], ]

    for i in range(2):
        ax[i].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i].coastlines (resolution='10m', color='black', linewidth=1)
        # sea level pressure
        min_slp = 900
        max_slp = 1020
        bounds = np.linspace(min_slp, max_slp, 7)
        #start_level = max( np.amin(slp_smooth[0,:,:]), np.amin(slp_smooth[1,:,:]) )+0.1
        #end_level = min( np.amax(slp_smooth[0,:,:]), np.amax(slp_smooth[1,:,:]) )
        #step = (end_level -start_level)/5
        #level = np.arange(start_level, end_level,step)
        #slp_contour = ax[i].contour(lon,lat,slp_smooth[i,:,:],levels=level,cmap='Greys_r',vmin=min_slp,vmax=max_slp,transform=ccrs.PlateCarree())
        slp_contourf = ax[i].contourf(lon,lat,slp[i,:,:],cmap='magma',vmin=min_slp,vmax=max_slp,levels=bounds,transform=ccrs.PlateCarree())
        # Mark the best track
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
            ax[i].scatter(tc_lon,tc_lat, 20, 'black', marker='*',transform=ccrs.PlateCarree())
            ax[i].text(tc_lon-1,tc_lat-0.5,tc_slp,fontsize=10)
        #ax[i].add_patch(patches.Polygon(path,facecolor='none',edgecolor='black',linewidth=1.5 ))

    # Adding the colorbar
    cbaxes = fig.add_axes([0.01, 0.1, 0.02, 0.8])
    wind_bar = fig.colorbar(slp_contourf,cax=cbaxes,fraction=0.046, pad=0.04) #Make a colorbar for the ContourSet returned by the contourf call.
    wind_bar.set_clim( vmin=min_slp, vmax=max_slp )
    wind_bar.ax.set_ylabel('slp (hPa)')
    wind_bar.ax.tick_params(labelsize='12')

    # Title
    ax[0].set_title( 'Xb--min slp: '+str("{0:.3f}".format(np.min( slp[0,:,:] )))+' hPa',  fontweight='bold') #, fontsize=12)
    ax[1].set_title( 'Xa--min slp: '+str("{0:.3f}".format(np.min( slp[1,:,:] )))+' hPa',  fontweight='bold') #, fontsize=12)
    fig.suptitle(Storm+': '+Exper_name+'('+DAtime+')', fontsize=12, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2))
    for j in range(2):
        gl = ax[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.xlabels_bottom = True
        if j==0:
            gl.ylabels_left = True
            gl.ylabels_right = False
        else:
            gl.ylabels_left = False
            gl.ylabels_right = False
        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 12}
        gl.ylabel_style = {'size': 12}

    des = plot_dir+DAtime+'_slp_original.png'
    plt.savefig( des, dpi=300 )
    print('Saving the figure: ', des)
    plt.close()


if __name__ == '__main__':

    Storm = 'HARVEY'
    Exper_name = 'JerryRun/IR_WSM6'
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    Plot_slp = True

    # Time range set up
    start_time_str = '201708241200'
    end_time_str = '201708251200'
    Consecutive_times = True

    if not Consecutive_times:
        DAtimes = ['201708251200',]
        #DAtimes = ['201708230000','201708230600','201708231200','201708231800','201708240000','201708240600','201708241200']
        #DAtimes = ['201709031200','201709031800','201709040000','201709040600','201709041200','201709041800','201709050000']
        #DAtimes = ['201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Plot low-level circulation
    if Plot_slp:
        start_time=time.process_time()
        # Loop through each DAtime/analysis
        for DAtime in DAtimes:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
            print('Reading WRF background and analysis: ', wrf_dir)
            DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
            # ------ Plot -------------------
            plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/study_slp/'
            plotdir_exists = os.path.exists( plot_dir )
            if plotdir_exists == False:
                os.mkdir(plot_dir)
                plot_slp( Storm, Exper_name, DAtime, wrf_dir, plot_dir )
            else:
                plot_slp( Storm, Exper_name, DAtime, wrf_dir, plot_dir )        
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')

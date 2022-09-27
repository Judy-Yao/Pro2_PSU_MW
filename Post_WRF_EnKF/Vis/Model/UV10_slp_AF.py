#!/usr/bin/env python3

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import pickle
import netCDF4 as nc
from wrf import getvar
# It might be possible that you are not able to conda install wrf-var with a pretty new python version
# Solution:
# 1. conda create -n $PYTHON34_ENV_NAME python=3.4 anaconda 
# 2. conda activate python=3.4 (use wrf-python in this python environment)
import USAF
import math
import scipy as sp
import scipy.ndimage
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import copy

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

def plot_UV10_slp( wrfout,plot_dir, dict_AF_masked):

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
    ax.scatter(lon_minslp, lat_minslp, s=5, marker='*', edgecolors='red', transform=ccrs.PlateCarree())
    # Plot rectangle of domain 3
    #d03_name = wrfout.replace('d02','d03')
    #[d03_lon_min, d03_lon_max, d03_lat_min, d03_lat_max] = d03_domain( d03_name )
    #rec_x = [d03_lon_max, d03_lon_max, d03_lon_min, d03_lon_min, d03_lon_max]
    #rec_y = [d03_lat_max, d03_lat_min, d03_lat_min, d03_lat_max, d03_lat_max]
    #ax.plot( rec_x,rec_y,color='darkgreen',linewidth=1,marker='.',transform=ccrs.PlateCarree())
    # Plot aircraft observation
    if dict_AF_masked is not None:
        min_height = np.amin( dict_AF_masked['gpsa'][:] )
        max_height = np.amax( dict_AF_masked['gpsa'][:] )
        AFAF = ax.scatter(dict_AF_masked['lon'][:],dict_AF_masked['lat'][:],3,dict_AF_masked['gpsa'][:],cmap='Greens_r',transform=ccrs.PlateCarree()) 
        cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.8])
        AFAF_bar = fig.colorbar(AFAF,cax=cbaxes,fraction=0.046, pad=0.04)
        AFAF_bar.ax.set_ylabel('Flight Height (m)')
        AFAF_bar.ax.tick_params(labelsize=7)
    # Title
    wrfout_head_tail = os.path.split( wrfout )
    ax.set_title(wrfout_head_tail[1].replace('wrfout_d03_',' '),  fontweight='bold')

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

    plt.savefig( plot_dir+wrfout_head_tail[1]+'.png', dpi=300 )
    print('Saving the figure: ', plot_dir+wrfout_head_tail[1]+'.png')
    plt.close()

#def plot_d02_d03(wrf_dir, plot_dir, USAF_per_day):
#    wrfout_list = sorted(glob.glob( wrf_dir+'/wrfout_d02_2017*00' ))
def plot_d03(wrf_dir, plot_dir, USAF_per_day):
    wrfout_list = sorted(glob.glob( wrf_dir+'/wrfout_d03_2017*00' ))
    wrfout_all_time = []
    for wrfout in wrfout_list:
        wrfout_head_tail_all =  os.path.split( wrfout ) 
        wrfout_all_time.append( wrfout_head_tail_all[1][11:21].replace('-','') ) # eg.,20170919
    
    # Separate wrffiles of days and read the USAFAF file for each day
    for wrf_time in USAF_per_day.keys(): 
        print('On the day: ', wrf_time)
        bool_idx = [it == wrf_time for it in wrfout_all_time ]
        idx_wrf_at_day = np.where(bool_idx)[0]

        attr_interest_usaf = ['GMT','GPSA','LAT','LON']
        if USAF_per_day[wrf_time] is not None:
            dict_AF_all_today, uni_hh = USAF.read_USAF(wrf_time, USAF_per_day[wrf_time], attr_interest_usaf)
        else:
            dict_AF_all_today = None
    
        wrfout_filter = [wrfout_list[it] for it in idx_wrf_at_day]
        for wrfout in wrfout_filter:
            print( '--------- Plotting ', wrfout, '------------' )
            # check if UASA obs exisits for this hour
            wrfout_head_tail_all =  os.path.split( wrfout )
            wrf_hh = int(wrfout_head_tail_all[1][22:24])
            if dict_AF_all_today is not None:
                print('USAF obs exists on this day!')
                if wrf_hh in uni_hh:
                    print('USAF obs exists at this hour!')
                    # mask time
                    Nminutes = 15
                    dict_AF_masked = USAF.mask_time( dict_AF_all_today, wrf_hh, Nminutes)
                else:
                    print('USAF AF obs does not exist at this hour!')
                    dict_AF_masked = None
            else:
                print('USAF AF obs does not exist on this day!')
                dict_AF_masked = None
            # plot
            plot_UV10_slp( wrfout, plot_dir, dict_AF_masked )

class Info_USAF_clt:
    '''This object is designed to collect information for each wrf day'''
    def __init__(self, filedir, beg_time, end_time):
        self.filedir = filedir
        self.beg_time = beg_time
        self.end_time = end_time

if __name__ == '__main__':
    Storm = 'MARIA'
    Exper_name = 'newWRF_MW_THO'
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/' 

    DFtimes = ['201709170600','201709171200','201709171800',]
    DFend = '201709210000'

    USAF_list = glob.glob( small_dir + Storm + '/USAF/*.txt' )
    USAF_name_list = [ os.path.split(ifile)[1] for ifile in USAF_list ]

    for DF_start in DFtimes:
        print('The deterministic forecast begins at '+ DF_start)
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/wrf_df/'+DF_start+'/'
        plot_dir = wrf_dir + 'UV10_slp/'

        # calculate days between the DF_start and the DF_end
        time_diff_days = datetime.strptime(DFend[0:8],"%Y%m%d") - datetime.strptime(DF_start[0:8],"%Y%m%d")
        times_between_dt = [datetime.strptime(DF_start[0:8],"%Y%m%d") + timedelta(days=t) for t in list(range(0, time_diff_days.days+1, 1))]
        times_between_str = [datetime.strftime(it,"%Y%m%d") for it in times_between_dt]
        print( 'Days between the start and the end of the deterministic forecast: ', times_between_str )

        # ---- check if USAF obs exists for days within this period and collect it if it exists -----
        print('Check if USAF obs exists....')
        USAF_per_day_v1 = {}
        # simply check by examing the filename
        for it_str in times_between_str:
            USAF_name_exist = []

            if_USAF_exist = [] # if_USAF_exist: collect the USAF file for each day_DF
            for USAF_file in USAF_name_list:
                    if_USAF_exist.append( it_str in USAF_file ) # if_USAF_exist: [False, True

            if True in if_USAF_exist:
                idx_USAF_name = np.where(if_USAF_exist)[0]
                USAF_list_filter = sorted([USAF_name_list[it] for it in idx_USAF_name])
                for ivalue in USAF_list_filter:
                     USAF_name_exist += small_dir + Storm + '/USAF/' + ivalue,
                USAF_per_day_v1[it_str] = USAF_name_exist
                # # add a trailing comma to the end to allow looping the file when only one file exists
            else:
                USAF_per_day_v1[it_str] = None
        print('By examining the file name, USAF data exist on: ', USAF_per_day_v1)

        # further check by examing the content of each file
        USAF_per_day = {} #copy.deepcopy( USAF_per_day_v1 ) #shallow copy only creates a new refernce !!
        for it in USAF_per_day_v1.keys():
            USAF_per_day[it] =  []

        print('Further check if the file content spans two days...')
        for it in USAF_per_day_v1.keys():

            print('On the day ', it)
            # if no obs on this day
            if USAF_per_day_v1[it] is None:
                USAF_per_day[it] = None
                continue
            
            # ------- if obs exists on this day ---------------
            for idx in range(0,len(USAF_per_day_v1[it])):
                idx_t0 = times_between_str.index(it)
                t1_key = times_between_str[idx_t0+1]
                if_span_days,info =  USAF.check_timespan( USAF_per_day_v1[it][idx],it,t1_key )

                if if_span_days:
                    info_usaf_t0 = Info_USAF_clt( USAF_per_day_v1[it][idx],info[it][0], info[it][1] )
                    info_usaf_t1 = Info_USAF_clt( USAF_per_day_v1[it][idx],info[t1_key][0], info[t1_key][1] )
                    USAF_per_day[it].append( info_usaf_t0 )
                    USAF_per_day[t1_key].append( info_usaf_t1 )

                else: # obs doesn't span two days
                    info_usaf = Info_USAF_clt( USAF_per_day_v1[it][idx],info[it][0], info[it][1] )
                    USAF_per_day[it].append( info_usaf )
        
        print('End the further check!')

        # Enter the main function
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            plot_d03( wrf_dir, plot_dir, USAF_per_day )
        else:
            plot_d03( wrf_dir, plot_dir, USAF_per_day )


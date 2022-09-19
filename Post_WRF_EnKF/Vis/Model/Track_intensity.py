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

matplotlib.rcParams['lines.linewidth'] = 1
matplotlib.rcParams['lines.markersize'] = 2
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['font.size'] = 6

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True


# Read the post-storm best-track file of the storm at all times
def read_bestrack(Storm):

    Best_track_file = os.listdir('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/' + Storm + '/Post_Storm_btk')
    with open('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/' + Storm + '/Post_Storm_btk/' + Best_track_file[0]) as f:
        all_lines = f.readlines()
        
    # Process all of records to our format/unit 
    time_all = []
    lat_all = []
    lon_all = []
    maxV_all = []
    minP_all = []
    for line in all_lines:
        # split one record into different parts
        split_line = line.split()
        # Read time
        time_all.append(split_line[2].replace(',','') + '00')
        # Read latitude
        lat_line = split_line[6].replace(',','')
        if 'N' in lat_line:
            lat_all.append(float(lat_line.replace('N',''))/10)
        else:
            lat_all.append(0-float(lat_line.replace('S',''))/10)
        # Read longitute
        lon_line = split_line[7].replace(',','')
        if 'W' in lon_line:
            lon_all.append(0-float(lon_line.replace('W',''))/10)
        else:
            lon_all.append(float(lon_line.replace('E',''))/10)
        # Read max wind
        maxV_all.append(float(split_line[8].replace(',',''))*0.51444) # knots to m/s
        # Read min sea level pressure
        minP_all.append(float(split_line[9].replace(',',''))) # mb

    dict_bestrack = {'time': time_all, 'lat': lat_all, 'lon': lon_all, 'max_ws': maxV_all, 'min_slp': minP_all}
    return dict_bestrack


# Return records in the complete best-track file at times of interest 
def btk_in_duration(Storm, DF_start, DF_end):
    
    # Calculate the duration from the start to the end of deterministic forecast in hour
    DF_diff = datetime.strptime(DF_end,"%Y%m%d%H%M") - datetime.strptime(DF_start,"%Y%m%d%H%M") 
    DF_diff_hour = DF_diff.total_seconds() / 3600 
    # List the times (in string format) every 6 hours in the duration 
    time_interest_dt = [datetime.strptime(DF_start,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(DF_diff_hour), 6))]
    time_interest_str = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    # Read the complete best-track file
    dict_all_btk = read_bestrack(Storm)    
    # Get the indices in the best-track file corresponded to the times of interest
    idx_in_btk = []
    
    lat_btk = []
    lon_btk = []
    max_ws_btk = []
    min_slp_btk = []
    for time_str in time_interest_str:
        boolean_compare = [t_btk == time_str for t_btk in dict_all_btk['time'][:]]
        if any(boolean_compare):
            idx = int(np.where(boolean_compare)[0][0])
            idx_in_btk.append(idx) 
            lat_btk.append(dict_all_btk['lat'][idx])
            lon_btk.append(dict_all_btk['lon'][idx])
            max_ws_btk.append(dict_all_btk['max_ws'][idx])
            min_slp_btk.append(dict_all_btk['min_slp'][idx])

    dict_btk = {'time': time_interest_str, 'lat': lat_btk, 'lon': lon_btk, 'max_ws': max_ws_btk, 'min_slp': min_slp_btk}
    
    return dict_btk


# Get track/intensity from model output
def read_wrfout(directory, filename_pickle=None, force_reload=False):
    if filename_pickle is None:
        filename_pickle = directory + '/HPI.pickle'
    if not os.path.exists( filename_pickle ) or force_reload:
        filenames = sorted(glob.glob(directory + '/wrfout_d03_*') )
        nstep = len(filenames)
        HPI = {}
        HPI['time'] = [datetime(2017,1,1) for i in range(nstep)]
        HPI['max_ws'] = np.zeros( nstep )
        HPI['min_slp'] = np.zeros( nstep )
        HPI['lat'] = np.zeros( nstep )
        HPI['lon'] = np.zeros( nstep )
        for ifile in range(nstep):
            filename = filenames[ifile]
            print(filename)
            with Dataset(filename) as ncid:
                start_time = datetime.strptime( ncid.SIMULATION_START_DATE, '%Y-%m-%d_%H:%M:%S')
                dtime = timedelta( minutes= float(ncid.variables['XTIME'][:][0]) )
                time_file = start_time + dtime
                HPI['time'][ifile] = time_file.strftime("%Y%m%d%H%M") 
                ws = np.sqrt( ncid.variables['U10'][:]**2 + ncid.variables['V10'][:]**2 )
                HPI['max_ws'][ifile] = np.max( ws )
                slp = getvar(ncid, 'slp')
                HPI['min_slp'][ifile] = np.min( slp )
                
                # smooth to find min slp location
                slp_smooth = sp.ndimage.filters.gaussian_filter(slp, [11, 11] )
                idx = np.nanargmin( slp_smooth )
                HPI['lat'][ifile] = ncid.variables['XLAT'][:].flatten()[idx]
                HPI['lon'][ifile] = ncid.variables['XLONG'][:].flatten()[idx]

        # write to pickle
        with open(filename_pickle, 'wb') as f:
            pickle.dump(HPI,f)
    else:
        with open(filename_pickle, 'rb') as f:
            HPI = pickle.load(f)
    return HPI


def plot_one( ax0, ax1, ax2,  state, linestyle, label, step=1):
    dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in state['time']]
    #hours = np.array([date.hour for date in dates])
    #idx = np.mod(hours.astype(int), step) == 0
    ax0.plot(state['lon'], state['lat'], linestyle, label=label, transform=ccrs.PlateCarree())
    ax1.plot_date(dates, state['min_slp'], linestyle, label=label)
    ax2.plot_date(dates, state['max_ws'], linestyle, label=label)


def plot_hpi(Strom, domain_range, output_dir=None):
    reload_data = False
    step = 1

    # Set up figure
    fig = plt.figure( figsize=(12,4), dpi=150 )
    gs = fig.add_gridspec(1,3)

    ax0 = fig.add_subplot( gs[0,0],  projection=ccrs.PlateCarree()) 
    ax0.set_extent( domain_range,  crs=ccrs.PlateCarree())
    ax0.coastlines( resolution='10m', color='black',linewidth=0.5 )
    ax1 = fig.add_subplot( gs[0,1] )
    ax2 = fig.add_subplot( gs[0,2] )

<<<<<<< HEAD:Post_WRF_EnKF/Vis/Model/Track_intensity_IRMW.py
    # Plot HPI from post-storm analysis
    Btk_start = '201709030600'#'201708221200' #'201709030600'
    Btk_end = '201709090000'#'201708270000'#'201709090000'
    best_track = btk_in_duration(Storm, Btk_start, Btk_end)
    plot_one ( ax0, ax1, ax2, best_track,  'k-', 'Best track', step=step )

    # Plot HPI from deterministic forecasts
    DF_model_end  = '201709090000' #''201708270000'
    DF_init_times = sorted(os.listdir(  wrf_dir+'/'+Storm+'/'+Exper_name+'/wrf_df/' ))

    color_model = [ '#b8d5cd', '#8abaae', '#5ca08e', '#2e856e', '#006a4e']

    if read_HPI_wrfout == True:
        i = 0
        for init_time in DF_init_times:
            plot_one( ax0, ax1, ax2, read_wrfout(Storm, Exper_name, wrf_dir+'/'+Storm+'/'+Exper_name+'/wrf_df/'+init_time, force_reload=reload_data), color_model[i], 'IR+MW: '+ init_time, step=step)
            i = i+1
    elif read_HPI_wrfout == False:
        i = 0
        for init_time in DF_init_times:
            plot_one( ax0, ax1, ax2, read_rsl_error(Storm, Exper_name, wrf_dir+'/'+Storm+'/'+Exper_name+'/wrf_df/'+init_time, init_time, DF_model_end), color_model[i], 'IR+MW: '+ init_time, step=step)
            i = i+1
=======
    # Plot
    plot_one( ax0, ax1, ax2, read_wrfout('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/HARVEY/MW_THO/wrf_df/201708221200/mean/', force_reload=reload_data), 'b-', 'IR+MW: 201708221200', step=step)
    plot_one( ax0, ax1, ax2, read_wrfout('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/HARVEY/MW_THO/wrf_df/201708221800/mean/', force_reload=reload_data), 'g-', 'IR+MW: 201708221800', step=step)

    DF_start = '201708221200'
    DF_end = '201708270000'
    best_track = btk_in_duration(Storm, DF_start, DF_end)
    plot_one ( ax0, ax1, ax2, best_track,  'k-', 'Best track', step=step )
>>>>>>> parent of e20505b... Minor detail:Post_WRF_EnKF/Vis/Model/Track_intensity.py

    # Set labels
    lon_ticks = list(range(math.ceil(domain_range[0]), math.ceil(domain_range[1]), 2))
    lat_ticks = list(range(math.ceil(domain_range[2]), math.ceil(domain_range[3]), 2))
    gl = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
    gl.left_labels = True
    gl.bottom_label = True
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlabel_style = {'size': 6}
    gl.ylabel_style = {'size': 6}  
    
<<<<<<< HEAD:Post_WRF_EnKF/Vis/Model/Track_intensity_IRMW.py
    #ax1.set_xlim([datetime(2017, 8, 22, 12, 0, 0), datetime(2017, 8, 27)])
    #ax2.set_xlim([datetime(2017, 8, 22, 12, 0, 0), datetime(2017, 8, 27)])
    ax1.set_xlim([datetime(2017, 9, 3, 6, 0, 0), datetime(2017, 9, 9)])
    ax2.set_xlim([datetime(2017, 9, 3, 6, 0, 0), datetime(2017, 9, 9)])
=======
    ax1.set_xlim([datetime(2017, 8, 22, 12, 0, 0), datetime(2017, 8, 27)])
    ax2.set_xlim([datetime(2017, 8, 22, 12, 0, 0), datetime(2017, 8, 27)])
>>>>>>> parent of e20505b... Minor detail:Post_WRF_EnKF/Vis/Model/Track_intensity.py
    ax1.tick_params(axis='x', labelrotation=45)
    ax2.tick_params(axis='x', labelrotation=45)
    ax2.legend(frameon=False, loc='upper left')
    # Set titles
    ax0.set_title( 'Track' )
    ax1.set_title( 'MSLP' )
    ax2.set_title( 'Vmax' )
    
    if output_dir is not None:                             
<<<<<<< HEAD:Post_WRF_EnKF/Vis/Model/Track_intensity_IRMW.py
        plt.savefig('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/newWRF_MW_THO/Vis_analyze/Model/'+Storm+'_IRMW.png')
=======
        plt.savefig('./test.png')
>>>>>>> parent of e20505b... Minor detail:Post_WRF_EnKF/Vis/Model/Track_intensity.py
        #plt.savefig('%s/hpi_fcst_%s_%s.png' % (output_dir, '201708240000', '201708270000') )                                                
    else:
        plt.show()
    plt.clf()


if __name__ == '__main__':
<<<<<<< HEAD:Post_WRF_EnKF/Vis/Model/Track_intensity_IRMW.py
    Storm = 'IRMA'
    Exper_name = 'newWRF_MW_THO'
    wrf_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'

    read_HPI_wrfout = False

    lon_min = -85#-101#-85
    lon_max = -45#-86#-45
    lat_min = 12#16#12
    lat_max = 30#31#30
=======
    Storm = 'HARVEY'
    
    lon_min = -101
    lon_max = -85
    lat_min = 16
    lat_max = 31
>>>>>>> parent of e20505b... Minor detail:Post_WRF_EnKF/Vis/Model/Track_intensity.py
    domain_range = [lon_min, lon_max, lat_min, lat_max]

    plot_hpi( Storm, domain_range, output_dir = './out')

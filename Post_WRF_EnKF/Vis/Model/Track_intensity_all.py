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

matplotlib.rcParams['lines.linewidth'] = 2
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
def btk_in_duration(Storm, Btk_start, Btk_end):
    
    # Calculate the duration from the start to the end of deterministic forecast in hour
    Btk_diff = datetime.strptime(Btk_end,"%Y%m%d%H%M") - datetime.strptime(Btk_start,"%Y%m%d%H%M") 
    Btk_diff_hour = Btk_diff.total_seconds() / 3600 
    # List the times (in string format) every 6 hours in the duration 
    time_interest_dt = [datetime.strptime(Btk_start,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(Btk_diff_hour), 6))]
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
def read_wrfout(Storm, Exper_name, directory, filename_pickle=None, force_reload=False):
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


# Get track and intensity from ATCF part in rsl.error.0000
def read_rsl_error(Storm, Exper_name, directory, DF_start, DF_end):
    with open(directory + '/ATCF_rsl.error.0000') as f:
        all_lines = f.readlines()

    # Read and process records
    DF_start_real = [] # if use the option "time_to_move" in WRF, DF_start_real is not equal to DF_start
    time_all = []
    lat_all = []
    lon_all = []
    maxV_all = []
    minP_all = []

    num_line = 0
    for line in all_lines:
        # Split 
        split_line = line.split()
        # time
        time_wrf = split_line[1]
        time_dt = datetime.strptime( time_wrf,"%Y-%m-%d_%H:%M:%S" )
        time_all.append( time_dt.strftime( "%Y%m%d%H%M" ) )
        if num_line == 0:
            DF_start_real.append( time_dt.strftime( "%Y%m%d%H%M" ) )

        # lat
        lat_all.append( float(split_line[2]) )
        # lon
        lon_all.append( float(split_line[3]) )
        # minimum slp
        minP_all.append( float(split_line[4]) )
        # maximum wind speed
        maxV_all.append( float(split_line[5])*0.51444 )
        num_line = num_line + 1

    # Only select a subset of HPI records
    # Calculate the duration from the start to the end of deterministic forecast in hour
    DF_start = DF_start_real[0]
    DF_diff = datetime.strptime(DF_end,"%Y%m%d%H%M") - datetime.strptime(DF_start,"%Y%m%d%H%M")
    DF_diff_hour = DF_diff.total_seconds() / 3600
    # List the times (in string format) every hour in the duration 
    time_interest_dt = [datetime.strptime(DF_start,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(DF_diff_hour), 1))]
    time_interest_str = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    # Get the indices in the best-track file corresponded to the times of interest
    idx_sbs = []

    lat_sbs = []
    lon_sbs = []
    max_ws_sbs = []
    min_slp_sbs = []
    for time_str in time_interest_str:
        boolean_compare = [ eachT  == time_str for eachT in time_all ]
        if any(boolean_compare):
            idx = int( np.where(boolean_compare)[0] )
            idx_sbs.append(idx)
            lat_sbs.append( lat_all[idx] )
            lon_sbs.append( lon_all[idx] )
            max_ws_sbs.append( maxV_all[idx] )
            min_slp_sbs.append( minP_all[idx] )

    dict_model = {'time': time_interest_str, 'lat': lat_sbs, 'lon': lon_sbs, 'max_ws': max_ws_sbs, 'min_slp': min_slp_sbs} 

    return dict_model


def plot_one( ax0, ax1, ax2,  state, linestyle, label, step=1):
    dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in state['time']]
    #hours = np.array([date.hour for date in dates])
    #idx = np.mod(hours.astype(int), step) == 0
    ax0.plot(state['lon'], state['lat'], linestyle, label=label, transform=ccrs.PlateCarree())
    ax1.plot_date(dates, state['min_slp'], linestyle, label=label)
    ax2.plot_date(dates, state['max_ws'], linestyle, label=label)


def plot_hpi(Storm, wrf_dir, read_HPI_wrfout, domain_range, output_dir=None):
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

    # Plot HPI from post-storm analysis
    Btk_start = '201709030600'
    Btk_end = '201709090000'
    best_track = btk_in_duration(Storm, Btk_start, Btk_end)
    plot_one ( ax0, ax1, ax2, best_track,  'k-', 'Best track', step=step )

    # Plot HPI from deterministic forecasts
    DF_model_end  = '201709090000' #''201708270000'
    IR_init_times = sorted(os.listdir(  wrf_dir+'/'+Storm+'/newWRF_IR_only/wrf_df/' ))
    IRMW_init_times = sorted(os.listdir(  wrf_dir+'/'+Storm+'/newWRF_MW_THO/wrf_df/' ))
    
    IR_color = [ '#FFB5A6','#FF998B','FF7E72','E36359','C44841','A62C2B']
    IRMW_color = [ '#b8d5cd', '#8abaae', '#5ca08e', '2e856e', '006a4e']

    if read_HPI_wrfout == True:
        i = 0
        for init_time in DF_init_times:
            plot_one( ax0, ax1, ax2, read_wrfout(Storm, wrf_dir+'/'+Storm+'/newWRF_MW_THO/wrf_df/'+init_time, force_reload=reload_data), color_model[i], 'IR+MW: '+ init_time, step=step)
            i = i+1

    elif read_HPI_wrfout == False:
        i = 0
        for init_time in IR_init_times:
            plot_one( ax0, ax1, ax2, read_rsl_error(Storm, 'newWRF_IR_only', wrf_dir+'/'+Storm+'/newWRF_IR_only/wrf_df/'+init_time, init_time, DF_model_end), IR_color[i], 'IR: '+ init_time, step=step)
            i = i+1
        
        i = 0
        for init_time in IRMW_init_times:
            plot_one( ax0, ax1, ax2, read_rsl_error(Storm, 'newWRF_MW_THO', wrf_dir+'/'+Storm+'/newWRF_MW_THO/wrf_df/'+init_time, init_time, DF_model_end), IRMW_color[i], 'IR+MW: '+ init_time, step=step)
            i = i+1 



    # Set labels
    lon_ticks = list(range(math.ceil(domain_range[0]), math.ceil(domain_range[1]), 4))
    lat_ticks = list(range(math.ceil(domain_range[2]), math.ceil(domain_range[3]), 2))
    gl = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
    gl.ylabels_left = True
    gl.xlabels_bottom = True
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlabel_style = {'size': 6}
    gl.ylabel_style = {'size': 6}  
    
    ax1.set_xlim([datetime(2017, 9, 3, 6, 0, 0), datetime(2017, 9, 9)])
    ax2.set_xlim([datetime(2017, 9, 3, 6, 0, 0), datetime(2017, 9, 9)])
    ax1.tick_params(axis='x', labelrotation=45)
    ax2.tick_params(axis='x', labelrotation=45)
    ax2.legend(frameon=False, loc='upper left')
    # Set titles
    ax0.set_title( 'Track' )
    ax1.set_title( 'MSLP' )
    ax2.set_title( 'Vmax' )
    
    if output_dir is not None:                             
        plt.savefig('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/newWRF_MW_THO/Vis_analyze/Model/'+Storm+'_all.png')
        #plt.savefig('%s/hpi_fcst_%s_%s.png' % (output_dir, '201708240000', '201708270000') )                                                
    else:
        plt.show()
    plt.clf()


if __name__ == '__main__':
    Storm = 'IRMA'
    wrf_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'

    read_HPI_wrfout = False

    lon_min = -85
    lon_max = -45
    lat_min = 12
    lat_max = 30
    domain_range = [lon_min, lon_max, lat_min, lat_max]

    plot_hpi( Storm, wrf_dir, read_HPI_wrfout, domain_range, output_dir = './out')

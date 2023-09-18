#!/bin/bash
# Credit to Yinghui LV 
# Author: Zhu (Judy) Yao. 2022-2023.


# ------------------------------------------------------------------------------------------------------
#               Import Modules 
# ------------------------------------------------------------------------------------------------------
# modules to interact with the OS (operating system)
import os,fnmatch 
import glob
# modules to perform scientific calculations
import math
import numpy as np
import scipy as sp
from datetime import datetime, timedelta
# modules to process netcdf files (such as WRF output)
from netCDF4 import Dataset
from wrf import getvar 
    # It might be possible that you are not able to conda install wrf-var with a pretty new python version
    # Solution:
    # 1. conda create -n $PYTHON34_ENV_NAME python=3.4 anaconda 
    # 2. conda activate python=3.4 (use wrf-python in this python environment)
# modules to plot figures
import matplotlib
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs # map drawing
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
# ------------------------------------------------------------------------------------------------------

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['lines.linewidth'] = 2.5#1.5
matplotlib.rcParams['lines.markersize'] = 2.5
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['font.size'] = 15#6

# ------------------------------------------------------------------------------------------------------
#           Object: Post-storm analysis from National Hurricane Center (NHC)
# ------------------------------------------------------------------------------------------------------

# Read the post-storm best-track file of the storm at all times
def read_bestrack( Storm ):

    """
    This function reads out all records in the post-storm best-track file of the storm analyzed by NHC

    Args:
        Storm (:obj:`str`): name of the storm. E.g., 'MARIA'
    
    Returns:
        A dictionary.
    """

    # Identify the file
    btk_dir = small_dir + Storm + '/Post_Storm_btk/'
    Best_track_file = os.listdir( btk_dir )

    # Read all lines in the file into the memory
    # Note: records are only available at synoptic time (00,06,12,18UTC)
    with open( btk_dir + Best_track_file[0]) as f:
        all_lines = f.readlines()

    # Process all of records to formats/units used by us
    # ---------------------------------------------------------------------------- 
    time_all = []
    lat_all = []
    lon_all = []
    maxV_all = []
    minP_all = []
    # Loop thru each record
    for line in all_lines:
        # split one record into different parts/strings
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
        minP_all.append(float(split_line[9].replace(',',''))) # hPa

    # Assemble the dictionary
    # -------------------------
    dict_bestrack = {'time': time_all, 'lat': lat_all, 'lon': lon_all, 'max_ws': maxV_all, 'min_slp': minP_all}
    return dict_bestrack


# Read records only at times of interest from the complete best-track file
def btk_in_duration(Storm, Btk_start, Btk_end, hour_step):

    """
    This function reads and stores a subset of records from the complete post-storm best-track file

    Args:
        Storm (:obj:`str`): name of the storm. E.g., 'MARIA'
        Btk_start (:obj:`str`): start point in the time period of interest to this study. E.g., '201709160000'
        Btk_end (:obj:`str`): end point in the time period of interest to this study. E.g., '201709210000'
        hour_step (:obj:`float`): how often to subset the recors. E.g., 6 (hours)
 
    Returns:
        A dictionary.
    """

    # Calculate the duration from the start to the end of time period of interest, in hour
    Btk_diff = datetime.strptime(Btk_end,"%Y%m%d%H%M") - datetime.strptime(Btk_start,"%Y%m%d%H%M")
    Btk_diff_hour = Btk_diff.total_seconds() / 3600
    # List the times (in string format) every 6 hours in the duration 
    time_interest_dt = [datetime.strptime(Btk_start,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(Btk_diff_hour+1), hour_step))]
    time_interest_str = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    # Read the complete best-track file
    dict_all_btk = read_bestrack(Storm)

    # Reads out the subset of records corresponding to time of interest
    # --------------------------------------------------------------------
    idx_in_btk = [] # indices in the best-track file corresponded to the times of interest
    lat_btk = []
    lon_btk = []
    max_ws_btk = []
    min_slp_btk = []

    # Loop thru each time string of interest
    for time_str in time_interest_str:
        boolean_compare = [t_btk == time_str for t_btk in dict_all_btk['time'][:]]
        if any(boolean_compare):
            idx = int(np.where(boolean_compare)[0][0])
            idx_in_btk.append(idx)
            lat_btk.append(dict_all_btk['lat'][idx])
            lon_btk.append(dict_all_btk['lon'][idx])
            max_ws_btk.append(dict_all_btk['max_ws'][idx])
            min_slp_btk.append(dict_all_btk['min_slp'][idx])

    # Assemble the dictionary
    # -------------------------
    dict_btk = {'time': time_interest_str, 'lat': lat_btk, 'lon': lon_btk, 'max_ws': max_ws_btk, 'min_slp': min_slp_btk}

    return dict_btk


# ------------------------------------------------------------------------------------------------------
#            Object: Deterministic Forecast; Operation: Read and Process  
# ------------------------------------------------------------------------------------------------------

# Read out a storm forecast's HPI (Hurricane Potential Intensity) from ATCF part in rsl.error.0000
def read_rsl_error(Storm, Exper_name, directory, DF_start, DF_end):
    
    """
    This function reads out hourly HPI of the forecast within the time period of interest. Before you use this function,
    make sure you have collected data from rsl.error.0000 file. 
        E.g., issue 'grep ATCF rsl.erro.0000 > ATCF_rsl.error.0000' on your terminal window

    Args:
        Storm (:obj:`str`): name of the storm. E.g., 'MARIA'
        Exper_name (:obj:`str`): name of the experiment.  E.g., 'IR_THO'
        directory (:obj:`str`): absolute path to the directory where the rsl.error.0000 resides.
        DF_start (:obj:`str`): initialization time of the forecast.Only available at 00, 06, 12, 18 UTC.  E.g., '201709160000'
        DF_end (:obj:`str`): end time of the forecast. E.g., '201709210000'
    
    Returns:
        A dictionary.
    """

    # Check if the file exists & If so, read all lines in the file into the memory
    # ----------------------------------------------------------------------------
    if os.path.exists( directory + '/ATCF_rsl.error.0000' ):
        with open(directory + '/ATCF_rsl.error.0000') as f:
            all_lines = f.readlines()
    else:
        return None        

    # Read and process each line/record
    # Note: HPI of the forecast is written to the file every 15 minutues (in forecast time) 
    # -------------------------------------------------------------------------------------
    DF_start_real = []  # if the value of "time_to_move" in the namelist.input  is not set as 0 in WRF, DF_start_real will not be equal to DF_start.
                        # E.g., Set "time_to_move = 60,60,60" and the initialization time (DF_start) of the experiment is at '201709160000'.
                        #       the WRF program will not output HPI data to rsl.error.0000 until '201709160100'.
                        #       In this case, DF_start_real is '201709160100'.
    time_all = []
    lat_all = []
    lon_all = []
    maxV_all = []
    minP_all = []

    num_line = 0
    # Loop thru each line
    for line in all_lines:
        # Split the line into a list of strings 
        split_line = line.split()
        # time
        time_wrf = split_line[1]
        time_dt = datetime.strptime( time_wrf,"%Y-%m-%d_%H:%M:%S" ) #2017-09-16_03:00:00
        time_all.append( time_dt.strftime( "%Y%m%d%H%M" ) ) #201709160300
        if num_line == 0:
            DF_start_real.append( time_dt.strftime( "%Y%m%d%H%M" ) )
        # lat
        lat_all.append( float(split_line[2]) )
        # lon
        lon_all.append( float(split_line[3]) )
        # minimum slp (hPa)
        minP_all.append( float(split_line[4]) )
        # maximum wind speed (m/s)
        maxV_all.append( float(split_line[5])*0.51444 )
        num_line = num_line + 1

    # Only reads out hourly records (namely, one of every four records)
    # --------------------------------------------------------------------
    # Calculate the distance from the real start point to the end of the forecast in hour
    DF_diff = datetime.strptime(DF_end,"%Y%m%d%H%M") - datetime.strptime(DF_start_real[0],"%Y%m%d%H%M")
    DF_diff_hour = DF_diff.total_seconds() / 3600
    # List these times (in string format) every hour
    time_interest_dt = [datetime.strptime(DF_start_real[0],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(DF_diff_hour), 1))]
    time_interest_str = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Reads out the subset of records corresponding to time of interest
    idx_sbs = []
    lat_sbs = []
    lon_sbs = []
    max_ws_sbs = []
    min_slp_sbs = []

    # Loop thru each time string of interest
    for time_str in time_interest_str:
        # identify the location/index of the time of interest in the ATCF_rsl.error.0000
        boolean_compare = [ eachT  == time_str for eachT in time_all ]
        if any(boolean_compare):
            if len( np.where(boolean_compare)[0] ) > 1:
                idx = int( np.where(boolean_compare)[0][0] )
            else:
                idx = int( np.where(boolean_compare)[0] )
            # Store HPI at this time
            idx_sbs.append(idx)
            lat_sbs.append( lat_all[idx] )
            lon_sbs.append( lon_all[idx] )
            max_ws_sbs.append( maxV_all[idx] )
            min_slp_sbs.append( minP_all[idx] )

    # Calculate the distance between the next 6x time and the DF_start_real
    # ---------------------------------- idea ---------------------------------- 
    # This piece of info helps identify if time_to_move is set as 0 or not without reading namelist.input.
    # E.g., the initialization time (DF_start) is 201709160000, the first ATCF record is at 201709160100 (DF_start_real).
    #       We can infer the next synoptic time (00,06,12,18UTC) from the initialization time is 201709160600 (start_next6X).
    #       If the difference between start_next6X and the DF_start_real is between 0 and 6 (i.e.,1,2,3,4,5),
    #       it means the time_to_move option is not 0 (ie., DF_start is not equal to DF_start_real).
    #       Under such a circumstance, the plotting program should be tweaked specially.
    
    start_next6X = datetime.strptime(DF_start,"%Y%m%d%H%M") + timedelta(hours=6)
    Diff_start_next6x = (start_next6X - datetime.strptime(DF_start_real[0],"%Y%m%d%H%M")).total_seconds() / 3600

    # Assemble the dictionary
    # -------------------------
    dict_model = {'Diff_start_next6x': int(Diff_start_next6x), 'time': time_interest_str, 'lat': lat_sbs, 'lon': lon_sbs, 'max_ws': max_ws_sbs, 'min_slp': min_slp_sbs} 
    return dict_model


# ------------------------------------------------------------------------------------------------------
#            Object: Deterministic Forecast; Operation: Plot
# ------------------------------------------------------------------------------------------------------
# Plot one set of data
def plot_one( ax0, ax1, ax2, state, color, line, line_width, label, steps=6):

    """
    This function plot one set of data (E.g., best-track data or one forecast data) on the figure. 
    Multiple calling of this function results in all of these data visualized on the figure.

    Args:
        ax0/ax1/ax2 (:obj:`Axes`): an individual plot within a figure.
        state (:obj:`dictionary`): HPI data read from either the best-track data or forecast data.   
        color (:obj:`Hexadecimal RGB values`): line color.
        line : line style.
        line_width (:obj:`float`): line width.
        label (:obj:`string`): indicates if the data is best-track or forecast.
        steps (:obj:`int`): how often to plot the data. E.g., every 6 hours
    
    Returns:
        A finished individual plot within the figure.
    """
 
    if label == 'Best track':
        times = state['time']
        lon = state['lon']
        lat = state['lat']
        x_min_slp = state['min_slp']
        x_max_ws = state['max_ws']

        ax0.plot(lon, lat, marker='o',markersize=3, color=color,linewidth=3, label=label, linestyle=line, transform=ccrs.PlateCarree())
        dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]
        ax1.plot_date(dates, x_min_slp, color, label=label, linewidth=3)
        ax2.plot_date(dates, x_max_ws, color, label=label, linewidth=3)

        for it in times:
            if it[8:10] == '00':
                boolean_compare = [ it  == eachT for eachT in times ]
                idx = int( np.where(boolean_compare)[0] )
                ax0.scatter( lon[idx], lat[idx],s=5, marker='o',edgecolor="white",transform=ccrs.PlateCarree())
                ax0.annotate(it[6:8], xy=(lon[idx], lat[idx]), color=color, xycoords='data', transform=ccrs.PlateCarree())
        
    else:
        Diff_start_next6x = state['Diff_start_next6x']
    
        if Diff_start_next6x != 6 and Diff_start_next6x != 0: # it means time_to_move in namelist is not equal to 0
            start_next6x = datetime.strptime(state['time'][0],"%Y%m%d%H%M") + timedelta(hours=Diff_start_next6x)
            boolean_compare = [ start_next6x.strftime("%Y%m%d%H%M") == time_str for time_str in state['time'][:] ]
            idx_next6x = int( np.where(boolean_compare)[0] )

            times = state['time'][idx_next6x::steps]
            lon = state['lon'][idx_next6x::steps]
            lat = state['lat'][idx_next6x::steps]
            x_min_slp = state['min_slp'][idx_next6x::steps]
            x_max_ws = state['max_ws'][idx_next6x::steps]
            
            times.insert(0,state['time'][0])
            lon.insert(0,state['lon'][0])
            lat.insert(0,state['lat'][0])
            x_min_slp.insert(0,state['min_slp'][0])
            x_max_ws.insert(0,state['max_ws'][0])
        else:
            times = state['time'][::steps]
            lon = state['lon'][::steps]
            lat = state['lat'][::steps]
            x_min_slp = state['min_slp'][::steps]
            x_max_ws = state['max_ws'][::steps]
        
        ax0.plot(lon, lat, marker='o', markersize=3, color=color,linewidth=line_width, label=label, linestyle=line, transform=ccrs.PlateCarree())
        for it in times:
            if it[8:10] == '00':
                boolean_compare = [ it  == eachT for eachT in times ]
                idx = int( np.where(boolean_compare)[0] )
                ax0.scatter( lon[idx], lat[idx],s=5, marker='o',edgecolor="white",transform=ccrs.PlateCarree())    
                ax0.annotate(it[6:8], xy=(lon[idx], lat[idx]), color=color, xycoords='data', transform=ccrs.PlateCarree())

        dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times] 
        ax1.plot_date(dates, x_min_slp, color, label=label, linestyle=line)
        ax2.plot_date(dates, x_max_ws, color, label=label, linestyle=line)


# Plot track, min sea level pressure, maximum wind speed of the storm in a forecast
def plot_hpi_df( wrf_dir, Storm, Exper, domain_range ):

    """
    This function identifies the available forecasts of one or multiple experiments and visualize them. 

    Args:
        wrf_dir (:obj:`str`): parent directory of the experiments
        Storm (:obj:`str`): name of the storm. E.g., 'MARIA'   
        Exper (:obj:`str`): a list of experiment names.
        domain_range (:obj:`float`): domain range to plot for track visualization.
 
    Returns:
        A figure saved on the specified directory.
    """

    # Set the end point of the period to plot
    # -----------------------------------------------------------------
    if Storm == 'JOSE':
        DF_model_end  = '201709100000'
    else:
        DF_model_end  = None

    # Identify the available forecasts with their initialization time as keywords
    # ---------------------------------------------------------------------------
    Exper_content_lbl = {}
    for iExper in Exper:
        if iExper is not None:
            if os.path.exists( wrf_dir+'/'+Storm+'/'+iExper+'/wrf_df/' ):
                Exper_content_lbl[iExper] = sorted(fnmatch.filter(os.listdir( wrf_dir+'/'+Storm+'/'+iExper+'/wrf_df/' ),'20*')) 
                # manually set 
                #Exper_content_lbl[iExper] = ['201709160600','201709161200','201709161800','201709170000','201709170600'] 
            else:
                Exper_content_lbl[iExper] = None
        else:
             Exper_content_lbl[iExper] = None

    # Set up figure
    # -----------------------------------------------------------------
    # Create a figure
    fig = plt.figure( figsize=(22,6.5), dpi=150 ) 
    # Control the layout of all plot elements 
    gs = fig.add_gridspec(2,7)
    # Set up the layout for Track plot
    ax0 = fig.add_subplot( gs[0:,0:3],  projection=ccrs.PlateCarree())
    ax0.set_extent( domain_range,  crs=ccrs.PlateCarree())
    ax0.coastlines( resolution='10m', color='black',linewidth=0.5 )
    # Set up the layout for Intensity plots
    ax1 = fig.add_subplot( gs[:,3:5] )
    ax2 = fig.add_subplot( gs[:,5:7] )

    # Plot best-track data/post-storm analysis from National Hurricane Center (NHC)
    # -----------------------------------------------------------------
    # Set start and end of the period for plotting
    if Storm == 'JOSE':
        Btk_start = '201709050000'
        Btk_end = '201709100000'
    else:
        pass
    
    # Plot HPI from post-storm analysis
    best_track = btk_in_duration(Storm, Btk_start, Btk_end, hour_step=6)
    plot_one ( ax0, ax1, ax2, best_track,  'black', '-', 3, 'Best track' )

    # Plot forecasts of one or two experiments
    # -----------------------------------------------------------------
    # Customize color maps for forecasts of different experiments (only support <= 2 experiments)
    Color1 = ["#c23728","#e14b31","#de6e56","#e1a692","#786028","#a57c1b","#d2980d","#ffb400","#503f3f","#6d4b4b","#a86464","#e27c7c"] #redish
    Color2 = ["#115f9a", "#1984c5", "#22a7f0", "#48b5c4", "#48446e", "#5e569b", "#776bcd", "#9080ff","#3c4e4b", "#466964", "#599e94", "#6cd4c5"] #blueish
    Color_set = {'c0':Color1, 'c1':Color2}
    # Customize linestyles
    Line_types = ['-','-']
   
    # Customize labels 
    Labels = ['THO:',] 

    # Plot HPI for each deterministic forecasts
    iExper = 0
    # Loop thru each experiments
    for key in Exper_content_lbl: 
        print('Experiment: ', key)
        if Exper_content_lbl[key] is not None:
            ic = 0
            # Loop thru each forecast of this experiment
            for it in Exper_content_lbl[key]:
                print(wrf_dir+'/'+Storm+'/'+key+'/wrf_df/'+it)
                HPI_model = read_rsl_error(Storm, key, wrf_dir+'/'+Storm+'/'+key+'/wrf_df/'+it, it, DF_model_end)
                plot_one( ax0, ax1, ax2, HPI_model, Color_set['c'+str(iExper)][ic], Line_types[iExper], 1.5, Labels[iExper]+it, steps=6 )
                ic = ic + 1
        else:
            print('No available data!')
        iExper = iExper + 1

    # Set ticks/labels for track subplot
    lon_ticks = list(range(math.ceil(domain_range[0])-2, math.ceil(domain_range[1])+2, 4))
    lat_ticks = list(range(math.ceil(domain_range[2])-2, math.ceil(domain_range[3])+2, 2))
    gl = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
    gl.ylabels_left = True
    gl.xlabels_bottom = True
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlabel_style = {'size': 15}
    gl.ylabel_style = {'size': 15}  
   
    # Set ticks/labels for intensity subplots
    ax1.set_xlim([datetime(int(Btk_start[0:4]), int(Btk_start[4:6]), int(Btk_start[6:8]), int(Btk_start[8:10])), datetime(int(Btk_end[0:4]), int(Btk_end[4:6]), int(Btk_end[6:8]), int(Btk_end[8:10]))])
    ax2.set_xlim([datetime(int(Btk_start[0:4]), int(Btk_start[4:6]), int(Btk_start[6:8]), int(Btk_start[8:10])), datetime(int(Btk_end[0:4]), int(Btk_end[4:6]), int(Btk_end[6:8]), int(Btk_end[8:10]))])
    ax1.set_ylim([900,1020])     
    ax2.set_ylim([10,80])   
    ax1.tick_params(axis='x', labelrotation=30,labelsize=12)
    ax2.tick_params(axis='x', labelrotation=30,labelsize=12)
    ax2.legend(bbox_to_anchor=(1.5, 1.0),frameon=True,loc='upper right',fontsize='10')
    
    # Set titles
    ax0.set_title( 'Track',fontsize = 15 )
    ax1.set_title( 'MSLP (hPa)',fontsize = 15 )
    ax2.set_title( 'Vmax ($\mathregular{ms^{-1}}$)',fontsize = 15 )
    fig.suptitle('THO',fontsize = 15)

    # Save the figure
    des_name = Storm+'_forecast.png'
    #des_name = small_dir+Storm+'/'+Exper[0]+'/Vis_analyze/Model/'+Storm+'_forecast_IR_MW.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name+'!' )

if __name__ == '__main__':
    
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/' # where forecasts (input data) are 
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/' # where to output the figure

    # configuration
    Storm = 'JOSE'
    Exper_name = ['IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900',]

    # Pre-set the domain for track plotting
    if Storm == 'JOSE':
        lon_min = -65
        lon_max = -35
        lat_min = 0
        lat_max = 25
    else:
        pass

    domain_range = [lon_min, lon_max, lat_min, lat_max]

    # Call the function to plot 
    plot_hpi_df(  big_dir, Storm, Exper_name, domain_range )



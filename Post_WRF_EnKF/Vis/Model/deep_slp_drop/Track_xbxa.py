#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os,fnmatch # functions for interacting with the operating system
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
import matplotlib.dates as mdates
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from itertools import chain
import time

import Util_data as UD

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['lines.linewidth'] = 2.5#1.5
matplotlib.rcParams['lines.markersize'] = 2.5
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['font.size'] = 20#6


big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'


# ------------------------------------------------------------------------------------------------------
#           Object: Analyses (wrf_enkf_output_d03_mean); Operation: Read and Process  
# ------------------------------------------------------------------------------------------------------

# Get analyses' estimated HPI 
def read_HPI_model( Storm, Exper_name, file_kind, DAtimes ):

    DAtime_str = []
    max_wind = []
    min_slp = []
    lat_storm = []
    lon_storm = []
       
    # Read through analyses of the whole ensemble
    for it_dir in DAtimes:
        # get DA time
        DAtime_str.append( os.path.split(it_dir)[1])

        filename = it_dir + '/' + file_kind
        if not os.path.exists( filename ):
            continue
        with Dataset(filename) as ncid:
            #print('Reading ', filename)
            # maximum wind
            ws = np.sqrt( ncid.variables['U10'][:]**2 + ncid.variables['V10'][:]**2 )
            max_wind.append( np.max( ws ))
            # minimum sea level pressure
            slp = UD.compute_slp( ncid )
            min_slp.append( np.min( slp )) 
            slp_smooth = sp.ndimage.filters.gaussian_filter(slp, [11, 11] )
            idx = np.nanargmin( slp_smooth )
            lat_storm.append( ncid.variables['XLAT'][:].flatten()[idx]) 
            lon_storm.append( ncid.variables['XLONG'][:].flatten()[idx]) 
            
    HPI_model = {'time':DAtime_str, 'lat': lat_storm, 'lon':lon_storm, 'max_ws': max_wind, 'min_slp':min_slp}

    return HPI_model

# ------------------------------------------------------------------------------------------------------
#            Operation: Plot
# ------------------------------------------------------------------------------------------------------
def plot_one( ax0, ax1, ax2,  state, color, style, line_width, label, steps=1):
 
    if label == 'Best track':
        times = state['time']
        lon = state['lon']
        lat = state['lat']
        x_min_slp = state['min_slp']
        x_max_ws = state['max_ws']

        ax0.plot(lon, lat, marker='o',markersize=3, color=color,linewidth=3, label=label, linestyle=style, transform=ccrs.PlateCarree())
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
        # Plot track for xb and xa
        for key in state.keys():
            lon = state[key]['lon'][::steps]
            lat = state[key]['lat'][::steps]
            if key == 'wrf_enkf_input_d03_mean':
                ax0.plot(lon, lat, marker='o', markersize=3, color=color,linewidth=line_width, label='Xb'+label, linestyle='-', transform=ccrs.PlateCarree())
            else:
                ax0.plot(lon, lat, marker='o', markersize=3, color=color,linewidth=line_width, label='Xa'+label, linestyle='--', transform=ccrs.PlateCarree())
        # Plot zipped HPI of xb and xa
        times = state['wrf_enkf_input_d03_mean']['time'][::steps]
        dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]
        dates_zip = list( chain.from_iterable( zip(dates,dates)) )

        xb_min_slp = state['wrf_enkf_input_d03_mean']['min_slp'][::steps]
        xa_min_slp = state['wrf_enkf_output_d03_mean']['min_slp'][::steps]
        min_slp_zip = list( chain.from_iterable( zip(xb_min_slp,xa_min_slp)) )

        xb_max_ws = state['wrf_enkf_input_d03_mean']['max_ws'][::steps]
        xa_max_ws = state['wrf_enkf_output_d03_mean']['max_ws'][::steps]
        max_ws_zip = list( chain.from_iterable( zip(xb_max_ws,xa_max_ws)) )

        ax1.plot_date(dates_zip, min_slp_zip, color, linewidth=line_width )
        ax2.plot_date(dates_zip, max_ws_zip, color, label=label, linewidth=line_width )


def Plot_HPI_xbxa(  domain_range ):

    # Set the end point of the period to investigate
    if Storm == 'HARVEY':
        DF_model_end  = '201708251200'
    elif Storm == 'IRMA':
        DF_model_end  = '201709050000'
    elif Storm == 'MARIA':
        DF_model_end  = '201709180000'
    elif Storm == 'JOSE':
        DF_model_end  = '201709070000'
    else:
        DF_model_end  = None

    # Set up figure
    fig = plt.figure( figsize=(22,12), dpi=300 ) #12,4 20,6
    gs = fig.add_gridspec(4,7) # 2,7

    ax0 = fig.add_subplot( gs[:2,0:3],  projection=ccrs.PlateCarree())
    ax0.set_extent( domain_range,  crs=ccrs.PlateCarree())
    ax0.coastlines( resolution='10m', color='black',linewidth=0.5 )
    ax1 = fig.add_subplot( gs[:2,3:5] )
    ax2 = fig.add_subplot( gs[:2,5:7] )
    ax3 = fig.add_subplot( gs[2:,1:6] )

    # Set start and end of the period
    if Storm == 'HARVEY':
        Btk_start = '201708221200' # '201709161800' #'201709030600'
        Btk_end = '201708251200' # '201709210000' #'201709090000'
    elif Storm == 'IRMA':
        Btk_start = '201709030000'
        Btk_end = '201709050000'
    elif Storm == 'MARIA':
        Btk_start = '201709160000'#'201709160000'
        Btk_end = '201709180000'
    elif Storm == 'JOSE':
        Btk_start = '201709050000'
        Btk_end = '201709070000'
    else:
        pass

    # Plot HPI from post-storm analysis
    best_track = UD.btk_in_duration(Storm, Btk_start, Btk_end, hour_step=6)
    plot_one ( ax0, ax1, ax2, best_track,  'black', '-', 3, 'Best track' )

    Colors = ['#e6194B','#4363d8',] #'#748b97'
    Line_types = ['-','-']
    Labels = MP

    # Plot HPI from model
    for iExper in Exper_names:
        print('Plotting '+iExper)
        idx_exper = Exper_names.index( iExper )
        HPI_models = {}
        file_kinds = ['wrf_enkf_input_d03_mean','wrf_enkf_output_d03_mean']
        for ifk in file_kinds:
            idx = file_kinds.index( ifk )
            # Find the times of interest
            DAtimes_dir =  sorted(glob.glob(big_dir+Storm+'/'+Exper_name+'/fc/20*') )
            DAtimes_dir.pop(0) #remove the first directory (spin-up)
            HPI_models[ifk] = read_HPI_model( Storm, iExper, ifk, DAtimes_dir )
        plot_one( ax0, ax1, ax2, HPI_models, Colors[idx_exper], '-',  3.5, Labels[idx_exper], steps=1 )
        # Plot increment
        times = HPI_models['wrf_enkf_input_d03_mean']['time']
        dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]
        xb_min_slp = HPI_models['wrf_enkf_input_d03_mean']['min_slp']
        xa_min_slp = HPI_models['wrf_enkf_output_d03_mean']['min_slp']
        incre_slp = np.array(xa_min_slp) - np.array(xb_min_slp)
        #ax3.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        #ax3.xaxis.set_major_locator(mdates.DayLocator())
        #ax3.semilogy(dates, abs(incre_slp), color=Colors[idx_exper], linewidth=3, linestyle='-')
        ax3.plot_date(dates, incre_slp, color=Colors[idx_exper], linewidth=5, linestyle='-' )

    # Set ticks/labels for track subplot
    lon_ticks = list(range(math.ceil(domain_range[0])-2, math.ceil(domain_range[1])+2, 4))
    lat_ticks = list(range(math.ceil(domain_range[2])-2, math.ceil(domain_range[3])+2, 2))
    gl = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=0.5, color='gray', linestyle='--')
    gl.ylabels_left = True
    gl.xlabels_bottom = True
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlabel_style = {'size': 15}
    gl.ylabel_style = {'size': 20}  
   
    # Set ticks for intensity subplots
    ax_date = [ax1,ax2,ax3]
    for iax in ax_date:
        iax.set_xlim([datetime(int(Btk_start[0:4]), int(Btk_start[4:6]), int(Btk_start[6:8]), int(Btk_start[8:10])), datetime(int(Btk_end[0:4]), int(Btk_end[4:6]), int(Btk_end[6:8]), int(Btk_end[8:10]))])
   
    # Set legend
    ax2.legend(bbox_to_anchor=(1.45, 1.0),frameon=True,loc='upper right',fontsize='18')

    # Set limitation and ticks
    ax1.set_ylim([900,1020])     #([940, 1015])
    ax2.set_ylim([10,80])   #([10,60])
    ax3.set_ylim([10,-100])  
    ax1.tick_params(axis='x', labelrotation=12, labelsize=12)
    ax2.tick_params(axis='x', labelrotation=12,labelsize=12)
    ax3.tick_params(axis='x', labelrotation=12,labelsize=20)
    ax1.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax1.xaxis.set_minor_locator(mdates.HourLocator(interval=6))
    ax2.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax2.xaxis.set_minor_locator(mdates.HourLocator(interval=6))
    ax3.xaxis.set_major_locator(mdates.HourLocator(interval=6))
    ax3.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d-%H'))

    ax3.set_ylabel('Increment of MSLP (hPa)',fontsize = 20)
    ax1.grid(True,color='grey',linewidth=0.5, linestyle='--')
    ax2.grid(True,color='grey',linewidth=0.5, linestyle='--')
    ax3.grid(True,color='grey',linewidth=0.5, linestyle='--')

    # Set titles
    ax0.set_title( 'Track' )
    ax1.set_title( 'MSLP (hPa)' )
    ax2.set_title( 'Vmax ($\mathregular{ms^{-1}}$)')
    fig.suptitle( Storm+': '+DA )

    des_name = small_dir+Storm+'/'+Exper_names[0]+'/Vis_analyze/Model/'+Storm+'_'+DA+'_xbxa_HPI.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name+'!' )

if __name__ == '__main__':
    
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'JOSE'
    DA = 'IR'
    MP = ['THO','WSM6']

    # Pre-set the domain for DF forecast
    if Storm == 'HARVEY':
        lon_min = -100
        lon_max = -85
        lat_min = 16
        lat_max = 30
    elif Storm == 'IRMA':
        lon_min = -58
        lon_max = -44
        lat_min = 14
        lat_max = 24
    elif Storm == 'MARIA':
        lon_min = -58
        lon_max = -44#-45
        lat_min = 5# 
        lat_max = 15
    elif Storm == 'JOSE':
        lon_min = -48
        lon_max = -34#-45
        lat_min = 8
        lat_max = 18
    else:
        pass
    domain_range = [lon_min, lon_max, lat_min, lat_max]
    # -------------------------------------------------------  

    # Create experiment names
    Exper_names = []
    for imp in MP:
        Exper_names.append( UD.generate_one_name( Storm,DA,imp ) )  

    start_time=time.process_time()
    Plot_HPI_xbxa( domain_range )
    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')










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
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import time

import Util_data as UD

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['lines.linewidth'] = 2.5#1.5
matplotlib.rcParams['lines.markersize'] = 2.5
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['font.size'] = 15#6


# ------------------------------------------------------------------------------------------------------
#           Object: Analyses (wrf_enkf_output_d03_mean); Operation: Read and Process  
# ------------------------------------------------------------------------------------------------------

# Get analyses' estimated HPI 
def read_HPI_analyses(Storm, Exper_name, wrf_dir, filename_analyses=None, force_reload=False):
  

    #start_als = '201708231800'
    #end_als = '201708241200'
    #als_diff = datetime.strptime(end_als,"%Y%m%d%H%M") - datetime.strptime(start_als,"%Y%m%d%H%M")
    #als_diff_hour = als_diff.total_seconds() / 3600
    #time_interest_dt = [datetime.strptime(start_als,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(als_diff_hour), 1))]
    #DAtimes =  [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    #DAtimes_dir = [wrf_dir+Storm+'/'+Exper_name+'/fc/'+ it for it in DAtimes]
    
    DAtimes_dir =  sorted(glob.glob(wrf_dir+Storm+'/'+Exper_name+'/fc/20*') )
    # remove the first directory (spin-up)
    DAtimes_dir.pop(0)

    DAtime_str = []
    max_wind = []
    min_slp = []
    lat_storm = []
    lon_storm = []
       
    # Read through analyses of the whole ensemble
    for it_dir in DAtimes_dir:
        # get DA time
        DAtime_str.append( os.path.split(it_dir)[1])

        filename = it_dir+'/wrf_enkf_output_d03_mean'
        if not os.path.exists( filename ):
            continue
        #print('Reading analysis: '+ filename)
        with Dataset(filename) as ncid:
            print('Reading analysis: ', filename)
            # maximum wind
            ws = np.sqrt( ncid.variables['U10'][:]**2 + ncid.variables['V10'][:]**2 )
            max_wind.append( np.max( ws ))
            # minimum sea level pressure
            slp = getvar(ncid, 'slp')
            min_slp.append( np.min( slp )) 
            # location of the minimum slp
            slp_smooth = sp.ndimage.filters.gaussian_filter(slp, [11, 11] )
            idx = np.nanargmin( slp_smooth )
            lat_storm.append( ncid.variables['XLAT'][:].flatten()[idx]) 
            lon_storm.append( ncid.variables['XLONG'][:].flatten()[idx]) 

    max_wind = np.array(max_wind)
    min_slp = np.array(min_slp)
    lat_storm = np.array(lat_storm)
    lon_storm = np.array(lon_storm)

    HPI_analyses = {'Diff_start_next6x': 6,'time':DAtime_str, 'lat': lat_storm, 'lon':lon_storm, 'max_ws': max_wind, 'min_slp':min_slp}
    print('Location of analysis at ', HPI_analyses['time'][0], ': ', HPI_analyses['lon'][0], ' ', HPI_analyses['lat'][0])

    return HPI_analyses

# ------------------------------------------------------------------------------------------------------
#            Object: Deterministic Forecast; Operation: Read and Process  
# ------------------------------------------------------------------------------------------------------

# Get DF forecast's HPI from model output
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


# Read out a forecast's hourly HPI from ATCF part in rsl.error.0000
def read_rsl_error(Storm, Exper_name, directory, DF_start, DF_end):
    
    """
    This function reads out hourly HPI of the forecast (wrf files) within the time period of interest. Before you use this function,
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

    # Check if the file exists
    if os.path.exists( directory + '/ATCF_rsl.error.0000' ):
        with open(directory + '/ATCF_rsl.error.0000') as f:
            all_lines = f.readlines()
    else:
        return None        

    # Read and process records
    DF_start_real = []  # if the value of "time_to_move" is not set as 0 in WRF, DF_start_real is not equal to DF_start.
                        # E.g., Set "time_to_move = 60,60,60" and the initialization time (DF_start) is at '201709160000'.
                        #       the WRF program will not output HPI data to rsl.error.0000 until '201709160100'.
                        #       In this case, DF_start_real is '201709160100'.
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
        time_dt = datetime.strptime( time_wrf,"%Y-%m-%d_%H:%M:%S" ) #2017-09-16_03:00:00
        time_all.append( time_dt.strftime( "%Y%m%d%H%M" ) ) #201709160300
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

    # Only reads out hourly HPI records
    # - Calculate the distance from the real start point to the end of the forecast in hour
    DF_diff = datetime.strptime(DF_end,"%Y%m%d%H%M") - datetime.strptime(DF_start_real[0],"%Y%m%d%H%M")
    DF_diff_hour = DF_diff.total_seconds() / 3600
    # - List these times (in string format) every hour
    time_interest_dt = [datetime.strptime(DF_start_real[0],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(DF_diff_hour+1), 1))]
    time_interest_str = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Reads out a subset of records corresponding to times of interest
    idx_sbs = []
    lat_sbs = []
    lon_sbs = []
    max_ws_sbs = []
    min_slp_sbs = []

    for time_str in time_interest_str[:-1]:
        boolean_compare = [ eachT  == time_str for eachT in time_all ]
        if any(boolean_compare):
            if len( np.where(boolean_compare)[0] ) > 1:
                idx = int( np.where(boolean_compare)[0][0] )
            else:
                idx = int( np.where(boolean_compare)[0] )
            idx_sbs.append(idx)
            lat_sbs.append( lat_all[idx] )
            lon_sbs.append( lon_all[idx] )
            max_ws_sbs.append( maxV_all[idx] )
            min_slp_sbs.append( minP_all[idx] )

    # Add the last record (e.g.,2017-09-07_23:45:00) to approximate o' clock (e.g.,2017-09-08_00:00:00) 
    lat_sbs.append( lat_all[-1] )
    lon_sbs.append( lon_all[-1] )
    max_ws_sbs.append( maxV_all[-1] )
    min_slp_sbs.append( minP_all[-1] )  

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

    lat_sbs = np.array(lat_sbs)
    lon_sbs = np.array(lon_sbs)
    max_ws_sbs = np.array(max_ws_sbs)
    min_slp_sbs = np.array(min_slp_sbs)

    # Assemble the dictionary
    dict_model = {'Diff_start_next6x': int(Diff_start_next6x), 'time': time_interest_str, 'lat': lat_sbs, 'lon': lon_sbs, 'max_ws': max_ws_sbs, 'min_slp': min_slp_sbs} 
    return dict_model


# ------------------------------------------------------------------------------------------------------
#            Operation: Plot HPI values
# ------------------------------------------------------------------------------------------------------

def plot_one_hpi( ax0, ax1, ax2,  state, color, line, line_width, label, steps=6, btk=None):
 
    # Best track
    if label == 'Best track':
        # Read data
        times = state['time']
        lon = state['lon']
        lat = state['lat']
        x_min_slp = state['min_slp']
        x_max_ws = state['max_ws']

        # track
        ax0.plot(lon, lat, marker='o',markersize=3, color=color,linewidth=3, label=label, linestyle=line, transform=ccrs.PlateCarree())
        dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]
        for it in times:
            if it[8:10] == '00':
                boolean_compare = [ it  == eachT for eachT in times ]
                idx = int( np.where(boolean_compare)[0] )
                ax0.scatter( lon[idx], lat[idx],s=5, marker='o',edgecolor="white",transform=ccrs.PlateCarree())
                ax0.annotate(it[6:8], xy=(lon[idx], lat[idx]), color=color, xycoords='data', transform=ccrs.PlateCarree())
        # intensity
        ax1.plot_date(dates, x_min_slp, color, label=label, linewidth=3)
        ax2.plot_date(dates, x_max_ws, color, label=label, linewidth=3)

    # Model simulation
    else:
        # Read data (and process it if it is special)
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
            
            times = np.insert(times,0,state['time'][0])
            lon = np.insert(lon,0,state['lon'][0])
            lat = np.insert(lat,0,state['lat'][0])
            x_min_slp = np.insert(x_min_slp,0,state['min_slp'][0])
            x_max_ws = np.insert(x_max_ws,0,state['max_ws'][0])
        else:
            times = state['time'][::steps]
            lon = state['lon'][::steps]
            lat = state['lat'][::steps]
            x_min_slp = state['min_slp'][::steps]
            x_max_ws = state['max_ws'][::steps]

        # track        
        ax0.plot(lon, lat, marker='o', markersize=3, color=color,linewidth=line_width, label=label, linestyle=line, transform=ccrs.PlateCarree())
        for it in times:
            if it[8:10] == '00':
                boolean_compare = [ it  == eachT for eachT in times ]
                idx = int( np.where(boolean_compare)[0] )
                ax0.scatter( lon[idx], lat[idx],s=5, marker='o',edgecolor="white",transform=ccrs.PlateCarree())    
                ax0.annotate(it[6:8], xy=(lon[idx], lat[idx]), color=color, xycoords='data', transform=ccrs.PlateCarree())
                #ax0.text((state['lon'][idx], state['lat'][idx],it[6:8], fontsize=4,transform=ccrs.PlateCarree())

        # intensity
        dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]
        ax1.plot_date(dates, x_min_slp, color=color, label=label, linestyle=line, linewidth=line_width,)
        ax2.plot_date(dates, x_max_ws, color=color, label=label, linestyle=line, linewidth=line_width,)


def plot_hpi_df( Config ):

    wrf_dir = Config[0]
    Storm = Config[1]
    Exper = Config[2]
    mem_id = Config[4]
    read_fc_wrfout = Config[5]
    Plot_analyses = Config[6]

    # Set the end point of the period to investigate
    if Storm == 'HARVEY':
        DF_model_end  = '201708270000'
    elif Storm == 'IRMA':
        DF_model_end  = '201709080000'
    elif Storm == 'MARIA':
        DF_model_end  = '201709210000'
    elif Storm == 'JOSE':
        DF_model_end  = '201709100000'
    else:
        DF_model_end  = None

    # Identify the available forecasts with either their initialization times or their id members
    Exper_content_lbl = {}
    for iExper in Exper:
        if iExper is not None:
            #Exper_content_lbl[iExper] = sorted(fnmatch.filter(os.listdir( wrf_dir+'/'+Storm+'/'+iExper+'/' ),'20*input')) 
            if os.path.exists( wrf_dir+'/'+Storm+'/'+iExper+'/wrf_df/' ):
                Exper_content_lbl[iExper] = sorted(fnmatch.filter(os.listdir( wrf_dir+'/'+Storm+'/'+iExper+'/wrf_df/' ),'20*'))  
            else:
                Exper_content_lbl[iExper] = None
        else:
             Exper_content_lbl[iExper] = None

    # Set up figure
    fig = plt.figure( figsize=(22,6.5), dpi=150 ) #12,4 20,6
    gs = fig.add_gridspec(2,7) # 1,3
    ax0 = fig.add_subplot( gs[0:,0:3],  projection=ccrs.PlateCarree())
    ax0.set_extent( domain_range,  crs=ccrs.PlateCarree())
    ax0.coastlines( resolution='10m', color='black',linewidth=0.5 )
    ax1 = fig.add_subplot( gs[:,3:5] )
    ax2 = fig.add_subplot( gs[:,5:7] )

    # Set start and end of the period
    if Storm == 'HARVEY':
        Btk_start = '201708221200' # '201709161800' #'201709030600'
        Btk_end = '201708270000' # '201709210000' #'201709090000'
    elif Storm == 'IRMA':
        Btk_start = '201709030000'
        Btk_end = '201709080000' #201709080000
    elif Storm == 'MARIA':
        Btk_start = '2017091600'#'201709160000'
        Btk_end = '201709210000'
    elif Storm == 'JOSE':
        Btk_start = '201709050000'
        Btk_end = '201709100000'
    else:
        pass
    
    # Plot HPI from post-storm analysis
    best_track = UD.btk_in_duration(Storm, Btk_start, Btk_end, hour_step=6)
    plot_one_hpi( ax0, ax1, ax2, best_track,  'black', '-', 5, 'Best track')

    # Customize color maps
    if distinct_colors:
        Color1 = ["#c23728","#e14b31","#de6e56","#e1a692","#786028","#a57c1b","#d2980d","#ffb400","#503f3f","#6d4b4b","#a86464","#e27c7c"] #redish
        Color2 = ["#115f9a", "#1984c5", "#22a7f0", "#48b5c4", "#48446e", "#5e569b", "#776bcd", "#9080ff","#3c4e4b", "#466964", "#599e94", "#6cd4c5"] #blueish
    else:
        red_cm = cm.Reds
        blue_cm = cm.Blues
        num_colors = 14
        discretize_red = ListedColormap(red_cm(np.linspace(0,1,num_colors)))
        discretize_blue = ListedColormap(blue_cm(np.linspace(0,1,num_colors)))
        Color1 = discretize_red.colors[5:]
        Color2 = discretize_blue.colors[5:]
    Color_set = {'c0':Color1, 'c1':Color2}
    Ana_color = ['#eea990','#748b97'] #'#748b97'
    Line_types = ['-','-']
    
    # Customize labels ###### Chnage it every time !!!!!!!!!!!!!!! 
    #Labels = ['GTS+HPI(hroi:300km): ']
    Labels = ['WSM6:','THO:']
    #Labels = ['Stp2-Intel17 ','Eps-Intel19 ']
    Ana_labels = ['Ref Ans','MD Ans' ]
    #Ana_labels = ['Stampede2 Analysis', 'Expanse Analysis']

    # Plot HPI for each deterministic forecasts
    if read_fc_wrfout == True:
        i = 0
        for init_time in DF_init_times:
            #plot_one_hpi( ax0, ax1, ax2, read_wrfout(Storm, wrf_dir+'/'+Storm+'/newWRF_MW_THO/wrf_df/'+init_time, force_reload=reload_data), color_model[i], 'IR+MW: '+ init_time, step=step)
            pass
            i = i+1

    elif read_fc_wrfout == False:
        iExper = 0
        for key in Exper_content_lbl: 
            print('Experiment: ', key)
            if Exper_content_lbl[key] is not None:
                ic = 0
                if Plot_analyses == True:
                    print('Plotting the analyses...')
                    HPI_analyses = read_HPI_analyses(Storm, key, wrf_dir)
                    plot_one_hpi( ax0, ax1, ax2, HPI_analyses, Ana_color[iExper], '--', 2.5, Ana_labels[iExper], steps=1 )
                for it in Exper_content_lbl[key]:
                    print('Plotting ', it)
                    #HPI_model = read_rsl_error(Storm, key, wrf_dir+'/'+Storm+'/'+key+'/'+it, '201709030000', DF_model_end)
                    #print(wrf_dir+'/'+Storm+'/'+key+'/wrf_df/'+it)
                    HPI_model = read_rsl_error(Storm, key, wrf_dir+'/'+Storm+'/'+key+'/wrf_df/'+it, it, DF_model_end)
                    plot_one_hpi( ax0, ax1, ax2, HPI_model, Color_set['c'+str(iExper)][ic], Line_types[iExper], 2.5, Labels[iExper]+it, steps=6 )
                    #if 'THO' in it:
                    #    plot_one_hpi( ax0, ax1, ax2, HPI_model,'#B1D866', '-', 5, 'conv_THO', steps=1 )
                    #else:
                    #    plot_one_hpi( ax0, ax1, ax2, HPI_model,'#FFED13', '-', 5, 'conv_WSM6', steps=1 )
                    ic = ic + 1
            else:
                print('No available data!')
            iExper = iExper + 1

    # Set ticks/labels for track subplot
    lon_ticks = list(range(math.ceil(domain_range[0])-2, math.ceil(domain_range[1])+3, 4))
    lat_ticks = list(range(math.ceil(domain_range[2])-2, math.ceil(domain_range[3])+2, 2))
    gl = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.left_labels = True
    gl.bottom_labels = True
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlabel_style = {'size': 15}
    gl.ylabel_style = {'size': 15}  
   
    # Set ticks/other attributes for intensity subplots
    ax1.set_xlim([datetime(int(Btk_start[0:4]), int(Btk_start[4:6]), int(Btk_start[6:8]), int(Btk_start[8:10])), datetime(int(Btk_end[0:4]), int(Btk_end[4:6]), int(Btk_end[6:8]), int(Btk_end[8:10]))])
    ax2.set_xlim([datetime(int(Btk_start[0:4]), int(Btk_start[4:6]), int(Btk_start[6:8]), int(Btk_start[8:10])), datetime(int(Btk_end[0:4]), int(Btk_end[4:6]), int(Btk_end[6:8]), int(Btk_end[8:10]))])
    ax1.set_ylim([900,1020])     #([940, 1015])
    ax2.set_ylim([10,80])   #([10,60])
    ax1.tick_params(axis='x', labelrotation=30,labelsize=12)
    ax2.tick_params(axis='x', labelrotation=30,labelsize=12)
    ax2.legend(bbox_to_anchor=(1.5, 1.0),frameon=True,loc='upper right',fontsize='10') #10
    ax1.grid(True,linewidth=1, color='gray', alpha=0.5, linestyle='-')
    ax2.grid(True,linewidth=1, color='gray', alpha=0.5, linestyle='-')

    # Set titles
    ax0.set_title( 'Track',fontsize = 15 )
    ax1.set_title( 'MSLP (hPa)',fontsize = 15 )
    ax2.set_title( 'Vmax ($\mathregular{ms^{-1}}$)',fontsize = 15 )
    fig.suptitle('conv',fontsize = 15)

    # Save figure
    des_name = small_dir+Storm+'/'+Exper[1]+'/Vis_analyze/Model/'+Storm+'_forecast_conv_WSM6_THO.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


# ------------------------------------------------------------------------------------------------------
#            Operation: Plot absolute errors and their mean
# ------------------------------------------------------------------------------------------------------
# Sample number are the same across different forecasts
def SameL_error_eachInit(Storm,wrf_dir,exper,best_track,exper_inits,DF_model_end,run_hrs):
    
    dict_Error = {}
    run_times = list(range(0, run_hrs+1, 6))
    # Loop over initialization time for each experiment
    for init_t in exper_inits:
        error_perInit = np.zeros((3,len(run_times))) # track,min slp,max wind
        error_perInit[:] = np.nan
        # list forecast times every 6 hours
        fc_times_dt = [datetime.strptime(init_t,"%Y%m%d%H%M") + timedelta(hours=t) for t in run_times]
        fc_times_str = [time_dt.strftime("%Y%m%d%H%M") for time_dt in fc_times_dt] 
        # read from wrf forecasts
        if not os.path.exists( wrf_dir+'/'+Storm+'/'+exper+'/wrf_df/'+init_t ):
            continue
        # calculate the error!!!
        Diff_start_next6x = HPI_model['Diff_start_next6x']
        for it in fc_times_str:
            idx_t = fc_times_str.index(it)
            # special treatment to t0
            bool_t0 = Diff_start_next6x < 6
            if idx_t == 0 and bool_t0:
                continue # value will be np.nan
            # index of the data point in best_track
            bool_btk = [ it == btk_str for btk_str in best_track['time'][:] ]
            idx_btk = int( np.where(bool_btk)[0] ) 
            # index of the data point in model
            bool_model = [ it == model_str for model_str in HPI_model['time'][:] ]
            #if not any(bool_model): # speical treatment to data at the end time
            #    continue
            #else:
            idx_m = int( np.where(bool_model)[0] )
            # error of track
            distance_track = UD.mercator_distance(best_track['lon'][idx_btk],best_track['lat'][idx_btk],HPI_model['lon'][idx_m],HPI_model['lat'][idx_m]) # km
            error_perInit[0,idx_t] = distance_track
            # error of min slp
            error_perInit[1,idx_t] = HPI_model['min_slp'][idx_m]-best_track['min_slp'][idx_btk] 
            # error of max wind
            error_perInit[2,idx_t] = HPI_model['max_ws'][idx_m]-best_track['max_ws'][idx_btk]
        # assemble the dictionary
        dict_Error[init_t] = error_perInit
        # accumulate the abs error
        if exper_inits.index(init_t) == 0:
            abs_error_track = error_perInit[0,:].reshape(len(run_times),1)
            abs_error_slp = abs(error_perInit[1,:]).reshape(len(run_times),1)
            abs_error_maxws = abs(error_perInit[2,:]).reshape(len(run_times),1)
        else:
            abs_error_track = np.concatenate((abs_error_track, abs(error_perInit[0,:]).reshape(len(run_times),1)), axis=1)
            abs_error_slp = np.concatenate((abs_error_slp,abs(error_perInit[1,:]).reshape(len(run_times),1)), axis=1)
            abs_error_maxws = np.concatenate((abs_error_maxws,abs(error_perInit[2,:]).reshape(len(run_times),1)), axis=1)

    # calculate the mean
    assert abs_error_track.shape[1] == len(exper_inits)  
    abs_error_mean = np.zeros((3,len(run_times)))
    abs_error_mean[0,:] = np.nanmean(abs_error_track,axis=1)
    abs_error_mean[1,:] = np.nanmean(abs_error_slp,axis=1)
    abs_error_mean[2,:] = np.nanmean(abs_error_maxws,axis=1)
    dict_Error['mean_abs_error'] = abs_error_mean

    return dict_Error 

# Sample number are different across different forecasts 
# I.e., sample numbers are decided by the start and end of the forecast time
def DiffL_error_eachInit(Storm,wrf_dir,exper,best_track,exper_inits,DF_model_end):

    dict_Error = {}
    # Loop over initialization time for each experiment
    for init_t in exper_inits:
		# generate the forecast time series every 6 hours
        fc_times_str = UD.generate_times(init_t,DF_model_end,6) 
        run_times = len(fc_times_str)
		# initialize error array
        error_perInit = np.zeros((4,run_times)) # position_distance,position_azimuth relative to the best track, min slp,max wind
        error_perInit[:] = np.nan
        # read from wrf forecasts
        if not os.path.exists( wrf_dir+'/'+Storm+'/'+exper+'/wrf_df/'+init_t ):
            continue
        HPI_model = read_rsl_error(Storm,exper,wrf_dir+'/'+Storm+'/'+exper+'/wrf_df/'+init_t, init_t, DF_model_end)
        # calculate the error!!!
        Diff_start_next6x = HPI_model['Diff_start_next6x']

        for it in fc_times_str:
            idx_t = fc_times_str.index(it)
            # special treatment to t0
            bool_t0 = Diff_start_next6x < 6
            if idx_t == 0 and bool_t0:
                continue # value will be np.nan
            # index of the data point in best_track
            bool_btk = [ it == btk_str for btk_str in best_track['time'][:] ]
            idx_btk = int( np.where(bool_btk)[0] )
            # index of the data point in model
            bool_model = [ it == model_str for model_str in HPI_model['time'][:] ]
            #if not any(bool_model): # speical treatment to data at the end time
            #    continue
            #else:
            idx_m = int( np.where(bool_model)[0] )
            # error of track
            dist, azi_btk = UD.twoP_inverse(best_track['lat'][idx_btk],best_track['lon'][idx_btk],HPI_model['lat'][idx_m],HPI_model['lon'][idx_m]) # km
            error_perInit[0,idx_t] = dist #km
            error_perInit[1,idx_t] = UD.azi_geo_to_math(azi_btk)  #math azimuth:degree
            # error of min slp
            error_perInit[2,idx_t] = HPI_model['min_slp'][idx_m]-best_track['min_slp'][idx_btk]
            # error of max wind
            error_perInit[3,idx_t] = HPI_model['max_ws'][idx_m]-best_track['max_ws'][idx_btk]
        # assemble the dictionary
        dict_Error[init_t] = error_perInit 

    return dict_Error

def plot_abs_err( Config ):

    wrf_dir = Config[0]
    Storm = Config[1]
    Exper = Config[2]
    mem_id = Config[4]
    read_fc_wrfout = Config[5]
    Plot_analyses = Config[6]

    # Set the end point of the period to investigate
    if Storm == 'HARVEY':
        DF_model_end  = '201708270000'
    elif Storm == 'IRMA':
        DF_model_end  = '201709080000'
    elif Storm == 'MARIA':
        DF_model_end  = '201709210000'
    elif Storm == 'JOSE':
        DF_model_end  = '201709100000'
    else:
        DF_model_end  = None

    # Identify the available forecasts with either their initialization times or their id members
    Exper_content_lbl = {}
    for iExper in Exper:
        if iExper is not None:
            if os.path.exists( wrf_dir+'/'+Storm+'/'+iExper+'/wrf_df/' ):
                Exper_content_lbl[iExper] = sorted(fnmatch.filter(os.listdir( wrf_dir+'/'+Storm+'/'+iExper+'/wrf_df/' ),'20*'))
            else:
                Exper_content_lbl[iExper] = None
        else:
             Exper_content_lbl[iExper] = None
   
    # Read best-track data
    if Storm == 'HARVEY':
        Btk_start = '201708221200' # '201709161800' #'201709030600'
        Btk_end = '201708270000' # '201709210000' #'201709090000'
    elif Storm == 'IRMA':
        Btk_start = '201709030000'
        Btk_end = '201709080000' #201709080000
    elif Storm == 'MARIA':
        Btk_start = '2017091600'#'201709160000'
        Btk_end = '201709210000'
    elif Storm == 'JOSE':
        Btk_start = '201709050000'
        Btk_end = '201709100000'
    else:
        pass
    best_track = UD.btk_in_duration(Storm, Btk_start, Btk_end, hour_step=6)
 
    # Calculate the error with lead times 
    Exper_absError_lbl = {}
    for iExper in Exper:
        print(iExper)
        if Exper_content_lbl[iExper] is not None:
            Exper_absError_lbl[iExper] = SameL_error_eachInit(Storm,wrf_dir,iExper,best_track,Exper_content_lbl[iExper],DF_model_end,run_hrs=fc_run_hrs)
        else:
             Exper_absError_lbl[iExper] = None 

    # ------------- Plot ---------------------------
    fig, ax=plt.subplots(1,3,sharex='all', figsize=(15,6.5), dpi=400)
   
    # Customize color maps
    if distinct_colors:
        Color1 = ["#c23728","#e14b31","#de6e56","#e1a692","#786028","#a57c1b","#d2980d","#ffb400","#503f3f","#6d4b4b","#a86464","#e27c7c"] #redish
        Color2 = ["#115f9a", "#1984c5", "#22a7f0", "#48b5c4", "#48446e", "#5e569b", "#776bcd", "#9080ff","#3c4e4b", "#466964", "#599e94", "#6cd4c5"] #blueish
    else:
        red_cm = cm.Reds
        blue_cm = cm.Blues
        green_cm = cm.Greens
        num_colors = 14
        discretize_red = ListedColormap(red_cm(np.linspace(0,1,num_colors)))
        discretize_blue = ListedColormap(blue_cm(np.linspace(0,1,num_colors)))
        discretize_green = ListedColormap(green_cm(np.linspace(0,1,num_colors)))
        Color1 = discretize_red.colors[5:]
        Color2 = discretize_blue.colors[5:]
        Color3 = discretize_green.colors[5:]
    Color_set = {'c0':Color1, 'c1':Color2,'c2':Color3}
    mean_color = ['#D70040','#0F52BA','#26580F'] 
    Line_types = ['-','-',]
    mean_labels = ['MD_mem','MD_mean']
    #mean_labels = ['IC:Ref;FC:Ref','IC:MD;FC:MD','IC:MD;FC:Ref']
    # X axis: leading times
    lead_t = list(range(0, fc_run_hrs+1, 6))

    for exper in Exper:
        iE = Exper.index( exper )
        err = Exper_absError_lbl[exper]
        # loop thru each forecast
        ic = 0
        for init_t in Exper_content_lbl[exper]:   
            ax[0].plot(lead_t,err[init_t][0,:],color=Color_set['c'+str(iE)][ic], linestyle=Line_types[iE],alpha=0.6)
            ax[1].plot(lead_t,abs(err[init_t][1,:]),color=Color_set['c'+str(iE)][ic], linestyle=Line_types[iE],alpha=0.6)
            ax[2].plot(lead_t,abs(err[init_t][2,:]),color=Color_set['c'+str(iE)][ic], linestyle=Line_types[iE],alpha=0.6)
            ic = ic + 1
        # plot the mean
        ax[0].plot(lead_t,err['mean_abs_error'][0,:],color=mean_color[iE],linestyle='-',linewidth=6)
        ax[1].plot(lead_t,err['mean_abs_error'][1,:],color=mean_color[iE],linestyle='-',linewidth=6)
        ax[2].plot(lead_t,err['mean_abs_error'][2,:],color=mean_color[iE],linestyle='-',linewidth=6,label=mean_labels[iE])
    
    # Details
    ax[2].legend(frameon=True,loc='upper right',fontsize='17')
    ax[0].set_ylim([0,200])
    ax[1].set_ylim([0,50])
    ax[2].set_ylim([0,50])
    for i in range(3):
        ax[i].set_xticks( lead_t[::2] )
        ax[i].set_xlim([0,fc_run_hrs])
        ax[i].tick_params(axis='x',labelsize=20)
        ax[i].grid(True,linewidth=1, color='gray', alpha=0.5, linestyle='-')

    # Set titles and axis label
    fig.text(0.52, 0.02, 'Forecast Lead Time (hour)', ha='center',fontsize=20)
    ax[0].set_ylabel('Track Error (km)',fontsize=15)
    ax[1].set_ylabel( 'MSLP Error (hPa)',fontsize=15)
    ax[2].set_ylabel( 'Vmax Error ($\mathregular{ms^{-1}}$)',fontsize=15)
    ax[0].set_title( 'Track (km)',fontsize = 20 )
    ax[1].set_title( 'MSLP (hPa)',fontsize = 20 )
    ax[2].set_title( 'Vmax ($\mathregular{ms^{-1}}$)',fontsize = 20 )
    fig.suptitle(Storm+'~'+DA+': Mean Absolute Error',fontsize = 20)

    # Save figure
    des_name = small_dir+Storm+'/'+Exper[1]+'/Vis_analyze/Model/'+Storm+'_forecast_abs_error_IR_TuneWSM6_mem_mean.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name+'!' )


if __name__ == '__main__':
    
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # Configuration
    Storm = 'IRMA'
    MP = ['WSM6','THO']
    DA = 'conv'
    DF_model_start = '20170822180000' # Default value of DF_model_start. Especially useful when dealing with ensemble forecast
    mem_id = 'mean' # Default value of member id. Especially useful when dealing with deterministic forecast
    read_fc_wrfout = False # Feature that decides the way of reading HPI from model files

    distinct_colors = False
    Plot_analyses = False # Feature that plots the analyses of an experiment
    Plot_HPI = True
    Plot_abs_error = False
    fc_run_hrs = 60
    
    # Generate experiment names
    Exper_name = []
    for imp in MP:
        Exper_name.append( UD.generate_one_name( Storm,DA,imp ) )

    wrf_dir = big_dir
    Config = [wrf_dir, Storm, Exper_name, DF_model_start, mem_id, read_fc_wrfout, Plot_analyses]

    # Pre-set the domain for DF forecast
    if Storm == 'HARVEY':
        lon_min = -102
        lon_max = -85
        lat_min = 16
        lat_max = 31
    elif Storm == 'IRMA':
        lon_min = -72
        lon_max = -45
        lat_min = 10
        lat_max = 30
    elif Storm == 'MARIA':
        lon_min = -68
        lon_max = -43#-45
        lat_min = 5# 
        lat_max = 25
    elif Storm == 'JOSE':
        lon_min = -65
        lon_max = -35#-45
        lat_min = 0
        lat_max = 25
    domain_range = [lon_min, lon_max, lat_min, lat_max]

    if Plot_HPI:
        plot_hpi_df( Config )

    if Plot_abs_error:
        start_time=time.process_time()
        plot_abs_err( Config ) 
        end_time = time.process_time()
        print('time needed: ', end_time-start_time, ' seconds')


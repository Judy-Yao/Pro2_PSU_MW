#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

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
from math import modf
import matlab.engine
import scipy as sp
import scipy.ndimage
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from Track_intensity_all import read_bestrack
from mpl_toolkits.axes_grid1 import make_axes_locatable


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

def plot_UV10_slp( wrfout,plot_dir, dict_AF_masked, dict_btk, if_btk_exist, idx_btk):

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
    # Mark the best track
    if if_btk_exist:
        ax.scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 8, 'green', marker='*',transform=ccrs.PlateCarree()) 
    #ax.scatter(lon_minslp, lat_minslp, s=5, marker='*', edgecolors='red', transform=ccrs.PlateCarree())
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
    ax.set_title(wrfout_head_tail[1].replace('wrfout_d03_',' '),  fontweight='bold', fontsize=10)

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
def plot_d03(Storm, wrf_dir, plot_dir, USAF_per_day):
    
    # Best-track
    dict_btk = read_bestrack(Storm)

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
            wrfout_head_tail_all =  os.path.split( wrfout ) 
            nowt_str = datetime.strptime(wrfout_head_tail_all[1][11:27],"%Y-%m-%d_%H:%M")
            # Find the best-track position
            btk_dt = [datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]
            bool_match = [nowt_str == it for it in btk_dt]
            if True in bool_match:
                if_btk_exist = True
                idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
            else:
                if_btk_exist = False
            
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
            plot_UV10_slp( wrfout, plot_dir, dict_AF_masked, dict_btk, if_btk_exist, idx_btk )



def UV10_slp_AF( Storm, Exper_name, big_dir, small_dir ):

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
            plot_d03( Storm, wrf_dir, plot_dir, USAF_per_day )
        else:
            plot_d03( Storm, wrf_dir, plot_dir, USAF_per_day )

# Interpolate model points to one obs location
def ijk_obs(lon_x, lat_x, H_x, lon_f, lat_f, H_f ):

    nx = 297
    ny =297
    nz = 42
    
    print('Geolocation of the obs is: ', lon_f, lat_f, H_f)
    # find the four nearest horizontal model points around the obs location
    if (lon_x[0,0] > lon_f) or (lon_x[0,-1] < lon_f) or (lat_x[0,0] > lat_f) or (lat_x[-1,0] < lat_f):
        return None
        print('The obs is outside the domain!')

    for i in range( nx ):
        if lon_x[0,i] > lon_f:
            if i == 0:
                print('The i+1 is on the left border!')
                #i_PlusOne = i
                #i_MinusOne = i 
                #print('Find the i and i+1!', i_MinusOne, i_PlusOne )
                return None
            else:
                i_PlusOne = i
                i_MinusOne = i - 1
                print('Find the i and i+1!', i_MinusOne, i_PlusOne )
                break

    for j in range( ny ):
        if lat_x[j,0] > lat_f:
            if j == 0:
                print('The i+1 is on the bottom border!')
                #j_PlusOne = j
                #j_MinusOne = j
                #print('Find the j and j+1!', j_MinusOne, j_PlusOne )
                return None
            else:
                j_PlusOne = j
                j_MinusOne = j  - 1
                print('Find the j and j+1!', j_MinusOne, j_PlusOne )
                break

    # convert the location of the obs in the geo coordinate to the ijk coordinate in the horizontal plane
    i_prime = (lon_f - lon_x[0,i_MinusOne])/(lon_x[0,i_PlusOne]-lon_x[0,i_MinusOne])*(i_PlusOne-i_MinusOne)+i_MinusOne
    j_prime = (lat_f - lat_x[j_MinusOne,0])/(lat_x[j_PlusOne,0]-lat_x[j_MinusOne,0])*(j_PlusOne-j_MinusOne)+j_MinusOne
    #print('The location of the obs in the ith direction is: ', i_prime)
    #print('The location of the obs in the jth direction is: ', j_prime)

    # --- interpolate the height profile of model points to the obs location ---
    H_left_bot = H_x[:,j_PlusOne,i_MinusOne]
    H_left_up = H_x[:,j_PlusOne,i_MinusOne]
    H_right_bot = H_x[:,j_PlusOne,i_PlusOne]
    H_right_up = H_x[:,j_PlusOne,i_PlusOne]
    
    # perform the bilinear interpolation over an unit area
    area = (i_PlusOne-i_MinusOne)*(j_PlusOne-j_MinusOne)
    if area != 1:
        raise Exception('The area is not a unit square!') # raise an error and stop the program
    
    fra_i_prime = modf(i_prime)[1]
    fra_j_prime = modf(j_prime)[1]
    H_ij_prime = H_left_bot*(1-fra_i_prime)*(1-fra_j_prime) + H_right_bot*fra_i_prime*(1-fra_j_prime) + H_left_up*fra_j_prime*(1.0-fra_i_prime) + H_right_up*fra_i_prime*fra_j_prime
    
    # find the two nearest vertical model points around the obs location
    if (H_ij_prime[0] > H_f) or (H_ij_prime[-1] < H_f):
        return None
        print('The obs is outside the domain!')
    for z in range( nz ):
        if H_ij_prime[z] > H_f:
            if z == 0:
                print('The z+1 is on the bottom border!')
                #z_PlusOne = z
                #z_MinusOne = z
                #print('Find the z and z+1!', z_MinusOne, z_PlusOne )
                return None
            else:
                z_PlusOne = z
                z_MinusOne = z_PlusOne - 1
                print('Find the z and z+1!', z_MinusOne, z_PlusOne )
                break

    # convert the location of the obs in the geo coordinate to the ijk coordinate in the vertical plane
    z_prime = (H_f - H_ij_prime[z_MinusOne])/(H_ij_prime[z_PlusOne]-H_ij_prime[z_MinusOne])*(z_PlusOne-z_MinusOne)+z_MinusOne

    print('The location of the obs in the ijk coordinate is: ', i_prime, j_prime, z_prime)
    return [i_prime, j_prime, z_prime]

# Plot hourly flight track + variable values of model output at XX UTC to the flight level points (lon,lat,height) within XX UTC +/- 10 minutes
def Plot_compare_track_wind_td( Storm, small_dir, df_dir, plot_dir ):

    # List the USFA files to investigate
    USAF_list = sorted(glob.glob( small_dir + Storm + '/USAF/201708*' ))
    # Define attributes of interest to read
    attrs_itt = ['GMT','GPSA','LAT','LON','WSpd','TD'] # GMT time / GPS Altimeter (height of the air plane)/ latitude / longitude/ wind speed / dew point temperature
    
    # List all of wrf files and extract their times
    wrfout_list = sorted(glob.glob( df_dir+'/wrfout_d03_2017*' ))  
    wrfout_all_time_str = []
    for wrfout in wrfout_list:
        wrfout_head_tail_all =  os.path.split( wrfout ) 
        wrfout_all_time_str.append( datetime.strftime(datetime.strptime(wrfout_head_tail_all[1][11:27],"%Y-%m-%d_%H:%M"),"%Y%m%d%H") ) #e.g., '2017091810'
    
    eng = matlab.engine.start_matlab() # start a new matlab process
    # Loop through each mission
    for imission in USAF_list:
        print( '--------- The mission is ', imission, '--------------' )
        # Read the obs from the mission
        mission_head_tail = os.path.split( imission )
        time_mission = mission_head_tail[1][0:8] # e.g., '20170918'
        dict_AF_mission,uni_hhdd = USAF.read_USAF_mission( time_mission, imission, attrs_itt )

        # Best-track
        dict_btk = read_bestrack(Storm)
        btk_dt = [datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]

        # Loop through each unique day and hour
        for it_str in uni_hhdd:
            print('Dealing with time: ', it_str) 
            bool_idx = [it == it_str for it in wrfout_all_time_str ]
            if np.any( bool_idx ) == True:
                idx_wrf = int(np.where(bool_idx)[0]) # identify the location of the wrfout
            else:
                print('No wrfout is available at this time!')
                continue
           
            # Get flight-level track 
            it_dt = datetime.strptime(it_str,"%Y%m%d%H")
            Nminutes = 10
            dict_AF_hour = USAF.mask_time_of_t( dict_AF_mission, it_dt, Nminutes )
            if dict_AF_hour == None:
                print('No AF obs is available at this time span!')
                continue
           
            # Find the best-track position
            bool_match = [it_dt == it for it in btk_dt]
            if True in bool_match:
                if_btk_exist = True
                idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
            else:
                if_btk_exist = False

            # Read the wrf variables
            print('Check model data......')
            ncdir = nc.Dataset( wrfout_list[idx_wrf] )
            lat_x = ncdir.variables['XLAT'][0,:,:]#[::10,::10]
            lon_x = ncdir.variables['XLONG'][0,:,:]#[::10,::10]

            u = ncdir.variables['U'][0,:,:,:]
            u_mass = (u[:,:,0:-1:1] + u[:,:,1::1])/2 # interpolate staggered point to the mass grid 
            v = ncdir.variables['V'][0,:,:,:]
            v_mass = (v[:,0:-1:1,:] + v[:,1::1,:])/2 # interpolate staggered point to the mass grid 
            ws_x =  ((u_mass ** 2 + v_mass ** 2) ** 0.5)#[::4,::10,::10]

            PHB = ncdir.variables['PHB'][0,:,:,:]
            PH = ncdir.variables['PH'][0,:,:,:]
            GP = PHB+PH
            GP_mass = (GP[0:-1:1,:,:] + GP[1::1,:,:])/2 # interpolate staggered point to the mass grid 
            Height_mass_x = (GP_mass/9.8)#[::4,::10,::10] # Convert the geopotential to height
            
            td_x = getvar(ncdir, 'td', units='degC').values #[::4,::10,::10].values
           
            print('min_lon_x: ', np.amin(lon_x),'max_lon_x: ', np.amax(lon_x))
            print('min_lat_x: ', np.amin(lat_x),'max_lat_x: ', np.amax(lat_x)) 
            print('min_ws_x: ', np.amin(ws_x),'max_ws_x: ', np.amax(ws_x))
            print('min_height_x: ', np.amin(Height_mass_x),'max_height_x: ', np.amax(Height_mass_x))
            print('min_td_x: ', np.amin(td_x),'max_td_x: ', np.amax(td_x))
 
            # Interpolate the model values in the space that is centered around the flight-level obs
            lon_f = dict_AF_hour['LON']
            lat_f = dict_AF_hour['LAT']
            gpsa_f = dict_AF_hour['GPSA']

            print('min_lon_f: ', np.amin(lon_f),'max_lon_f: ', np.amax(lon_f))
            print('min_lat_f: ', np.amin(lat_f),'max_lat_f: ', np.amax(lat_f))
            print('min_height_f: ', np.amin(gpsa_f),'max_height_f: ', np.amax(gpsa_f))

            # Roughly check if obs is within the model domain
            if (np.amax(lon_f) < np.amin(lon_x)) or (np.amin(lon_f) > np.amax(lon_x)):
                print('The flight obs is not within the model domain!')
                continue
            if (np.amax(lat_f) < np.amin(lat_x)) or (np.amin(lat_f) > np.amax(lat_x)):
                print('The flight obs is not within the model domain!')
                continue

            idx_itp = []
            ws_AFspace = []
            td_AFspace = []
            for iobs in range( len(dict_AF_hour['LON']) ):
                print('------------------------------------------------')
                print('GMT: ', dict_AF_hour['GMT'][iobs])
                # get the ijk of obs
                list_ijk_obs = ijk_obs( lon_x, lat_x, Height_mass_x, lon_f[iobs], lat_f[iobs], gpsa_f[iobs] )
                if list_ijk_obs is None:
                    continue           

                # collect the ijk coordinates of the 8 points  
                i_PlusOne = np.ceil( list_ijk_obs[0] )
                i_MinusOne = np.floor( list_ijk_obs[0] )
                j_PlusOne = np.ceil( list_ijk_obs[1] )
                j_MinusOne = np.floor( list_ijk_obs[1] )
                k_PlusOne = np.ceil( list_ijk_obs[2] )
                k_MinusOne = np.floor( list_ijk_obs[2] )
                i_8_points = [i_MinusOne,i_PlusOne,i_PlusOne,i_MinusOne,i_MinusOne,i_PlusOne,i_PlusOne,i_MinusOne] 
                j_8_points = [j_MinusOne,j_MinusOne,j_PlusOne,j_PlusOne,j_MinusOne,j_MinusOne,j_PlusOne,j_PlusOne]
                k_8_points = [k_MinusOne,k_MinusOne,k_MinusOne,k_MinusOne,k_PlusOne,k_PlusOne,k_PlusOne,k_PlusOne]
                # griddata interpolate
                ws_8_points = []
                td_8_points = []
                for it in range(len(i_8_points)):
                    ws_8_points.append( ws_x[int(k_8_points[it]), int(j_8_points[it]), int(i_8_points[it])] )
                    td_8_points.append( td_x[int(k_8_points[it]), int(j_8_points[it]), int(i_8_points[it])] )
    
                ws_ith_obs = eng.griddata( matlab.double(i_8_points), matlab.double(j_8_points),  matlab.double(k_8_points), matlab.double( ws_8_points), matlab.double( list_ijk_obs[0:1] ), matlab.double( list_ijk_obs[1:2] ), matlab.double( list_ijk_obs[2:3] ) )
                print('interpolated wind:', np.array(ws_ith_obs) )
                print('mean wind:', np.mean( ws_8_points ) )
                td_ith_obs = eng.griddata( matlab.double(i_8_points), matlab.double(j_8_points),  matlab.double(k_8_points), matlab.double( td_8_points), matlab.double( list_ijk_obs[0:1] ), matlab.double( list_ijk_obs[1:2] ), matlab.double( list_ijk_obs[2:3] ) )
                print('interpolated td:', np.array(td_ith_obs) )
                print('mean td:', np.mean( td_8_points ) )

                if np.array(ws_ith_obs).size > 0:
                    idx_itp.append( iobs )
                    ws_AFspace.append( np.array(ws_ith_obs).tolist() )
                    td_AFspace.append( np.array(td_ith_obs).tolist() )
                else:
                    continue
           
            print('Check interpolated model date......')
            print('min_ws: ', np.amin(ws_AFspace[:]),'max_ws: ', np.amax(ws_AFspace[:])) 
            #print('min_td: ', np.amin(td_AFspace[:]),'max_td: ', np.amax(td_AFspace[:])) 

            #  --------- Plot --------------------
            fig = plt.figure( figsize=(8,9.5), dpi=150 )
            gs = fig.add_gridspec(5,2)
            
            #  --------- Plot the flight level track  --------------------
            ax0 = fig.add_subplot( gs[0:3,:],  projection=ccrs.PlateCarree())
            
            lon_min = -100#-71#np.amin( dict_AF_mission['lon'] )
            lon_max = -85#-57#np.amax( dict_AF_mission['lon'] )
            lat_min = 15#10#np.amin( dict_AF_mission['lat'] )
            lat_max = 31#20#np.amax( dict_AF_mission['lat'] )
            gpsa_min = 0
            gpsa_max = np.amax( dict_AF_mission['GPSA'] )

            ax0.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
            ax0.coastlines (resolution='10m', color='black', linewidth=1)
            #  Plot the track for the whole mission
            gpsa_min = 0
            gpsa_max = 10000
            ax0.scatter(dict_AF_mission['LON'], dict_AF_mission['LAT'], 1, 'grey', cmap='jet', vmin=gpsa_min, vmax=gpsa_max,transform=ccrs.PlateCarree())
            AF = ax0.scatter(dict_AF_hour['LON'], dict_AF_hour['LAT'], 2, dict_AF_hour['GPSA'], cmap='jet', vmin=gpsa_min, vmax=gpsa_max, transform=ccrs.PlateCarree())
            AF_bar = fig.colorbar(AF,ax=ax0,shrink=0.9)
            AF_bar.ax.set_ylabel('Flight Height (m)')
            AF_bar.ax.tick_params(labelsize=7)
            # Mark the best track
            if if_btk_exist:
                ax0.scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 5, 'red', marker='*',transform=ccrs.PlateCarree())

            # Set labels
            lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
            lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
            gl = ax0.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
            gl.xlabels_top = False
            gl.xlabels_bottom = True
            gl.ylabels_left = True
            gl.ylabels_right = False
            gl.ylocator = mticker.FixedLocator(lat_ticks)
            gl.xlocator = mticker.FixedLocator(lon_ticks)
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 10}
            gl.ylabel_style = {'size': 10}

            # Set titles
            ax0.set_title('Flight Track at '+ it_str,fontweight="bold",fontsize='12')

            #  --------- Compare the flight-level obs with model output  --------------------
            time_start = it_dt - timedelta(minutes=Nminutes)
            time_end = it_dt + timedelta(minutes=Nminutes)

            ax1 = fig.add_subplot( gs[3,:] ) 
            ax1.plot_date( dict_AF_hour['GMT'], dict_AF_hour['WSpd']*0.51444,  color='black', linewidth=1, label='Flight' ) 
            ax1.plot_date( dict_AF_hour['GMT'][idx_itp], np.array( ws_AFspace ),  color='blue', linewidth=1, label='IR+MW')

            ax2 = fig.add_subplot( gs[4,:] )
            ax2.plot_date( dict_AF_hour['GMT'], dict_AF_hour['TD']+273.15,  color='black', linewidth=1, label='Flight' )
            ax2.plot_date( dict_AF_hour['GMT'][idx_itp], np.array( td_AFspace )+273.15, color='blue', linewidth=1, label='IR+MW')

            # Set labels
            ax1.tick_params(left = True, right = False , labelleft = True ,
                                            labelbottom = False, bottom = False, labelsize='9')
            ax2.tick_params( labelsize='9')
            # Set titles
            ax1.set_title('Flight-level Wind (m/s)',fontweight="bold",fontsize='12')
            ax2.set_title('Flight-level Dew Point (k)',fontweight="bold",fontsize='12')

            plt.savefig( plot_dir+it_str+'ws_td_approx.png', dpi=300 )
            print('Saving the figure: ',  plot_dir+it_str+'_track_ws_td_approx.png')
            plt.close()


    eng.quit()



# For the wrfout in each deterministic forecast, plot the aircraft flight track and the comparison of the wind speed/dwt 
def DF_compare_track_wind_dew( Storm, Exper_name, big_dir, small_dir ):

    DFtimes = ['201708221200','201708221800','201708230000','201708230600', '201708231200']#,'201709171200','201709171800',]
    DFend = '201708270000'


    for DF_start in DFtimes:
        print('The deterministic forecast begins at '+ DF_start)
        df_dir = big_dir+Storm+'/'+Exper_name+'/wrf_df/'+DF_start+'/'
        plot_dir = df_dir + 'compare_flight_level/'
        
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            Plot_compare_track_wind_td( Storm, small_dir, df_dir, plot_dir )
        else:
            Plot_compare_track_wind_td( Storm, small_dir, df_dir, plot_dir )


class Info_USAF_clt:
    '''This object is designed to collect information for each wrf day'''
    def __init__(self, filedir, beg_time, end_time):
        self.filedir = filedir
        self.beg_time = beg_time
        self.end_time = end_time


if __name__ == '__main__':
    Storm = 'HARVEY'
    Exper_name = 'newWRF_MW_THO'
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/' 

    #UV10_slp_AF( Storm, Exper_name, big_dir, small_dir )
    DF_compare_track_wind_dew( Storm, Exper_name, big_dir, small_dir )

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
import matlab.engine
import scipy as sp
import scipy.ndimage
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time

import Util_data as UD

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



# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot UV10_slp_USAF data
# ------------------------------------------------------------------------------------------------------


def plot_UV10_slp( Storm, DAtime, wrfout, plot_dir, dict_AF_masked):

    # Read storm center
    dict_btk = UD.read_bestrack(Storm)
    # Find the best-track position
    btk_dt = [it_str for it_str in dict_btk['time'] ]#[datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]
    bool_match = [DAtime == it for it in btk_dt]
    if True in bool_match:
        if_btk_exist = True
        idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
    else:
        if_btk_exist = False

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
        gpsa_min = 0 #np.amin( dict_AF_masked['gpsa'][:] )
        gpsa_max = 10000# np.amax( dict_AF_masked['gpsa'][:] )
        AFAF = ax.scatter(dict_AF_masked['LON'][:],dict_AF_masked['LAT'][:],3,dict_AF_masked['GPSA'][:],cmap='Greens_r',transform=ccrs.PlateCarree(),vmin=gpsa_min, vmax=gpsa_max)
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

    plt.savefig( plot_dir+wrfout_head_tail[1]+'_'+DAtime+'.png', dpi=300 )
    print('Saving the figure: ', plot_dir+wrfout_head_tail[1]+'_'+DAtime+'.png')
    plt.close()

def UV10_slp_AF( Storm, Exper_name, DAtimes, big_dir, small_dir, USAF_Tspan ):

    # Attributes in USAF of interest
    attr_interest_usaf = ['GMT','GPSA','LAT','LON']

    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
        print('Reading WRF analysis: ', wrf_dir)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )

        # ------- check if at DAtime USAF obs exists --------------
        Use_USAF = []
        for imission, iTspan in USAF_Tspan.items():
            USAF_start = iTspan[0]
            USAF_end = iTspan[1]
            if DAtime_dt <= USAF_end and DAtime_dt >= USAF_start: # there might be a bug here !
                Use_USAF.append( imission )
                break # terminate the loop once the mission is found
            else:
                Use_USAF.append( None )

        if all(i is None for i in Use_USAF): 
            print('No USAT obs exists!')
            dict_AF_masked = None
        else:
            imission = [i is not None for i in Use_USAF].index(True)
            mission_path = Use_USAF[imission]
            print('By examining the files, USAF data exist in : ', mission_path)
            dict_AF_all_today = USAF.read_USAF_mission(mission_path, attr_interest_usaf)
            # Mask time
            DA_hh = int( DAtime[8:10] )
            Nminutes = 10
            dict_AF_masked = USAF.mask_time( dict_AF_all_today, DA_hh, Nminutes)
        
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'//Vis_analyze/Model/UV10_slp_AF/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            plot_UV10_slp( Storm, DAtime, wrf_dir, plot_dir, dict_AF_masked )
        else:
            plot_UV10_slp( Storm, DAtime, wrf_dir, plot_dir, dict_AF_masked )



# ------------------------------------------------------------------------------------------------------
#           Operation: Interpolation from model variables to obs locations
# ------------------------------------------------------------------------------------------------------

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
    #H_left_bot = H_x[:,j_PlusOne,i_MinusOne]
    H_left_bot = H_x[:,j_MinusOne,i_MinusOne]
    H_left_up = H_x[:,j_PlusOne,i_MinusOne]
    #H_right_bot = H_x[:,j_PlusOne,i_PlusOne]
    H_right_bot = H_x[:,j_PlusOne,i_PlusOne]
    H_right_up = H_x[:,j_MinusOne,i_PlusOne]
    
    # perform the bilinear interpolation over an unit area
    area = (i_PlusOne-i_MinusOne)*(j_PlusOne-j_MinusOne)
    if area != 1:
        raise Exception('The area is not a unit square!') # raise an error and stop the program
    
    fra_i_prime = math.modf(i_prime)[1]
    fra_j_prime = math.modf(j_prime)[1]
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
    if ( z_prime - np.fix(z_prime)) == 0: # in case The underlying triangulation is empty - the points may be coplanar or collinear 
        z_prime = z_prime + 0.0000001

    print('The location of the obs in the ijk coordinate is: ', i_prime, j_prime, z_prime)
    return [i_prime, j_prime, z_prime]


def interp_v_to_obs_matlab( wrf_out, dict_AF_hour ):

    print("Initiate the function to interpolate model variables to obs locations...")
    start_time=time.process_time() 
    
    eng = matlab.engine.start_matlab() # start a new matlab process
 
    # Read the wrf variables
    print('Check wrf file: ', wrf_out)
    ncdir = nc.Dataset( wrf_out )
    lat_x = ncdir.variables['XLAT'][0,:,:]#[::10,::10]
    lon_x = ncdir.variables['XLONG'][0,:,:]#[::10,::10]

    u = ncdir.variables['U'][0,:,:,:]
    u_mass = (u[:,:,0:-1:1] + u[:,:,1::1])/2 # interpolate staggered point to the mass grid 
    v = ncdir.variables['V'][0,:,:,:]
    v_mass = (v[:,0:-1:1,:] + v[:,1::1,:])/2 # interpolate staggered point to the mass grid 
    ws_x =  ((u_mass ** 2 + v_mass ** 2) ** 0.5)#[::4,::10,::10]

    PHB = ncdir.variables['PHB'][0,:,:,:]
    PH = ncdir.variables['PH'][0,:,:,:]
    GP = PHB+PH #bottom_top, south_north, west_east
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
        return # exit the function
    if (np.amax(lat_f) < np.amin(lat_x)) or (np.amin(lat_f) > np.amax(lat_x)):
        print('The flight obs is not within the model domain!')
        return # exit the function

    idx_itp = []
    ws_AFspace = []
    td_AFspace = []
    for iobs in range( len(dict_AF_hour['LON']) ):
        print('------------------------------------------------')
        print('GMT: ', dict_AF_hour['GMT'][iobs])
        # get the ijk of obs
        list_ijk_obs = ijk_obs( lon_x, lat_x, Height_mass_x, lon_f[iobs], lat_f[iobs], gpsa_f[iobs] )
        if list_ijk_obs is None:
            print('Interpolation is not applicable here!')
            return

        # collect the ijk coordinates of the 8 points  
        i_PlusOne = int( np.ceil( list_ijk_obs[0] ) )
        i_MinusOne = int( np.floor( list_ijk_obs[0] ) )
        j_PlusOne = int( np.ceil( list_ijk_obs[1] ) )
        j_MinusOne = int( np.floor( list_ijk_obs[1] ) )
        k_PlusOne = int( np.ceil( list_ijk_obs[2] ) )
        k_MinusOne = int( np.floor( list_ijk_obs[2] ) )
        i_8_points = [i_MinusOne,i_PlusOne,i_PlusOne,i_MinusOne,i_MinusOne,i_PlusOne,i_PlusOne,i_MinusOne]
        j_8_points = [j_MinusOne,j_MinusOne,j_PlusOne,j_PlusOne,j_MinusOne,j_MinusOne,j_PlusOne,j_PlusOne]
        k_8_points = [k_MinusOne,k_MinusOne,k_MinusOne,k_MinusOne,k_PlusOne,k_PlusOne,k_PlusOne,k_PlusOne]
        # griddata interpolate
        ws_8_points = []
        td_8_points = []
        for it in range(len(i_8_points)):
            ws_8_points.append( ws_x[k_8_points[it], j_8_points[it], i_8_points[it]] )
            td_8_points.append( td_x[k_8_points[it], j_8_points[it], i_8_points[it]] )

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

    eng.quit()

    # Stack each list into an array
    all_times_str = [ datetime.strftime(it, "%Y%m%d%H%M%S") for it in dict_AF_hour['GMT'] ] # convert date object to string
    Wspd = [ it*0.51444 for it in dict_AF_hour['WSpd'] ]
    all_attrs = np.column_stack( (all_times_str, dict_AF_hour['LAT'], dict_AF_hour['LON'], dict_AF_hour['GPSA'], Wspd, ws_AFspace, dict_AF_hour['TD'], td_AFspace) )

    # ---- Write to file and save it to the disk ----
    header = ['GMT', 'Lat','Lon','GPSA','WS_obs','WS_interp','Td_obs','Td_interp']
    file_name = wrf_out.replace('wrf','Interp_USAF') #Hx_dir + "/mean_obs_res_d03" + DAtime + '_' +  sensor + '.txt'
    print('Saving interpolated variables to ', file_name)
    with open(file_name,'w') as f:
        # Add header 
        f.write('\t'.join( item.rjust(8) for item in header ) + '\n' )
        # Write the record to the file serially
        len_records = np.shape( all_attrs )[0]
        for irow in range( len_records ):
            #iformat = [ "{0:.4f}".format(item) for item in all_attrs[irow,1:] ]
            irecord =  [str(item) for item in all_attrs[irow,:] ]
            f.write('\t'.join( item.rjust(8) for item in irecord ) + '\n')

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')


def Pre_interp_v_to_obs_matlab( Storm, Exper_name, DAtimes, USAF_Tspan):

    # Attributes in USAF of interest
    attrs_itt = ['GMT','GPSA','LAT','LON','WSpd','TD'] 

    # Loop through each DAtime/analysis
    dict_interpolated = {}
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
        print('Reading WRF analysis: ', wrf_dir)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )

        # ------- check if at DAtime USAF obs exists --------------
        Use_USAF = []
        for imission, iTspan in USAF_Tspan.items():
            USAF_start = iTspan[0]
            USAF_end = iTspan[1]
            if DAtime_dt <= USAF_end and DAtime_dt >= USAF_start: # there might be a bug here !
                Use_USAF.append( imission )
                break # terminate the loop once the mission is found
            else:
                Use_USAF.append( None )

        if all(i is None for i in Use_USAF):
            print('No USAT obs exists!')
            dict_AF_masked = None
            continue
        else:
            imission = [i is not None for i in Use_USAF].index(True)
            mission_path = Use_USAF[imission]
            print('By examining the files, USAF data exist in : ', mission_path)
            dict_AF_all_today = USAF.read_USAF_mission(mission_path, attrs_itt)
            # Mask time
            DA_hh = int( DAtime[8:10] )
            Nminutes = 10 #!!!!!!!!!!!!!!!!!!!!
            dict_AF_masked = USAF.mask_time( dict_AF_all_today, DA_hh, Nminutes)
            if dict_AF_masked == None:
                print('No AF obs is available at this time span!')
                continue
            else:
                interp_v_to_obs_matlab( wrf_dir, dict_AF_masked ) 

    

# ------------------------------------------------------------------------------------------------------
#           Operation: Read and Plot Interpolated model variables at obs locations
# ------------------------------------------------------------------------------------------------------

# Read interpolated variables
def read_interpolation( Ifile ):

    print('Reading interpolated file: ', Ifile)
    time_obs = []
    lat_obs = []
    lon_obs = []
    gpsa_obs = []
    ws_obs = []
    interp_ws_obs = []
    td_obs = []
    interp_td_obs = []

    # Read records
    with open( Ifile ) as f:
        next(f)
        all_lines = f.readlines() 

    for line in all_lines:
        split_line = line.split()
        time_obs.append(  datetime.strptime(split_line[0], "%Y%m%d%H%M%S") )
        lat_obs.append(  float(split_line[1]) )        
        lon_obs.append( float(split_line[2]) )
        gpsa_obs.append( float(split_line[3]) )
        ws_obs.append( float(split_line[4]) )
        interp_ws_obs.append( float(split_line[5]) )
        td_obs.append( float(split_line[6]) )
        interp_td_obs.append( float(split_line[7]) )

    time_obs = np.array( time_obs )
    lat_obs = np.array( lat_obs )
    lon_obs = np.array( lon_obs )
    gpsa_obs = np.array( gpsa_obs )
    ws_obs = np.array( ws_obs )
    interp_ws_obs = np.array( interp_ws_obs )
    td_obs = np.array( td_obs )
    interp_td_obs = np.array( interp_td_obs )

    d_one = {'time_obs':time_obs, 'lat_obs':lat_obs, 'lon_obs':lon_obs, 'gpsa_obs': gpsa_obs, 'ws_obs':ws_obs, 'interp_ws_obs':interp_ws_obs, 'td_obs':td_obs, 'interp_td_obs':interp_td_obs}
    return d_one

# Plot hourly flight track + variable values of model output at XX UTC to the flight level points (lon,lat,height) within XX UTC +/- 10 minutes
def Plot_compare_track_wind_td( Storm, Exper_name, DAtime, idx_exist, dict_AF_mission, big_dir, plot_dir):
    
    # Read storm center
    dict_btk = UD.read_bestrack(Storm)
    # Find the best-track position
    btk_dt = [it_str for it_str in dict_btk['time'] ]#[datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]
    bool_match = [DAtime == it for it in btk_dt]
    if True in bool_match:
        if_btk_exist = True
        idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
    else:
        if_btk_exist = False

    # Read interpolated USAF obs info
    d_all = {}
    for iExper in Exper_name:
        if idx_exist[iExper][DAtime]:
            file_dir = big_dir+Storm+'/'+iExper+'/fc/'+DAtime+'/Interp_USAF_enkf_output_d03_mean'
            d_all[iExper] = read_interpolation( file_dir )
        else:
            d_all[iExper] = None
    
    #  -------------------------- Plot ---------------------------------------------
    fig = plt.figure( figsize=(8,9.5), dpi=150 )
    gs = fig.add_gridspec(5,2)
           
    #  --------- Plot the flight level track  --------------------
    ax0 = fig.add_subplot( gs[0:3,:],  projection=ccrs.PlateCarree())
           
    lon_min = -100#-71#np.amin( dict_AF_mission['lon'] )
    lon_max = -85#-57#np.amax( dict_AF_mission['lon'] )
    lat_min = 19#10#np.amin( dict_AF_mission['lat'] )
    lat_max = 31#20#np.amax( dict_AF_mission['lat'] )

    ax0.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax0.coastlines (resolution='10m', color='black', linewidth=1)
    #  Plot the track for the whole mission
    gpsa_min = 0
    gpsa_max = 10000
    ax0.scatter(dict_AF_mission['LON'], dict_AF_mission['LAT'], 1, 'grey', cmap='jet', vmin=gpsa_min, vmax=gpsa_max,transform=ccrs.PlateCarree())
    for iExper in Exper_name:
        if d_all[iExper] is not None:
            lon_obs = d_all[iExper]['lon_obs']
            lat_obs = d_all[iExper]['lat_obs']
            gpsa_obs = d_all[iExper]['gpsa_obs']
            AF = ax0.scatter(lon_obs, lat_obs, 2, gpsa_obs, cmap='jet', vmin=gpsa_min, vmax=gpsa_max, transform=ccrs.PlateCarree())
    AF_bar = fig.colorbar(AF,ax=ax0,shrink=0.9)
    AF_bar.ax.set_ylabel('Flight Height (m)', fontsize=12)
    AF_bar.ax.tick_params(labelsize=12)
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
    ax0.set_title('Flight Track at '+ DAtime,fontweight="bold",fontsize='12')

    #  --------- Compare the flight-level obs with model output  --------------------
    colors = ['red','blue']    
    labels = ['IR', 'IRMW']

    i = 0
    ax1 = fig.add_subplot( gs[3,:] )
    ax2 = fig.add_subplot( gs[4,:] )
    ax1.plot_date( d_all[iExper]['time_obs'], d_all[iExper]['ws_obs'],  color='black', linewidth=1, label='Flight')
    ax2.plot_date( d_all[iExper]['time_obs'], d_all[iExper]['td_obs'],  color='black', linewidth=1, label='Flight')
    for iExper in Exper_name:
        if d_all[iExper] is not None:
            ax1.plot_date( d_all[iExper]['time_obs'], d_all[iExper]['interp_ws_obs'],  color=colors[i], linewidth=1, label=labels[i])
            ax2.plot_date( d_all[iExper]['time_obs'], d_all[iExper]['interp_td_obs'], color=colors[i], linewidth=1, label=labels[i])
        i = i+1

    ax2.legend( bbox_to_anchor=(0.55,-0.15), fontsize = 12)
    # Set x/y ticks
    ax1.tick_params(left = True, right = False, bottom = False, labelbottom = False, labelsize=10)
    ax2.tick_params( labelsize=10)
    # Set titles
    ax1.set_title('Flight-level Wind (m/s)',fontweight="bold",fontsize=12)
    ax2.set_title('Flight-level Dew Point (c)',fontweight="bold",fontsize=12)

    plt.savefig( plot_dir+DAtime+'_track_ws_td.png', dpi=300 )
    print('Saving the figure: ',  plot_dir+DAtime+'_track_ws_td.png')
    plt.close()



# Plot the aircraft flight track and the comparison of the wind speed/dwt 
def DF_compare_track_wind_dew( Storm, Exper_name, DAtimes, big_dir, small_dir, USAF_Tspan ):

    # Set up directory
    plot_dir = small_dir+Storm+'/'+Exper_name[-1]+'/Vis_analyze/Model/USAFspace/'
    plotdir_exists = os.path.exists( plot_dir )
    if plotdir_exists == False:
        os.mkdir(plot_dir)

    # Get the whole mission data for all DAtimes if available
    attr_interest_usaf = ['GMT','GPSA','LAT','LON']
    dict_AF_allts = {}
    for DAtime in DAtimes:
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------- check if at DAtime USAF obs exists --------------
        Use_USAF = []
        for imission, iTspan in USAF_Tspan.items():
            USAF_start = iTspan[0]
            USAF_end = iTspan[1]
            if DAtime_dt <= USAF_end and DAtime_dt >= USAF_start: # there might be a bug here !
                Use_USAF.append( imission )
                break # terminate the loop once the mission is found
            else:
                Use_USAF.append( None )

        if all(i is None for i in Use_USAF):
            print('No USAT obs exists!')
            dict_AF_masked = None
        else:
            imission = [i is not None for i in Use_USAF].index(True)
            mission_path = Use_USAF[imission]
            print('By examining the files, USAF data exist in : ', mission_path)
            dict_AF_all_today = USAF.read_USAF_mission(mission_path, attr_interest_usaf)
            dict_AF_allts[DAtime] = dict_AF_all_today
    
    # Get the generated interpolation
    idx_exist = {}
    # Loop through each experiment
    for iExper in Exper_name:
        idx_exist[iExper] = {}
        for itime in DAtimes:
            file_dir = big_dir+Storm+'/'+iExper+'/fc/'+itime+'/Interp_USAF_enkf_output_d03_mean' 
            idx_lg = os.path.exists( file_dir ) 
            idx_exist[iExper][itime] = idx_lg

    # Loop through each DAtime/analysis
    for it in range(len(DAtimes)):
        if any([idx_exist[iExper][DAtimes[it]] for iExper in Exper_name]): # if the interpolation exists under any experiment on this time
            Plot_compare_track_wind_td( Storm, Exper_name, DAtimes[it], idx_exist, dict_AF_allts[DAtimes[it]], big_dir, plot_dir)  
        else:
            continue


if __name__ == '__main__':
    
    Storm = 'HARVEY'
    Exper_name = ['JerryRun/MW_THO/',]
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    Interp_to_obs = False
    Plot_UV10_slp_AF = True

    # Time range set up
    start_time_str = '201708230900'
    end_time_str = '201708231200'
    Consecutive_times = False

    if not Consecutive_times:
        DAtimes = ['201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt] 

    # Check the time_start and time_end for each USAF file
    USAF_Tspan = {}
    USAF_list = sorted(glob.glob( small_dir + Storm + '/USAF/201708*' ))
    for imission in USAF_list:
        print( '--------- The mission is ', imission, '--------------' )
        time_start, time_end = USAF.Time_range_USAF_mission( imission )
        USAF_Tspan[imission] = [time_start, time_end] # date object

    # Plot low-level circulation and AF location
    if Plot_UV10_slp_AF:
        for iExper in Exper_name:
            UV10_slp_AF( Storm, iExper, DAtimes, big_dir, small_dir, USAF_Tspan )

    # Plot interpolated model variables at USAF obs location
    # Interpolate variables in model resolution to obs location AND write it to a txt file
    #if Interp_to_obs:
    #    for iExper in Exper_name:
    #        print('------------ Interpolate variables at model resolution to obs location --------------')
    #        Pre_interp_v_to_obs_matlab( Storm, iExper, DAtimes, USAF_Tspan)

    #DF_compare_track_wind_dew( Storm, Exper_name, DAtimes, big_dir, small_dir, USAF_Tspan  )



















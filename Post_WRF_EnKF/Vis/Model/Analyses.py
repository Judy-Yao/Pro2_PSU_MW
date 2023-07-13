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
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time

import Util_data as UD

# setting font sizeto 30
plt.rcParams.update({'font.size': 15})

def d03_domain( wrfout_d03 ):
    ncdir = nc.Dataset(wrfout_d03, 'r')

    xlat = ncdir.variables['XLAT'][0,:,:]
    xlong = ncdir.variables['XLONG'][0,:,:]

    d03_lat_min = np.min( xlat.flatten() )
    d03_lat_max = np.max( xlat.flatten() )
    d03_lon_min = np.min( xlong.flatten() )
    d03_lon_max = np.max( xlong.flatten() )

    d03_list = [d03_lon_min, d03_lon_max, d03_lat_min, d03_lat_max]
    return d03_list

# !!!Judy wrote it and it was wrong!!!!
def cal_IC_water_way1( ncdir ):
    # read water vapor mixing ratio
    qvapor_3d = ncdir.variables['QVAPOR'][0,:,:,:] #bottom_top, south_north, west_east
    qvapor_3d_TtoB = np.flip( qvapor_3d, axis=0 )
    # calculated the water vapor mixting ratio at the grid center (staggered)
    qvapor_mass_center = (qvapor_3d_TtoB[1:,:,:] + qvapor_3d_TtoB[0:np.size(qvapor_3d,0)-1,:,:])/2
    # use the assumption (hydrostatic balance) to calculate IC water
    PB = ncdir.variables['PB'][0,:,:,:]
    P = ncdir.variables['P'][0,:,:,:]
    P_3d = PB + P
    P_3d_TtoB = np.flip( P_3d,axis=0 )
    delta_P = P_3d_TtoB[1:,:,:] - P_3d_TtoB[0:np.size(qvapor_3d,0)-1,:,:]
    vapor_content = qvapor_mass_center*delta_P;
    IC_water_file = np.sum(vapor_content,axis=0)/9.8
    return IC_water_file



# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot precipitable water per snapshot
# ------------------------------------------------------------------------------------------------------
def read_IC_water( wrf_dir ):

    # set file names for background and analysis
    mean_dir = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
    # read the shared variables: lat/lon
    ncdir = nc.Dataset( mean_dir[0] )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # calculate IC water
    IC_water = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] ) 
    for i in range(len(mean_dir)):
        file_name = mean_dir[i]
        ncdir = nc.Dataset( file_name )
        IC_water[i,:,:] = getvar(ncdir, 'pw')

    d_IC_water = {'lat':lat,'lon':lon,'IC_water':IC_water}
    return d_IC_water

def plot_IC_water( Storm, Exper_name, DAtime, wrf_dir, plot_dir ):

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

    #------ Read WRFout -------------------
    d_icwv = read_IC_water( wrf_dir ) 
    lat = d_icwv['lat']
    lon = d_icwv['lon']
    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )

    # ------------------ Plot -----------------------
    fig, axs=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(12,5), dpi=400)

    # Xb
    min_icwv = 30
    max_icwv = 65
    axs.flat[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    axs.flat[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    xb_wv = axs.flat[0].scatter(lon,lat,1.5,c=d_icwv['IC_water'][0,:,:],edgecolors='none', cmap='ocean_r', vmin=min_icwv, vmax=max_icwv,transform=ccrs.PlateCarree())
    # Mark the best track
    if if_btk_exist:
        axs.flat[0].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'white', marker='*',transform=ccrs.PlateCarree())

    # Xa
    axs.flat[1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    axs.flat[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    xa_wv = axs.flat[1].scatter(lon,lat,1.5,c=d_icwv['IC_water'][1,:,:],edgecolors='none', cmap='ocean_r', vmin=min_icwv, vmax=max_icwv,transform=ccrs.PlateCarree())
    # Mark the best track
    if if_btk_exist:
        axs.flat[1].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'white', marker='*',transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.12, 0.1, 0.5, 0.02])
    xwv_bar = fig.colorbar(xb_wv,ax=axs[0:2],orientation="horizontal", cax=caxes)
    xwv_bar.ax.tick_params()

    # Xa-Xb (increment)
    min_incre = -10
    max_incre = 10
    axs.flat[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    axs.flat[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    icwv_incre = d_icwv['IC_water'][1,:,:] - d_icwv['IC_water'][0,:,:]
    incre_wv = axs.flat[2].scatter(lon,lat,1.5,c=icwv_incre,edgecolors='none', cmap='bwr', vmin=min_incre, vmax=max_incre,transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.65, 0.1, 0.25, 0.02])
    cb_diff_ticks = np.linspace(min_incre, max_incre, 5, endpoint=True)
    cbar = fig.colorbar(incre_wv, ax=axs[2:], ticks=cb_diff_ticks, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params()

    #subplot title
    axs.flat[0].set_title('Xb--ICWV ($\mathregular{kgm^{-2}}$)',fontweight='bold')
    axs.flat[1].set_title('Xa--ICWV ($\mathregular{kgm^{-2}}$)',fontweight='bold')
    axs.flat[2].set_title('Xa-Xb--ICWV ($\mathregular{kgm^{-2}}$)', fontweight='bold')

    #title for all
    fig.suptitle(Storm+': '+Exper_name+'('+DAtime+')', fontsize=13, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2))
    for j in range(3):
        gl = axs.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
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

    plt.savefig( plot_dir+DAtime+'_ICWV.png', dpi=300 )
    print('Saving the figure: ', plot_dir+DAtime+'_ICWV.png')
    plt.close()


def IC_water( Storm, Exper_name, DAtimes, big_dir, small_dir ):

    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
        print('Reading WRF background and analysis: ', wrf_dir)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/IC_water/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            plot_IC_water( Storm, Exper_name, DAtime, wrf_dir, plot_dir )
        else:
            plot_IC_water( Storm, Exper_name, DAtime, wrf_dir, plot_dir )

# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot the evolution of min slp
# ------------------------------------------------------------------------------------------------------
def find_minSLP( wrfout ):
    ncdir = nc.Dataset( wrfout )
    slp = getvar(ncdir, 'slp')
    min_slp = np.amin( slp )
    slp_smooth = sp.ndimage.gaussian_filter(slp, [11,11])
    idx = np.nanargmin( slp_smooth )
    lat_minslp = ncdir.variables['XLAT'][:].flatten()[idx]
    lon_minslp = ncdir.variables['XLONG'][:].flatten()[idx]

    minSLP = [lat_minslp,lon_minslp,min_slp.values]
    return minSLP

def plot_slp_timeseries( small_dir, Storm, Expers, DAtimes, Evo_slp ):

    # Set up figure
    fig = plt.figure( figsize=(12,9), dpi=300 )
    ax = plt.subplot(1,1,1)
    dates = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]
    # Plot obs and mean
    print(Evo_slp['obs_slp'])
    ax.plot_date(dates, Evo_slp['obs_slp'], 'black', label='TCvital', linestyle='-')
    ax.plot_date(dates, Evo_slp['xa_slp'][0,:], 'red', label='conv_THO', linestyle='--')

    leg = plt.legend(loc='lower left',fontsize=15)
    # Set X/Y labels
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45,labelsize=15)
    ax.set_ylabel('Minimum Sea Level Pressure (hPa)',fontsize=15)
    ax.set_ylim( 900,980 ) #1010,1012

    #title_name = Storm+'('+Exper_name+')'+': mim SLP'
    title_name = Storm+': mim SLP'
    ax.set_title( title_name,fontweight="bold",fontsize='15' )
    #fig.suptitle('conv+HPI', fontsize=15, fontweight='bold')


    # Save the figure
    save_des = small_dir+Storm+'/'+Expers[0]+'/Vis_analyze/Model/minslp_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()

def Gather_slp( Storm, Expers, DAtimes, big_dir ):
    
    obs_minslp_lat = []
    obs_minslp_lon = []
    obs_minslp_value = []
    
    xa_minslp_lat = [[] for i in range(len(Expers))]
    xa_minslp_lon = [[] for i in range(len(Expers))]
    xa_minslp_value = [[] for i in range(len(Expers))]

    for DAtime in DAtimes:

        # collect min slp found from WRF output
        for Exper in Expers:
            wrf_dir = big_dir+Storm+'/'+Exper+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
            print('Reading the EnKF posterior mean from ', wrf_dir)
            list_wrfout = find_minSLP( wrf_dir )
            idx = Expers.index( Exper )
            xa_minslp_lat[idx].append( list_wrfout[0] )
            xa_minslp_lon[idx].append( list_wrfout[1] )
            xa_minslp_value[idx].append( list_wrfout[2] )

        # collect assimilated min slp obs from TCvital record in fort.10000
        diag_enkf = big_dir+Storm+'/'+Expers[0]+'/run/'+DAtime+'/enkf/d03/fort.10000'
        print('Reading the EnKF diagnostics from ', diag_enkf)
        enkf_minSlp = subprocess.run(["grep","slp",diag_enkf],stdout=subprocess.PIPE,text=True)
        list_enkf_minSlp = enkf_minSlp.stdout.split()
        # condition on the number of assimilated min slp
        if list_enkf_minSlp.count('slp') == 0 :
            raise ValueError('No min slp is assimilated!')
        elif list_enkf_minSlp.count('slp') == 1 :
            obs_minslp_lat.append( float(list_enkf_minSlp[1]) )    
            obs_minslp_lon.append( float(list_enkf_minSlp[2]) )  
            obs_minslp_value.append( float(list_enkf_minSlp[9])/100 )
        else : # at least two min slp records are assimilated
            print('At least two min slp obs are assimilated!')
            # find the index/location of 'slp' in fort.10000
            indices = [i for i ,e in enumerate(list_enkf_minSlp) if e == 'slp']
            # assemble a pair of coordinate for each 'slp'
            coor_pair = [ np.array([float(list_enkf_minSlp[it+1]),float(list_enkf_minSlp[it+2])]) for it in indices]
            # find the index of the current time
            idx_time = DAtimes.index( DAtime )
            # assemble the diagnosed min slp from an analysis
            model_minslp = np.array([xa_minslp_lat[0][idx_time],xa_minslp_lon[0][idx_time]] )
            # find the nearest TCvital min slp from the analysis
            distances = [ np.sum(np.square(coor_pair[it], model_minslp )) for it in range(len(coor_pair))]
            min_distances = [np.amin(distances) == it for it in distances]
            idx_coor = min_distances.index(True)
            # gather this TCvital min slp
            obs_minslp_lat.append( float(list_enkf_minSlp[indices[idx_coor]+1]) )
            obs_minslp_lon.append( float(list_enkf_minSlp[indices[idx_coor]+2]) )
            obs_minslp_value.append( float(list_enkf_minSlp[indices[idx_coor]+9])/100 )
 

    dict_minSLP = {'obs_slp':np.array(obs_minslp_value),'obs_lat':np.array(obs_minslp_lat),'obs_lon':np.array(obs_minslp_lon),'xa_slp':np.array(xa_minslp_value),'xa_lat':np.array(xa_minslp_lat), 'xa_lon':np.array(xa_minslp_lon)}
    return dict_minSLP
    

# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot UV10_slp per snapshot
# ------------------------------------------------------------------------------------------------------
def read_UV10_slp( wrf_dir ):

    # set file names for background and analysis
    mean_dir = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
    # read the shared variables: lat/lon
    ncdir = nc.Dataset( mean_dir[0] )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # read UV10 and slp
    slp_original = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    slp_smooth = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    U10 = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    V10 = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    windspeed = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    for i in range(len(mean_dir)):
        file_name = mean_dir[i]
        ncdir = nc.Dataset( file_name )
        # sea level pressure
        slp_original[i,:,:] = getvar(ncdir, 'slp')
        slp_smooth[i,:,:] = sp.ndimage.gaussian_filter(slp_original[i,:,:], [11,11])
        # UV10
        U10[i,:,:] = ncdir.variables['U10'][0,:,:]
        V10[i,:,:] = ncdir.variables['V10'][0,:,:]
        windspeed[i,:,:] = (U10[i,:,:] ** 2 + V10[i,:,:] ** 2) ** 0.5

    d_UV10_slp = {'lat':lat,'lon':lon,'slp':slp_original,'slp_smooth':slp_smooth,'U10':U10,'V10':V10,'windspeed':windspeed}
    return d_UV10_slp

def plot_UV10_slp( Storm, Exper_name, DAtime, wrf_dir, plot_dir ):

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
    d_field = read_UV10_slp( wrf_dir )
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

    for i in range(2):
        ax[i].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i].coastlines (resolution='10m', color='black', linewidth=1)
        # sea level pressure
        start_level = max( np.amin(slp_smooth[0,:,:]), np.amin(slp_smooth[1,:,:]) )+0.1
        end_level = min( np.amax(slp_smooth[0,:,:]), np.amax(slp_smooth[1,:,:]) )
        step = (end_level -start_level)/5
        level = np.arange(start_level, end_level,step)
        #slp_contour = ax[i].contour(lon,lat,slp_smooth[i,:,:],levels=level,cmap='Greys_r',vmin=min_slp,vmax=max_slp,transform=ccrs.PlateCarree())
        slp_contour = ax[i].contour(lon,lat,slp_smooth[i,:,:],levels=level,cmap='Greys_r',transform=ccrs.PlateCarree())
        plt.clabel(slp_contour,level,inline=True, fmt="%.3f", use_clabeltext=True) #, fontsize=12)
        # Wind at 10 meters
        wind_smooth = sp.ndimage.gaussian_filter( d_field['windspeed'][i,:,:], [2,2])
        min_wind = 0
        max_wind = 35
        bounds = np.linspace(min_wind, max_wind, 6)
        wind_contourf = ax[i].contourf(lon,lat,wind_smooth,cmap='hot_r',vmin=min_wind,vmax=max_wind,levels=bounds,extend='both',transform=ccrs.PlateCarree())
        ax[i].barbs(lon.flatten(), lat.flatten(), d_field['U10'][i,:,:].flatten(), d_field['V10'][i,:,:].flatten(), length=5, pivot='middle', color='royalblue', regrid_shape=20, transform=ccrs.PlateCarree())
        # Mark the best track
        if if_btk_exist:
            ax[i].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'green', marker='*',transform=ccrs.PlateCarree())

    # Adding the colorbar
    cbaxes = fig.add_axes([0.01, 0.1, 0.03, 0.8])
    wind_bar = fig.colorbar(wind_contourf,cax=cbaxes,fraction=0.046, pad=0.04) #Make a colorbar for the ContourSet returned by the contourf call.
    wind_bar.ax.set_ylabel('Wind Speed (m/s)')
    wind_bar.ax.tick_params(labelsize='15')

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

    plt.savefig( plot_dir+DAtime+'_UV10_slp.png', dpi=300 )
    print('Saving the figure: ', plot_dir+DAtime+'_UV10_slp.png')
    plt.close()

def UV10_slp( Storm, Exper_name, DAtimes, big_dir, small_dir ):

    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
        print('Reading WRF background and analysis: ', wrf_dir)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/UV10_slp/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            plot_UV10_slp( Storm, Exper_name, DAtime, wrf_dir, plot_dir )
        else:
            plot_UV10_slp( Storm, Exper_name, DAtime, wrf_dir, plot_dir )

# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot precipitable water per snapshot
# ------------------------------------------------------------------------------------------------------
def read_IC_water( wrf_dir ):

    # set file names for background and analysis
    mean_dir = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
    # read the shared variables: lat/lon
    ncdir = nc.Dataset( mean_dir[0] )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # calculate IC water
    IC_water = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    for i in range(len(mean_dir)):
        file_name = mean_dir[i]
        ncdir = nc.Dataset( file_name )
        IC_water[i,:,:] = getvar(ncdir, 'pw')

    d_IC_water = {'lat':lat,'lon':lon,'IC_water':IC_water}
    return d_IC_water

def plot_IC_water( Storm, Exper_name, DAtime, wrf_dir, plot_dir ):

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

    #------ Read WRFout -------------------
    d_icwv = read_IC_water( wrf_dir )
    lat = d_icwv['lat']
    lon = d_icwv['lon']
    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )

    # ------------------ Plot -----------------------
    fig, axs=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(12,5), dpi=400)

    # Xb
    min_icwv = 30
    max_icwv = 65
    axs.flat[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    axs.flat[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    xb_wv = axs.flat[0].scatter(lon,lat,1.5,c=d_icwv['IC_water'][0,:,:],edgecolors='none', cmap='ocean_r', vmin=min_icwv, vmax=max_icwv,transform=ccrs.PlateCarree())
    # Mark the best track
    if if_btk_exist:
        axs.flat[0].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'white', marker='*',transform=ccrs.PlateCarree())

    # Xa
    axs.flat[1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    axs.flat[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    xa_wv = axs.flat[1].scatter(lon,lat,1.5,c=d_icwv['IC_water'][1,:,:],edgecolors='none', cmap='ocean_r', vmin=min_icwv, vmax=max_icwv,transform=ccrs.PlateCarree())
    # Mark the best track
    if if_btk_exist:
        axs.flat[1].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'white', marker='*',transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.12, 0.1, 0.5, 0.02])
    xwv_bar = fig.colorbar(xb_wv,ax=axs[0:2],orientation="horizontal", cax=caxes)
    xwv_bar.ax.tick_params()

    # Xa-Xb (increment)
    min_incre = -10
    max_incre = 10
    axs.flat[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    axs.flat[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    icwv_incre = d_icwv['IC_water'][1,:,:] - d_icwv['IC_water'][0,:,:]
    incre_wv = axs.flat[2].scatter(lon,lat,1.5,c=icwv_incre,edgecolors='none', cmap='bwr', vmin=min_incre, vmax=max_incre,transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.65, 0.1, 0.25, 0.02])
    cb_diff_ticks = np.linspace(min_incre, max_incre, 5, endpoint=True)
    cbar = fig.colorbar(incre_wv, ax=axs[2:], ticks=cb_diff_ticks, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params()

    #subplot title
    axs.flat[0].set_title('Xb--ICWV ($\mathregular{kgm^{-2}}$)',fontweight='bold')
    axs.flat[1].set_title('Xa--ICWV ($\mathregular{kgm^{-2}}$)',fontweight='bold')
    axs.flat[2].set_title('Xa-Xb--ICWV ($\mathregular{kgm^{-2}}$)', fontweight='bold')    #title for all
    fig.suptitle(Storm+': '+Exper_name+'('+DAtime+')', fontsize=13, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2))
    for j in range(3):
        gl = axs.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
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

    plt.savefig( plot_dir+DAtime+'_ICWV.png', dpi=300 )
    print('Saving the figure: ', plot_dir+DAtime+'_ICWV.png')
    plt.close()

def IC_water( Storm, Exper_name, DAtimes, big_dir, small_dir ):

    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
        print('Reading WRF background and analysis: ', wrf_dir)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/IC_water/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            plot_IC_water( Storm, Exper_name, DAtime, wrf_dir, plot_dir )
        else:
            plot_IC_water( Storm, Exper_name, DAtime, wrf_dir, plot_dir )

# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot relative vorticity per snapshot
# ------------------------------------------------------------------------------------------------------
def read_rt_vo( wrf_dir ):

    P_of_interest = [850,700,]
    # set file names for background and analysis
    mean_dir = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
    # read the shared variables: lat/lon
    ncdir = nc.Dataset( mean_dir[0] )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # calculate relative vorticity and wind
    rtvo_p = np.zeros( [len(mean_dir), len(P_of_interest), np.size(lat,0), np.size(lat,1)] )
    U_p = np.zeros( [len(mean_dir), len(P_of_interest), np.size(lat,0), np.size(lat,1)] )
    V_p = np.zeros( [len(mean_dir), len(P_of_interest), np.size(lat,0), np.size(lat,1)] )
    #windspeed_p = np.zeros( [len(mean_dir), len(P_of_interest), np.size(lat,0), np.size(lat,1)] ) 
    for i in range(len(mean_dir)):
        file_name = mean_dir[i]
        ncdir = nc.Dataset( file_name )
        # pressure
        press = getvar( ncdir, 'pressure')
        # Wind at levels of interest
        Um = getvar( ncdir, 'ua') # U-component of wind on mass points
        Vm = getvar( ncdir, 'va') # V-component of wind on mass points
        # Absolute vorticity
        avo = getvar( ncdir, 'avo') # Absolute vorticity, units: 10-5 s-1
        # Earth vorticity
        Coriolis_sin = ncdir.variables['F'][0,:,:]/1e-5 #units: 10-5 s-1
        # Perform interpolation
        for ip in range(len(P_of_interest)):
            U_p[i,ip,:,:] = interplevel( Um,press,P_of_interest[ip] )
            V_p[i,ip,:,:] = interplevel( Vm,press,P_of_interest[ip] )
            #windspeed_p[i,ip,:,:] = (U_p[i,ip,:,:] ** 2 + V_p[i,ip,:,:] ** 2) ** 0.5
            rtvo_p[i,ip,:,:] = interplevel( avo,press,P_of_interest[ip] )- Coriolis_sin

    d_field = {'lat':lat,'lon':lon,'U':U_p,'V':V_p,'rtvo':rtvo_p,'P':P_of_interest}    
    return d_field

def plot_rt_vo( Storm, Exper_name, DAtime, wrf_dir, plot_dir ):


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

    #------ Read WRFout -------------------
    d_field = read_rt_vo( wrf_dir )
    P_of_interest = d_field['P']
    lat = d_field['lat']
    lon = d_field['lon']
    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(6.8,6.5), dpi=400)
    
    # customize the colormap
    color_intervals = [-5,-2.5,0.0,2.0,4.0,6.0,8.0,10.0]
    exist_cmap = plt.cm.rainbow
    colors = exist_cmap(np.linspace(0,1,len(color_intervals)))
    new_map = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals))

    for i in range(d_field['rtvo'].shape[0]): # files of interest
        for ip in range(d_field['rtvo'].shape[1]): # levels of interest
            ax[i,ip].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
            ax[i,ip].coastlines (resolution='10m', color='black', linewidth=1)
            # relative vorticity (shades)
            rtvo_smooth = d_field['rtvo'][i,ip,:,:]
            rtvo_smooth = sp.ndimage.gaussian_filter( d_field['rtvo'][i,ip,:,:], [2,2])
            min_rtvo = min(color_intervals)
            max_rtvo = max(color_intervals)
            bounds = color_intervals
            rtvo_contourf = ax[i,ip].contourf(lon,lat,rtvo_smooth,cmap=new_map,vmin=min_rtvo,vmax=max_rtvo,levels=bounds,extend='both',transform=ccrs.PlateCarree())
            # wind barbs
            lon_b = lon.flatten()
            lat_b = lat.flatten()
            U_b = d_field['U'][i,ip,:,:].flatten()
            V_b = d_field['V'][i,ip,:,:].flatten()
            ax[i,ip].barbs(lon_b, lat_b, U_b, V_b, length=5, pivot='middle', color='black', regrid_shape=13, transform=ccrs.PlateCarree())
            #ax[i,ip].barbs(lon_b[::20], lat_b[::20], U_b[::5], V_b[::5], length=5, pivot='middle', color='royalblue', regrid_shape=20, transform=ccrs.PlateCarree())
            # Mark the best track
            if if_btk_exist:
                ax[i,ip].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'green', marker='*',transform=ccrs.PlateCarree())

    # Adding the colorbar
    caxes = fig.add_axes([0.2, 0.05, 0.6, 0.02])
    rtvo_bar = fig.colorbar(rtvo_contourf,  orientation="horizontal", cax=caxes)
    rtvo_bar.ax.set_xlabel('Relative Vorticity (10-5 s-1)',fontsize=8)
    rtvo_bar.ax.tick_params(labelsize='8')

    #subplot title
    ax[0,0].set_title( 'Xb: '+str(P_of_interest[0])+' hPa', fontweight='bold')
    ax[0,1].set_title( 'Xb: '+str(P_of_interest[1])+' hPa', fontweight='bold')
    ax[1,0].set_title( 'Xa: '+str(P_of_interest[0])+' hPa', fontweight='bold')
    ax[1,1].set_title( 'Xa: '+str(P_of_interest[1])+' hPa', fontweight='bold')

    #title for all
    fig.suptitle(Storm+': '+Exper_name+'low level circulation', fontsize=8, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(4):
        gl = ax.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top = False
        gl.xlabels_bottom = True
        if j==0 or j==2:
            gl.ylabels_left = True
            gl.ylabels_right = False
        else:
            gl.ylabels_left = False
            gl.ylabels_right = False

        if j==2 or j==3:
            gl.xlabels_bottom = True
            gl.xlabels_top = False
        else:
            gl.xlabels_bottom = False
            gl.xlabels_top = False


        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 8}
        gl.ylabel_style = {'size': 8}

    # Save figures
    des_name = plot_dir + DAtime+'_low_level_cir.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)
    plt.close()


def relative_vo( Storm, Exper_name, DAtimes, big_dir, small_dir ):

    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
        print('Reading WRF background and analysis: ', wrf_dir)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/rt_vo/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            plot_rt_vo( Storm, Exper_name, DAtime, wrf_dir, plot_dir )
        else:
            plot_rt_vo( Storm, Exper_name, DAtime, wrf_dir, plot_dir )


if __name__ == '__main__':

    Storm = 'IRMA'
    Expers = ['IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900',]
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    
    Plot_UV10_slp = True
    Plot_IC_water = True
    Plot_minslp_evo = False
    Plot_rtvo = True


    # Time range set up
    start_time_str = '201708221200'
    end_time_str = '201708221800'
    Consecutive_times = False

    if not Consecutive_times:
        DAtimes = ['201709050000',]
        #DAtimes = ['201708230000','201708230600','201708231200','201708231800','201708240000','201708240600','201708241200']
        #DAtimes = ['201709031200','201709031800','201709040000','201709040600','201709041200','201709041800','201709050000']
        #DAtimes = ['201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Plot low-level circulation
    if Plot_UV10_slp:
        for iExper in Expers:
            start_time=time.process_time()
            UV10_slp( Storm, iExper, DAtimes, big_dir, small_dir )
            end_time = time.process_time()
            print ('time needed: ', end_time-start_time, ' seconds')

    # Plot circulation
    if Plot_rtvo:
        for iExper in Expers:
            start_time=time.process_time()
            relative_vo( Storm, iExper, DAtimes, big_dir, small_dir )
            end_time = time.process_time()
            print ('time needed: ', end_time-start_time, ' seconds')

    # Plot precipitable water (integral column of water vapor)
    if Plot_IC_water:
        for iExper in Expers:
            start_time=time.process_time()
            IC_water( Storm, iExper, DAtimes, big_dir, small_dir )
            end_time = time.process_time()
            print ('time needed: ', end_time-start_time, ' seconds')

    # Plot the evolution of minimum sea level pressure
    if Plot_minslp_evo:
        Evo_slp = Gather_slp( Storm, Expers, DAtimes, big_dir )
        plot_slp_timeseries( small_dir, Storm, Expers, DAtimes, Evo_slp )
    








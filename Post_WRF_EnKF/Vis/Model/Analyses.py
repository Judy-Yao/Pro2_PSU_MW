
import os # functions for interacting with the operating system
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
from wrf import getvar, interplevel
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
import subprocess
import metpy.calc as mpcalc
from metpy.units import units
import numpy.ma as ma

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
        gl.top_labels = False
        gl.bottom_labels = True
        if j==0:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False
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
#           Operation: Read, process, and plot slp per snapshot
# ------------------------------------------------------------------------------------------------------

def slp_plt( Storm, Exper_name, DAtime, wrf_dir, plot_dir ):

    # Read storm center
    dict_btk = UD.read_bestrack(small_dir,Storm)
    # Find the best-track position
    btk_dt = [it_str for it_str in dict_btk['time'] ]#[datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]
    bool_match = [DAtime == it for it in btk_dt]
    if True in bool_match:
        if_btk_exist = True
        idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
    else:
        if_btk_exist = False
    # ------ Read WRFout -------------------
    d_field = read_UV10_slp( DAtime, wrf_dir )
    lat = d_field['lat']
    lon = d_field['lon']
    slp = d_field['slp']
    slp_smooth = d_field['slp_smooth']
    lat_minslp = d_field['lat_minslp']
    lon_minslp = d_field['lon_minslp']
    mslp = d_field['min_slp']

    print( lon_minslp,lat_minslp )

    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )
    min_slp = np.min( slp )
    max_slp = np.max( slp )

    # ------ Plot Figure -------------------
    fig, ax=plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(12,6), dpi=400)

    min_slp = 990
    max_slp = 1015
    bounds = np.linspace(min_slp, max_slp, 6)
    for i in range(2):
        ax[i].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i].coastlines (resolution='10m', color='black', linewidth=1)
        slp_contourf = ax[i].contourf(lon,lat,slp[i,:,:],cmap='inferno',vmin=min_slp,vmax=max_slp,levels=bounds,extend='both',transform=ccrs.PlateCarree())
        # Mark the best track
        if if_btk_exist:
            ax[i].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'green', marker='*',transform=ccrs.PlateCarree())
        
        # Mask the simulated storm center
        ax[i].scatter(lon_minslp[i],lat_minslp[i],20, 'green', marker='s',transform=ccrs.PlateCarree())

    # Adding the colorbar
    caxes = fig.add_axes([0.2, 0.05, 0.6, 0.02])
    slp_bar = fig.colorbar(slp_contourf,  orientation="horizontal", cax=caxes)
    #rtvo_bar.ax.set_xlabel('Sea Level Pressure (hPa)',fontsize=8)
    #rtvo_bar.ax.tick_params(labelsize='10')

    # Title
    ax[0].set_title( 'Xb--min slp: '+str("{0:.3f}".format(mslp[0]))+' hPa',  fontweight='bold') #, fontsize=12)
    ax[1].set_title( 'Xa--min slp: '+str("{0:.3f}".format(mslp[1]))+' hPa',  fontweight='bold') #, fontsize=12)
    fig.suptitle(Storm+': '+Exper_name+'('+DAtime+')', fontsize=12, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2))
    for j in range(2):
        gl = ax[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.bottom_labels = True
        if j==0:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False
        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 12}
        gl.ylabel_style = {'size': 12}

    plt.savefig( plot_dir+DAtime+'_slp.png', dpi=300 )
    print('Saving the figure: ', plot_dir+DAtime+'_slp.png')
    plt.close()

def plot_slp( Storm, Exper_name, DAtimes, big_dir, small_dir ):

    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
        print('Reading WRF background and analysis: ', wrf_dir)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/Slp/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            slp_plt( Storm, Exper_name, DAtime, wrf_dir, plot_dir )
        else:
            slp_plt( Storm, Exper_name, DAtime, wrf_dir, plot_dir )


# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot UV10_slp per snapshot
# ------------------------------------------------------------------------------------------------------
def read_UV10_slp( DAtime, wrf_dir ):

    # set file names for background and analysis
    mean_dir = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
    # read the shared variables: lat/lon
    ncdir = nc.Dataset( mean_dir[0] )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # read UV10 and slp
    slp_original = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    slp_smooth = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    lat_minslp =  []
    lon_minslp = []
    minslp = []

    U10 = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    V10 = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    windspeed = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    for i in range(len(mean_dir)):
        file_name = mean_dir[i]
        ncdir = nc.Dataset( file_name )
        # sea level pressure
        slp = getvar(ncdir, 'slp')
        # original SLP
        slp_values = slp.values
        slp_values[slp_values > 1030] = np.nan
        slp_original[i,:,:] = slp_values
        # smoothed SLP
        slp_smt_values = sp.ndimage.gaussian_filter(slp, [11,11]) #[11,11]
        slp_smt_values[slp_smt_values > 1030] = np.nan
        slp_smooth[i,:,:] = slp_smt_values
        # simulated storm center
        if Storm == 'HARVEY': # Harvey is special with its location near land!

            # eyeball where the storm is
            if DAtime <= '201708221600': 
                lat_mask = lat <= 18
                lon_mask = lon <= -91.5
                mask = lat_mask | lon_mask
            elif (DAtime >= '201708221700') & (DAtime <= '201708222300'):
                lat_mask = lat <= 18
                lon_mask =  (lon >= -88) | (lon <= -91.5)
                mask = lat_mask | lon_mask
            else:
                mask = lat <= 18 
            slp_masked = ma.masked_array(slp_values, mask=mask)
            minslp.append( np.nanmin( slp_masked ) )

            slp_smooth_masked = ma.masked_array(slp_smt_values, mask=mask)
            idx = np.nanargmin( slp_smooth_masked )
            lat_minslp.append( lat.flatten()[idx] )
            lon_minslp.append( lon.flatten()[idx] )
        else:
            minslp.append( np.nanmin( slp_values ) )
            idx = np.nanargmin( slp_smt_values )
            lat_minslp.append( lat.flatten()[idx] )
            lon_minslp.append( lon.flatten()[idx] )

        # UV10
        U10[i,:,:] = ncdir.variables['U10'][0,:,:]
        V10[i,:,:] = ncdir.variables['V10'][0,:,:]
        windspeed[i,:,:] = (U10[i,:,:] ** 2 + V10[i,:,:] ** 2) ** 0.5

    d_UV10_slp = {'lat':lat,'lon':lon,'slp':slp_original,'slp_smooth':slp_smooth,'lat_minslp':lat_minslp,'lon_minslp':lon_minslp,'min_slp':minslp,'U10':U10,'V10':V10,'windspeed':windspeed}
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
    d_field = read_UV10_slp( DAtime, wrf_dir )
    lat = d_field['lat']
    lon = d_field['lon']
    slp = d_field['slp']
    slp_smooth = d_field['slp_smooth']
    lat_minslp = d_field['lat_minslp']
    lon_minslp = d_field['lon_minslp']
    mslp = d_field['min_slp']

    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )
    min_slp = np.nanmin( slp )
    max_slp = np.nanmax( slp )

    # ------ Plot Figure -------------------
    fig, ax=plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(12,6), dpi=400)

    for i in range(2):
        ax[i].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i].coastlines (resolution='10m', color='black', linewidth=1)
        # sea level pressure
        start_level = max( np.nanmin(slp_smooth[0,:,:]), np.nanmin(slp_smooth[1,:,:]) )+0.1
        end_level = min( np.nanmax(slp_smooth[0,:,:]), np.nanmax(slp_smooth[1,:,:]) )
        step = (end_level -start_level)/4
        level = np.arange(start_level, end_level,step)
        #slp_contour = ax[i].contour(lon,lat,slp_smooth[i,:,:],levels=level,cmap='Greys_r',vmin=min_slp,vmax=max_slp,transform=ccrs.PlateCarree())
        slp_contour = ax[i].contour(lon,lat,slp_smooth[i,:,:],levels=level,cmap='winter_r',transform=ccrs.PlateCarree())
        plt.clabel(slp_contour,level,inline=True, fmt="%i", use_clabeltext=True, fontsize=13)
        # Wind at 10 meters
        wind_smooth = sp.ndimage.gaussian_filter( d_field['windspeed'][i,:,:], [2,2])
        min_wind = 2#10
        max_wind = 24#60
        bounds = np.arange(min_wind,max_wind+1,2)
        wind_contourf = ax[i].contourf(lon,lat,wind_smooth,cmap='hot_r',vmin=min_wind,vmax=max_wind,levels=bounds,extend='both',transform=ccrs.PlateCarree())
        ax[i].barbs(lon.flatten(), lat.flatten(), d_field['U10'][i,:,:].flatten(), d_field['V10'][i,:,:].flatten(), length=5, pivot='middle', color='royalblue', regrid_shape=20, transform=ccrs.PlateCarree())

        # Mark the best track
        if if_btk_exist:
            fig.text(0.50,0.05,'Best-Track MSLP:'+str("{0:.3f}".format(dict_btk['min_slp'][idx_btk]))+' hPa', fontsize=12, ha='center', va='center',fontweight='bold')
            ax[i].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 40, 'green', marker='*',transform=ccrs.PlateCarree())

        # Mask the simulated storm center
        ax[i].scatter(lon_minslp[i],lat_minslp[i],40, 'green', marker='s',transform=ccrs.PlateCarree())

    # Adding the colorbar
    cbaxes = fig.add_axes([0.01, 0.1, 0.03, 0.8])
    wind_bar = fig.colorbar(wind_contourf,cax=cbaxes,fraction=0.046, pad=0.04) #Make a colorbar for the ContourSet returned by the contourf call.
    wind_bar.ax.set_ylabel('Wind Speed (m/s)')
    wind_bar.ax.tick_params(labelsize='13')

    # Title
    ax[0].set_title( 'Xb--min slp: '+str("{0:.3f}".format(mslp[0]))+' hPa',  fontweight='bold') #, fontsize=12)
    ax[1].set_title( 'Xa--min slp: '+str("{0:.3f}".format(mslp[1]))+' hPa',  fontweight='bold') #, fontsize=12)
    #ax[0].set_title( 'Xb--min slp: '+str("{0:.3f}".format(np.min( slp[0,:,:] )))+' hPa',  fontweight='bold') #, fontsize=12)
    #ax[1].set_title( 'Xa--min slp: '+str("{0:.3f}".format(np.min( slp[1,:,:] )))+' hPa',  fontweight='bold') #, fontsize=12)
    fig.suptitle(Storm+': '+Exper_name+'('+DAtime+')', fontsize=12, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2))
    for j in range(2):
        gl = ax[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.bottom_labels = True
        if j==0:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False
        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 12}
        gl.ylabel_style = {'size': 12}

    plt.savefig( plot_dir+DAtime+'_UV10_slp', dpi=300 )
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
#           Operation: Read, process, and plot accumulated grid scale precipitation per snapshot
# ------------------------------------------------------------------------------------------------------
def read_Precip( wrf_dir ):

    # set file names for background and analysis
    mean_dir = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
    # read the shared variables: lat/lon
    ncdir = nc.Dataset( mean_dir[0] )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # calculate IC water
    pp = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    for i in range(len(mean_dir)):
        file_name = mean_dir[i]
        ncdir = nc.Dataset( file_name )
        pp[i,:,:] = ncdir.variables['RAINNC'][0,:,:]

    d_pp = {'lat':lat,'lon':lon,'Precip':pp}
    return d_pp

def plot_pp( Storm, Exper_name, DAtime, wrf_dir, plot_dir ):

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
    d_field = read_Precip( wrf_dir )
    lat = d_field['lat']
    lon = d_field['lon']
    pp = d_field['Precip']
    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )
    min_pp = np.min( pp )
    max_pp = np.max( pp )

    # ------ Plot Figure -------------------
    fig, ax=plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(12,6), dpi=400)

    min_pp = 0
    max_pp = 40
    bounds = np.linspace(min_pp, max_pp, 6)
    for i in range(2):
        pp[pp == 0 ] = np.nan # set elements with 0 as np.nan
        ax[i].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i].coastlines(resolution='10m', color='black', linewidth=1)
        pp_contourf = ax[i].contourf(lon,lat,pp[i,:,:],cmap='ocean_r',vmin=min_pp,vmax=max_pp,levels=bounds,extend='max',transform=ccrs.PlateCarree())
        # Mark the best track
        if if_btk_exist:
            ax[i].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'purple', marker='*',transform=ccrs.PlateCarree())

    # Adding the colorbar
    caxes = fig.add_axes([0.2, 0.05, 0.6, 0.02])
    pp_bar = fig.colorbar(pp_contourf,  orientation="horizontal", cax=caxes)
    #rtvo_bar.ax.set_xlabel('Sea Level Pressure (hPa)',fontsize=8)
    #rtvo_bar.ax.tick_params(labelsize='10')

    # Title
    #ax[0].set_title( 'Xb--min slp: '+str("{0:.3f}".format(np.nanmin( slp[0,:,:] )))+' hPa',  fontweight='bold') #, fontsize=12)
    #ax[1].set_title( 'Xa--min slp: '+str("{0:.3f}".format(np.nanmin( slp[1,:,:] )))+' hPa',  fontweight='bold') #, fontsize=12)
    fig.suptitle(Storm+': '+Exper_name+'('+DAtime+')', fontsize=12, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2))
    for j in range(2):
        gl = ax[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.bottom_labels = True
        if j==0:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False
        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 12}
        gl.ylabel_style = {'size': 12}

    plt.savefig( plot_dir+DAtime+'_Precip', dpi=300 )
    print('Saving the figure: ', plot_dir+DAtime+'_Precip.png')
    plt.close()


def Precip( Storm, Exper_name, DAtimes, big_dir, small_dir ):

    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
        print('Reading WRF background and analysis: ', wrf_dir)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/Precip/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            plot_pp( Storm, Exper_name, DAtime, wrf_dir, plot_dir )
        else:
            plot_pp( Storm, Exper_name, DAtime, wrf_dir, plot_dir )


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
    dict_btk = UD.read_bestrack(small_dir, Storm)
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
        gl.top_labels = False
        gl.bottom_labels = True
        if j==0:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False
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
    color_intervals = [-5,-2.5,0.0,5.0,7.5,10.0,12.5,15.0,]
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

        gl.top_labels = False
        gl.bottom_labels = True
        if j==0 or j==2:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False

        if j==2 or j==3:
            gl.bottom_labels = True
            gl.top_labels = False
        else:
            gl.bottom_labels = False
            gl.top_labels = False


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


# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot divergence per snapshot
# ------------------------------------------------------------------------------------------------------
def read_diver( wrf_dir ):

    P_of_interest = [850,700,]
    # set file names for background and analysis
    mean_dir = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
    # read the shared variables: lat/lon
    ncdir = nc.Dataset( mean_dir[0] )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # calculate relative vorticity and wind
    dvg_p = np.zeros( [len(mean_dir), len(P_of_interest), np.size(lat,0), np.size(lat,1)] )
    U_p = np.zeros( [len(mean_dir), len(P_of_interest), np.size(lat,0), np.size(lat,1)] )
    V_p = np.zeros( [len(mean_dir), len(P_of_interest), np.size(lat,0), np.size(lat,1)] )
    for i in range(len(mean_dir)):
        file_name = mean_dir[i]
        ncdir = nc.Dataset( file_name )
        # pressure
        press = getvar( ncdir, 'pressure')
        # Wind at levels of interest
        Um = getvar( ncdir, 'ua') # U-component of wind on mass points
        Vm = getvar( ncdir, 'va') # V-component of wind on mass points
        Um = Um* units("m/s")
        Vm = Vm* units("m/s")
        # divergence
        dvg = mpcalc.divergence(Um, Vm, dx=3000* units("m"), dy=3000* units("m"))
        # Perform interpolation
        for ip in range(len(P_of_interest)):
            U_p[i,ip,:,:] = interplevel( Um,press,P_of_interest[ip] )
            V_p[i,ip,:,:] = interplevel( Vm,press,P_of_interest[ip] )
            dvg_p[i,ip,:,:] = interplevel( dvg,press,P_of_interest[ip] )
    d_field = {'lat':lat,'lon':lon,'U':U_p,'V':V_p,'dvg':dvg_p*1e5,'P':P_of_interest}
    return d_field

def plot_diver( Storm, Exper_name, DAtime, wrf_dir, plot_dir ):

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
    d_field = read_diver( wrf_dir )
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
    color_intervals = [-10,-5,0.0,5,10]
    #exist_cmap = plt.cm.rainbow
    #colors = exist_cmap(np.linspace(0,1,len(color_intervals)))
    #new_map = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals))

    for i in range(d_field['dvg'].shape[0]): # files of interest
        for ip in range(d_field['dvg'].shape[1]): # levels of interest
            ax[i,ip].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
            ax[i,ip].coastlines (resolution='10m', color='black', linewidth=1)
            # divergence (shades)
            #dvg_smooth = d_field['dvg'][i,ip,:,:]
            dvg_smooth = sp.ndimage.gaussian_filter( d_field['dvg'][i,ip,:,:], [2,2])
            min_dvg = min(color_intervals)
            max_dvg = max(color_intervals)
            bounds = color_intervals
            dvg_contourf = ax[i,ip].contourf(lon,lat,dvg_smooth,cmap='bwr',vmin=min_dvg,vmax=max_dvg,levels=bounds,extend='both',transform=ccrs.PlateCarree())
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
    dvg_bar = fig.colorbar(dvg_contourf,  orientation="horizontal", cax=caxes)
    dvg_bar.ax.set_xlabel('Divergence (10-5 s-1)',fontsize=8)
    dvg_bar.ax.tick_params(labelsize='8')

    #subplot title
    ax[0,0].set_title( 'Xb: '+str(P_of_interest[0])+' hPa', fontweight='bold')
    ax[0,1].set_title( 'Xb: '+str(P_of_interest[1])+' hPa', fontweight='bold')
    ax[1,0].set_title( 'Xa: '+str(P_of_interest[0])+' hPa', fontweight='bold')
    ax[1,1].set_title( 'Xa: '+str(P_of_interest[1])+' hPa', fontweight='bold')

    #title for all
    fig.suptitle(Storm+': '+Exper_name+'low level divergence', fontsize=8, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(4):
        gl = ax.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')

        gl.top_labels = False
        gl.bottom_labels = True
        if j==0 or j==2:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False

        if j==2 or j==3:
            gl.bottom_labels = True
            gl.top_labels = False
        else:
            gl.bottom_labels = False
            gl.top_labels = False


        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 8}
        gl.ylabel_style = {'size': 8}

    # Save figures
    des_name = plot_dir + DAtime+'_divergence.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)
    plt.close()


def diver( Storm, Exper_name, DAtimes, big_dir, small_dir ):

    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
        print('Reading WRF background and analysis: ', wrf_dir)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/divergence/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
        plot_diver( Storm, Exper_name, DAtime, wrf_dir, plot_dir )


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # -------- Configuration -----------------
    Storm = 'IRMA'
    DA = ['IR']   
    MP = 'THO' 

    Plot_Precip = False       # Accumulated total grid scale precipitation
    Plot_slp = False
    Plot_UV10_slp = False
    Plot_IC_water = True
    Plot_minslp_evo = False
    Plot_rtvo = False
    Plot_divergence = False
    # -----------------------------------------

    # Time range set up
    start_time_str = '201709030000'
    end_time_str = '201709030600'
    Consecutive_times = True

    if not Consecutive_times:
        DAtimes = ['201708221200','201708221800','201708230000','201708230600','201708231200',]
        #DAtimes = ['201709031200','201709031800','201709040000','201709040600','201709041200','201709041800','201709050000']
        #DAtimes = ['201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    ## Experiment name
    Expers= []
    for ida in DA:
        Expers.append( UD.generate_one_name( Storm,ida,MP ) )
    #Expers= ['IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900',]

    # Plot slp
    if Plot_slp:
        for iExper in Expers:
            start_time=time.process_time()
            plot_slp( Storm, iExper, DAtimes, big_dir, small_dir )
            end_time = time.process_time()
            print ('time needed: ', end_time-start_time, ' seconds')

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

    # Plot divergence
    if Plot_divergence:
        for iExper in Expers:
            start_time=time.process_time()
            diver( Storm, iExper, DAtimes, big_dir, small_dir )
            end_time = time.process_time()
            print ('time needed: ', end_time-start_time, ' seconds')


    # Plot precipitable water (integral column of water vapor)
    if Plot_IC_water:
        for iExper in Expers:
            start_time=time.process_time()
            IC_water( Storm, iExper, DAtimes, big_dir, small_dir )
            end_time = time.process_time()
            print ('time needed: ', end_time-start_time, ' seconds')

    # Plot accumulated total grid scale precipitation
    if Plot_Precip:
        for iExper in Expers:
            start_time=time.process_time()
            Precip( Storm, iExper, DAtimes, big_dir, small_dir )
            end_time = time.process_time()
            print ('time needed: ', end_time-start_time, ' seconds')

    # Plot the evolution of minimum sea level pressure
    if Plot_minslp_evo:
        Evo_slp = Gather_slp( Storm, Expers, DAtimes, big_dir )
        plot_slp_timeseries( small_dir, Storm, Expers, DAtimes, Evo_slp )
    








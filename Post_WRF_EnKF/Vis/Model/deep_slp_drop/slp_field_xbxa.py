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
import matplotlib.patches as patches
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
from wrf import getvar 
import numpy.ma as ma

from Track_xbxa import read_HPI_model
import Util_data as UD

def generate_times( start_time_str, end_time_str, interval ):

    time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
    time_diff_hour = time_diff.total_seconds() / 3600
    time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+interval, interval))]
    d_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return d_times


# Interpolate from 6-hour value to every hour value
def Interpolate_hourly( l_6hrs,d_obs6Hr,l_hrs,key):

    # group pairs of synoptic times 
    syno_ts = l_6hrs
    syno_pair = [[] for i in range(2)]
    for i in range(len(syno_ts)-1):
        syno_pair[0].append( syno_ts[i] )
        syno_pair[1].append( syno_ts[i+1] )

    # Linearly interpolate the Vmax between two synoptic times
    d_ObsHr = []

    for tt in l_hrs: #target time

        if any( hh in tt[8:10] for hh in ['00','06','12','18']):
            d_ObsHr.append( d_obs6Hr[tt][key]  )
        else:
            lg = [ syno_pair[0][ip] < tt < syno_pair[1][ip] for ip in range(len(syno_pair[0]))]
            idx = int( np.where(lg)[0] )

            # Calculate time differences
            time_diff_total = (datetime.strptime(syno_pair[1][idx],"%Y%m%d%H%M") - datetime.strptime(syno_pair[0][idx],"%Y%m%d%H%M")).total_seconds()
            time_diff_target = (datetime.strptime(tt,"%Y%m%d%H%M") - datetime.strptime(syno_pair[0][idx],"%Y%m%d%H%M")).total_seconds()

            # Calculate interpolation factor
            interp_factor = time_diff_target / time_diff_total

            # Perform linear interpolation for latitude and longitude
            d_ObsHr.append( d_obs6Hr[syno_pair[0][idx]][key] + interp_factor * (d_obs6Hr[syno_pair[1][idx]][key] - d_obs6Hr[syno_pair[0][idx]][key]) )

    return d_ObsHr




def read_slp( wrf_dir,DAtime ):

    # set file names for background and analysis
    mean_dir = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
    # read the shared variables: lat/lon
    ncdir = nc.Dataset( mean_dir[0] )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # read UV10 and slp
    slp_original = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    slp_smooth = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
    minslp  = np.zeros( [len(mean_dir), ] )
    lat_minslp = np.zeros( [len(mean_dir), ] )
    lon_minslp = np.zeros( [len(mean_dir), ] )

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
            slp_masked = ma.masked_array(slp_original[i,:,:], mask=mask)
            minslp[i] = np.nanmin( slp_masked )

            slp_smooth_masked = ma.masked_array(slp_smooth[i,:,:], mask=mask)
            idx = np.nanargmin( slp_smooth_masked )
            lat_minslp[i] =  lat.flatten()[idx]
            lon_minslp[i] = lon.flatten()[idx]
        else:
            minslp[i] = np.nanmin( slp_original[i,:,:] )
            idx = np.nanargmin( slp_smooth[i,:,:] )
            lat_minslp[i] =  lat.flatten()[idx]
            lon_minslp[i] = lon.flatten()[idx]


    d_slp = {'lat':lat,'lon':lon,'slp':slp_original,'slp_smooth':slp_smooth,'lat_minslp':lat_minslp,'lon_minslp':lon_minslp,'minslp':minslp}
    return d_slp



def plot_slp( Storm, Exper_name, DAtime, wrf_dir, plot_dir ):

    # Read location from TCvitals
    #if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
    #    tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir,Storm, DAtime)

    # ------ Read WRFout -------------------
    d_field = read_slp( wrf_dir, DAtime)
    lat = d_field['lat']
    lon = d_field['lon']
    slp = d_field['slp']
    slp_smooth = d_field['slp_smooth']
    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )

    # ------ Plot Figure -------------------
    fig, ax=plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(12,6), dpi=400)

    for i in range(2):
        ax[i].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i].coastlines (resolution='10m', color='black', linewidth=1)
        # sea level pressure
        min_slp = 940 #1000
        max_slp = 1000 #1015 
        bounds = np.linspace(min_slp, max_slp, 7)
        slp_contourf = ax[i].contourf(lon,lat,slp[i,:,:],cmap='magma',vmin=min_slp,vmax=max_slp,levels=bounds,transform=ccrs.PlateCarree())
        # location of min slp
        lon_minslp = d_field['lon_minslp']
        lat_minslp = d_field['lat_minslp']
        ax[i].scatter(lon_minslp,lat_minslp, 20, 'black', alpha=1,marker='o',transform=ccrs.PlateCarree())
        # Mark the best track
        #if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        it = DAtimes.index( DAtime )
        ax[i].scatter(d_tcvHr['lon'][it],d_tcvHr['lat'][it], 20, 'black', marker='*',alpha=0.6, transform=ccrs.PlateCarree())
        ax[i].text(d_tcvHr['lon'][it]-1,d_tcvHr['lat'][it]-0.5,round(d_tcvHr['mslp'][it],2),fontsize=12)
        #ax[i].add_patch(patches.Polygon(path,facecolor='none',edgecolor='black',linewidth=1.5 ))

    # Adding the colorbar
    cbaxes = fig.add_axes([0.01, 0.1, 0.02, 0.8])
    wind_bar = fig.colorbar(slp_contourf,cax=cbaxes,fraction=0.046, pad=0.04) #Make a colorbar for the ContourSet returned by the contourf call.
    #wind_bar.set_clim( vmin=min_slp, vmax=max_slp )
    wind_bar.ax.set_ylabel('slp (hPa)',fontsize=12)
    wind_bar.ax.tick_params(labelsize='12')

    # Title
    ax[0].set_title( 'Xb--min slp: '+str("{0:.3f}".format(d_field['minslp'][0]))+' hPa',  fontweight='bold', fontsize=13)
    ax[1].set_title( 'Xa--min slp: '+str("{0:.3f}".format(d_field['minslp'][1]))+' hPa',  fontweight='bold', fontsize=13)
    if deep_slp_incre:
        title_name = Storm+': '+Exper_name+' ('+DAtime+')'+'\nAbs of min SLP increment > '+str(incre_slp_th)+' hPa'
        title_name = title_name #+ ' kick5error'
        fig.suptitle(title_name,fontsize=12, fontweight='bold')
    else:
        fig.suptitle(Storm+': '+Exper_name+' ('+DAtime+')', fontsize=12, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for j in range(2):
        gl = ax[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
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

    des = plot_dir+DAtime+'_slp_hroi_100km.png'
    plt.savefig( des, dpi=300 )
    print('Saving the figure: ', des)
    plt.close()


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'CONV'
    MP = 'WSM6'

    # Time range set up
    start_time_str = '201709031800'
    end_time_str = '201709031800'
    Consecutive_times = True

    deep_slp_incre = True
    incre_slp_th = 0 # threshold of increment, unit:hpa 
    Plot_slp = True
    # -------------------------------------------------------   

    # Create experiment names
    Exper_name = UD.generate_one_name( Storm,DA,MP )

    if not Consecutive_times:
        #DAtimes = ['201709050000','201709050600','201709051200','201709051800','201709060000']
        #DAtimes = ['201708230000','201708230600','201708231200','201708231800','201708240000','201708240600','201708241200']
        #DAtimes = ['201709031200','201709031800','201709040000','201709040600','201709041200','201709041800','201709050000']
        DAtimes = ['201709160000','201709160600','201709161200','201709161800','201709170000',]
    else:
        DAtimes = generate_times( start_time_str, end_time_str, 1)
    l_6hrs = generate_times( start_time_str, end_time_str, 6)

    # Read TCvitals
    d_tcv6Hrs = {}
    for it in l_6hrs:
        d_tcv6Hrs[it] = {}
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir,Storm,it)
        d_tcv6Hrs[it]['lon'] = tc_lon
        d_tcv6Hrs[it]['lat'] = tc_lat
        d_tcv6Hrs[it]['mslp'] = tc_slp
    # interpolate to each hour
    d_tcvHr = {}
    d_tcvHr['time'] = DAtimes
    d_tcvHr['lon'] = Interpolate_hourly( l_6hrs,d_tcv6Hrs,DAtimes,'lon')
    d_tcvHr['lat'] = Interpolate_hourly( l_6hrs,d_tcv6Hrs,DAtimes,'lat')
    d_tcvHr['mslp'] = Interpolate_hourly( l_6hrs,d_tcv6Hrs,DAtimes,'mslp')

    # Plot low-level circulation
    if Plot_slp and deep_slp_incre:
        start_time=time.process_time()
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/deep_slp_incre/slp_field/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
        # ----- Read min slp from model-----------------
        HPI_models = {}
        DAtimes_dir = [big_dir+Storm+'/'+Exper_name+'/fc/'+it for it in DAtimes]
        file_kinds = ['wrf_enkf_input_d03_mean','wrf_enkf_output_d03_mean']
        for ifk in file_kinds:
            idx = file_kinds.index( ifk )
            HPI_models[ifk] = read_HPI_model( Storm, Exper_name, ifk, DAtimes_dir )
        incre_slp = np.array(HPI_models['wrf_enkf_output_d03_mean']['min_slp']) - np.array(HPI_models['wrf_enkf_input_d03_mean']['min_slp'])
        # Loop through each DAtime/analysis
        for DAtime in DAtimes:
            idx_t = DAtimes.index( DAtime )
            if abs(incre_slp[idx_t]) >= incre_slp_th:
                wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
                print('Reading WRF background and analysis: ', wrf_dir)
                DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
                plot_slp( Storm, Exper_name, DAtime, wrf_dir, plot_dir )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')

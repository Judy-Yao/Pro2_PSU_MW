#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
from wrf import getvar
# It might be possible that you are not able to conda install wrf-var with a pretty new python version
# Solution:
# 1. conda create -n $PYTHON34_ENV_NAME python=3.4 anaconda 
# 2. conda activate python=3.4 (use wrf-python in this python environment)
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
from Util_data import read_bestrack
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import subprocess
from itertools import chain

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
#           Operation: Read, process, and plot precipitable water
# ------------------------------------------------------------------------------------------------------
def read_IC_water( Expers,wrf_dirs, DAtime ):

    d_ICw_all = {}
    for iExper in Expers:
        idx_exper = Expers.index( iExper )
        if not meanOverEns:
            # set file names for background and analysis
            mean_dir = [wrf_dirs[idx_exper]+DAtime+'/wrf_enkf_input_d03_mean',] #wrf_dirs[idx_exper]+DAtime+'/wrf_enkf_output_d03_mean']
            # read the shared variables: lat/lon
            ncdir = nc.Dataset( mean_dir[0] )
            lat = ncdir.variables['XLAT'][0,:,:]
            lon = ncdir.variables['XLONG'][0,:,:]
            # calculate IC water
            IC_water = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
            for i in range(len(mean_dir)):
                file_name = mean_dir[i]
                ncdir = nc.Dataset( file_name )
                IC_water[i,:,:] = getvar( ncdir,'pw' )
            d_IC_water = {'lat':lat,'lon':lon,'IC_water':IC_water}
        else:
            mem_dirs = sorted( glob.glob(wrf_dirs[idx_exper]+DAtime+'/wrf_enkf_input_d03_0*') )
            # read the shared variables: lat/lon
            ncdir = nc.Dataset( mem_dirs[0] )
            lat = ncdir.variables['XLAT'][0,:,:]
            lon = ncdir.variables['XLONG'][0,:,:]
            IC_water = np.zeros( [len(mem_dirs), np.size(lat,0), np.size(lat,1)] )
            for mem_dir in mem_dirs:
                print('Reading '+mem_dir)
                file_name = mem_dir
                idx = mem_dirs.index( mem_dir )
                ncdir = nc.Dataset( file_name )
                IC_water[idx,:,:] = getvar( ncdir,'pw' )
           
            meanOverEns_ICw = np.mean(IC_water,axis=0)
            print(np.shape(meanOverEns_ICw))
            d_IC_water = {'lat':lat,'lon':lon,'IC_water':meanOverEns_ICw}

        d_ICw_all[iExper] = d_IC_water

    return d_ICw_all


def plot_IC_water( Storm, Expers, DAtime, wrf_dirs, plot_dir ):

    # Read storm center
    dict_btk = read_bestrack(Storm)
    # Find the best-track position
    btk_dt = [it_str for it_str in dict_btk['time'] ]#[datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]
    bool_match = [DAtime == it for it in btk_dt]
    if True in bool_match:
        if_btk_exist = True
        idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
    else:
        if_btk_exist = False

    #------ Read WRFout -------------------
    d_icwv = read_IC_water( Expers, wrf_dirs, DAtime )


    # Big Assumption: in the same domain!!!
    lat = d_icwv[Expers[0]]['lat']
    lon = d_icwv[Expers[0]]['lon']
    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )

    # ------------------ Plot -----------------------
    fig, axs=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(12,5), dpi=400)

    # Exper 1: Xb
    min_icwv = 0
    max_icwv = 75
    axs.flat[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    axs.flat[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    xb_wv = axs.flat[0].scatter(lon,lat,1.5,c=d_icwv[Expers[0]]['IC_water'][:,:],edgecolors='none', cmap='ocean_r', vmin=min_icwv, vmax=max_icwv,transform=ccrs.PlateCarree())
    # Mark the best track
    if if_btk_exist:
        axs.flat[0].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'white', marker='*',transform=ccrs.PlateCarree())

    # Exper 2: Xb
    axs.flat[1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    axs.flat[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    xa_wv = axs.flat[1].scatter(lon,lat,1.5,c=d_icwv[Expers[1]]['IC_water'][:,:],edgecolors='none', cmap='ocean_r', vmin=min_icwv, vmax=max_icwv,transform=ccrs.PlateCarree())
    # Mark the best track
    if if_btk_exist:
        axs.flat[1].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'white', marker='*',transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.12, 0.1, 0.5, 0.02])
    xwv_bar = fig.colorbar(xb_wv,ax=axs[0:2],orientation="horizontal", cax=caxes)
    xwv_bar.ax.tick_params()

    # Exper 1: Xb - Exper 2: Xb
    min_incre = -5
    max_incre = 5
    axs.flat[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    axs.flat[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    icwv_incre = d_icwv[Expers[0]]['IC_water'][:,:] - d_icwv[Expers[1]]['IC_water'][:,:]
    incre_wv = axs.flat[2].scatter(lon,lat,1.5,c=icwv_incre,edgecolors='none', cmap='bwr', vmin=min_incre, vmax=max_incre,transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.65, 0.1, 0.25, 0.02])
    cb_diff_ticks = np.linspace(min_incre, max_incre, 5, endpoint=True)
    cbar = fig.colorbar(incre_wv, ax=axs[2:], ticks=cb_diff_ticks, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params()

    #subplot title
    axs.flat[0].set_title('Xb1--ICWV ($\mathregular{kgm^{-2}}$)',fontweight='bold')
    axs.flat[1].set_title('Xb2--ICWV ($\mathregular{kgm^{-2}}$)',fontweight='bold')
    axs.flat[2].set_title('Xb1-Xb2--ICWV ($\mathregular{kgm^{-2}}$)', fontweight='bold')

    #title for all
    fig.suptitle(Storm+': '+'('+DAtime+')', fontsize=13, fontweight='bold')

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

    plt.savefig( plot_dir+DAtime+'_compare_ICWV.png', dpi=300 )
    print('Saving the figure: ', plot_dir+DAtime+'_compare_ICWV.png')
    plt.close()

def IC_water( Storm, Expers, DAtimes, big_dir, small_dir ):

    wrf_dirs = []
    for iExper in Expers:
        wrf_dirs.append( big_dir+Storm+'/'+iExper+'/fc/' )
    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        print('Reading WRF background and analysis at ', DAtime)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Expers[0]+'/Vis_analyze/Model/IC_water/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            plot_IC_water( Storm, Expers, DAtime, wrf_dirs, plot_dir )
        else:
            plot_IC_water( Storm, Expers, DAtime, wrf_dirs, plot_dir )


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
    ax.plot_date(dates, Evo_slp['obs_slp'], 'black', label='TCvital', linestyle='-',linewidth='4',)
    if slp_xa and not slp_xb: 
        ax.plot_date(dates, Evo_slp['xa_slp'][0,:], 'red', label='YES', linestyle='-',linewidth='4',)
        #ax.plot_date(dates, Evo_slp['xa_slp'][1,:], 'red', label='IR+MW-THO', linestyle='--',linewidth='4',)
        ax.plot_date(dates, Evo_slp['xa_slp'][1,:], 'blue', label='NO', linestyle='-',linewidth='4',)
        #ax.plot_date(dates, Evo_slp['xa_slp'][3,:], 'blue', label='IR+MW-WSM6', linestyle='--',linewidth='4',)
    if slp_xb and not slp_xa:
        ax.plot_date(dates, Evo_slp['xb_slp'][0,:], 'red', label='YES', linestyle='-',linewidth='4',)
        ax.plot_date(dates, Evo_slp['xb_slp'][1,:], 'blue', label='NO', linestyle='-',linewidth='4',)
    if slp_xa and slp_xb: # plot saw-tooth lines
        dates_zip = list( chain.from_iterable( zip(dates,dates)) )
        len_seg = len(dates_zip)-1
        colors = ['red','blue']
        labels = ['YES','NO']
        for iexper in range(len(Expers)):
            slp_zip = list( chain.from_iterable( zip(Evo_slp['xb_slp'][iexper,:],Evo_slp['xa_slp'][iexper,:]) ) )
            for i in range(1,len_seg):
                # specify which segment uses which line style
                if i % 2 == 0:
                    style = '-'
                else:
                    style = '--'  
                
                if i == 2:
                    ax.plot(dates_zip[i-1:i+1],slp_zip[i-1:i+1],linestyle=style,color=colors[iexper],linewidth='4',label=labels[iexper]) 
                else:
                    ax.plot(dates_zip[i-1:i+1],slp_zip[i-1:i+1],linestyle=style,color=colors[iexper],linewidth='4')


    leg = plt.legend(loc='lower left',fontsize=22)
    # Set X/Y labels
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45,labelsize=15)
    ax.set_ylabel('Minimum Sea Level Pressure (hPa)',fontsize=24)
    ax.set_ylim(900,990) #( 920,1020 ) #(970,1020)

    #title_name = Storm+'('+Exper_name+')'+': mim SLP'
    title_name = Storm+': mim SLP (hPa)'
    ax.set_title( title_name,fontweight="bold",fontsize='24' )
    #fig.suptitle('WSM6', fontsize=24, fontweight='bold')

    # Save the figure
    save_des = small_dir+Storm+'/'+Expers[0]+'/Vis_analyze/Model/minslp_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


def Gather_slp( Storm, Expers, DAtimes, big_dir ):

    # Initialize container
    obs_minslp_lat = []
    obs_minslp_lon = []
    obs_minslp_value = []
   
    if slp_xa:
        xa_minslp_lat = [[] for i in range(len(Expers))]
        xa_minslp_lon = [[] for i in range(len(Expers))]
        xa_minslp_value = [[] for i in range(len(Expers))]

    if slp_xb:
        xb_minslp_lat = [[] for i in range(len(Expers))]
        xb_minslp_lon = [[] for i in range(len(Expers))]
        xb_minslp_value = [[] for i in range(len(Expers))]

    # Obtain min slp
    for DAtime in DAtimes:

        # collect min slp found from WRF output
        for Exper in Expers:
            idx = Expers.index( Exper )
            if slp_xa:
                xa_dir = big_dir+Storm+'/'+Exper+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
                print('Reading the EnKF posterior mean from ', xa_dir)
                list_xa = find_minSLP( xa_dir )
                xa_minslp_lat[idx].append( list_xa[0] )
                xa_minslp_lon[idx].append( list_xa[1] )
                xa_minslp_value[idx].append( list_xa[2] )

            if slp_xb:
                xb_dir = big_dir+Storm+'/'+Exper+'/fc/'+DAtime+'/wrf_enkf_input_d03_mean'
                print('Reading the EnKF prior mean from ', xb_dir)
                list_xb = find_minSLP( xb_dir )
                xb_minslp_lat[idx].append( list_xb[0] )
                xb_minslp_lon[idx].append( list_xb[1] )
                xb_minslp_value[idx].append( list_xb[2] )

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

    # Assemble the dictionary
    dict_minSLP = {'obs_slp':np.array(obs_minslp_value),}
    if slp_xa:
        dict_minSLP['xa_slp'] = np.array(xa_minslp_value)
        dict_minSLP['xa_lat'] = np.array(xa_minslp_lat)
        dict_minSLP['xa_lon'] = np.array(xa_minslp_lon)
    if slp_xb:
        dict_minSLP['xb_slp'] = np.array(xb_minslp_value)
        dict_minSLP['xb_lat'] = np.array(xb_minslp_lat)
        dict_minSLP['xb_lon'] = np.array(xb_minslp_lon) 

    return dict_minSLP



if __name__ == '__main__':

    Storm = 'IRMA'
    Expers = ['IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900','IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900']
    #Expers = ['IR-J_DA+J_WRF+J_init-SP-intel17-THO-30hr-hroi900','IR+MW-J_DA+J_WRF+J_init-SP-intel17-THO-30hr-hroi900','IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900','IR+MW-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900',]
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    meanOverEns = False  # mean of H(ens) or H of mean(ens)
    Com_IC_water = False
    Com_minslp_evo = True

    slp_xa = True
    slp_xb = True

    # Time range set up
    start_time_str = '201709030000'
    end_time_str = '201709050000'
    Consecutive_times = True

    if not Consecutive_times:
        DAtimes = ['201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Plot precipitable water (integral column of water vapor)
    if Com_IC_water:
        start_time=time.process_time()
        IC_water( Storm, Expers, DAtimes, big_dir, small_dir )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')

    # Plot the evolution of minimum sea level pressure
    if Com_minslp_evo:
        Evo_slp = Gather_slp( Storm, Expers, DAtimes, big_dir )
        plot_slp_timeseries( small_dir, Storm, Expers, DAtimes, Evo_slp )

 


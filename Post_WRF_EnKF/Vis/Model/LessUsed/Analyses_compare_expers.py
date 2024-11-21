#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
from wrf import getvar
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
import subprocess
from itertools import chain
import numpy.ma as ma

from Util_data import read_bestrack
import Util_data as UD

# setting font size to 15
plt.rcParams.update({'font.size': 15})

def RMSE(simu, obs):
    return np.sqrt( ((simu - obs) ** 2).mean() )

def Bias(simu, obs):
    return  np.sum((simu - obs),0)/np.size(obs,0)

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

    min_slp = [lat_minslp,lon_minslp,min_slp.values]
    return min_slp

def plot_slp_timeseries( small_dir, Storm, Expers, DAtimes, Evo_slp ):

    # Set up figure
    fig = plt.figure( figsize=(10,9), dpi=300 )
    ax = plt.subplot(1,1,1)
    dates = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]
    # Plot obs and mean
    ax.plot_date(dates, Evo_slp['obs_slp'], 'black', label='TCvitals', linestyle='-',linewidth='4',)
    if slp_xa and not slp_xb: 
        ax.plot_date(dates, Evo_slp['xa_slp'][0,:], 'red', label='IR-THO', linestyle='-',linewidth='4',)
        ax.plot_date(dates, Evo_slp['xa_slp'][1,:], 'red', label='IR+MW-THO', linestyle='--',linewidth='4',)
        ax.plot_date(dates, Evo_slp['xa_slp'][2,:], 'blue', label='IR-WSM6', linestyle='-',linewidth='4',)
        ax.plot_date(dates, Evo_slp['xa_slp'][3,:], 'blue', label='IR+MW-WSM6', linestyle='--',linewidth='4',)
    if slp_xb and not slp_xa:
        ax.plot_date(dates, Evo_slp['xb_slp'][0,:], 'red', label='YES', linestyle='-',linewidth='4',)
        ax.plot_date(dates, Evo_slp['xb_slp'][1,:], 'blue', label='NO', linestyle='-',linewidth='4',)
    if slp_xa and slp_xb: # plot saw-tooth lines
        dates_zip = list( chain.from_iterable( zip(dates,dates)) )
        len_seg = len(dates_zip)
        # Customize labels
        labels = {}
        for imp in MP:
            if imp == 'THO': #'TuneWSM6':
                labels[imp] = {'xa':'THO_DA','xb':'THO_FC'} #{'xa':'MD_DA','xb':'MD_FC'}
            else:
                labels[imp] = {'xa':'WSM6_DA','xb':'WSM6_FC'} #{'xa':'Ref_DA','xb':'Ref_FC'}

        bias = {}
        rmse = {} 
        # Loop thru experiments
        for iexper in range(len(Expers)):
            # Calculate T-mean RMSE
            rmse[MP[iexper]] = {'xb':RMSE(Evo_slp['xb_slp'][iexper,:], Evo_slp['obs_slp'] ),'xa':RMSE(Evo_slp['xa_slp'][iexper,:], Evo_slp['obs_slp'] )}
            bias[MP[iexper]] = {'xb':Bias(Evo_slp['xb_slp'][iexper,:], Evo_slp['obs_slp'] ),'xa':Bias(Evo_slp['xa_slp'][iexper,:], Evo_slp['obs_slp'] )}

            slp_zip = list( chain.from_iterable( zip(Evo_slp['xb_slp'][iexper,:],Evo_slp['xa_slp'][iexper,:]) ) )
            for i in range(1,len_seg):
                # specify which segment uses which line style
                if i % 2 == 0:
                    line = '-'
                else:
                    line = '--'  

              # customize the line
                if MP[iexper] == 'THO':
                    color = 'blue'
                    if i == 1:
                        ax.plot(dates_zip[i-1:i+1],slp_zip[i-1:i+1],color,linestyle=line,linewidth='4',label=labels[MP[iexper]]['xa'])
                    elif i == 2:
                        ax.plot(dates_zip[i-1:i+1],slp_zip[i-1:i+1],color,linestyle=line,linewidth='4',label=labels[MP[iexper]]['xb']) #WSM6_Forecast
                    else:
                        ax.plot(dates_zip[i-1:i+1],slp_zip[i-1:i+1],color,linestyle=line,linewidth='4')
                else:
                    color = 'red'
                    if i == 1:
                        ax.plot(dates_zip[i-1:i+1],slp_zip[i-1:i+1],color,linestyle=line,linewidth='6',label=labels[MP[iexper]]['xa'])
                    elif i == 2:
                        ax.plot(dates_zip[i-1:i+1],slp_zip[i-1:i+1],color,linestyle=line,linewidth='4',label=labels[MP[iexper]]['xb'])
                    else:
                        ax.plot(dates_zip[i-1:i+1],slp_zip[i-1:i+1],color,linestyle=line,linewidth='4')
        # title
        bias_str = '\nBias: '
        rmse_str = '\nRMSE: '
        suptt = Storm+'; Over '+str(len(DAtimes))+' Cycles; Min SLP (hPa)'
        for imp in MP:
            bias_str= bias_str+labels[imp]['xb']+' '+'%.1f' %bias[imp]['xb']+'; ' 
            bias_str= bias_str+labels[imp]['xa']+' '+'%.1f' %bias[imp]['xa']+'; '
            rmse_str = rmse_str+labels[imp]['xb']+' '+'%.1f' %rmse[imp]['xb']+'; '  
            rmse_str= rmse_str+labels[imp]['xa']+' '+'%.1f' %rmse[imp]['xa']+'; '
        supt = suptt+ rmse_str + bias_str
        fig.suptitle( supt, fontsize=17, fontweight='bold')

    leg = plt.legend(loc='lower left',fontsize=22)
    # Set X/Y labels
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45,labelsize=15)
    ax.set_ylabel('Minimum Sea Level Pressure (hPa)',fontsize=24)
    ax.set_ylim(920,1020)  #( 900,1000 ) #(970,1020)

    #ax.set_title( 'mim SLP (hPa)',fontweight="bold",fontsize='15' )

    # Save the figure
    save_des = small_dir+Storm+'/'+Expers[1]+'/Vis_analyze/Model/minslp_'+DAtimes[0]+'_'+DAtimes[-1]+'_WSM6_THO.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


def Gather_slp( Storm, Expers, DAtimes, big_dir ):

    # Initialize container
    dict_minSLP = {}

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

    # Obtain min slp from WRF output
    for DAtime in DAtimes:
        for Exper in Expers:
            idx = Expers.index( Exper )
        #for DAtime in DAtimes:
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

    if slp_xa:
        dict_minSLP['xa_slp'] = np.array(xa_minslp_value)
        dict_minSLP['xa_lat'] = np.array(xa_minslp_lat)
        dict_minSLP['xa_lon'] = np.array(xa_minslp_lon)
    if slp_xb:
        dict_minSLP['xb_slp'] = np.array(xb_minslp_value)
        dict_minSLP['xb_lat'] = np.array(xb_minslp_lat)
        dict_minSLP['xb_lon'] = np.array(xb_minslp_lon)

    # Obtain assimilated min slp
    for DAtime in DAtimes:
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
            # find the index of the current time
            idx_time = DAtimes.index( DAtime )
            # assemble the diagnosed min slp from an analysis
            xa_ms_lat = xa_minslp_lat[0][idx_time]
            xa_ms_lon = xa_minslp_lon[0][idx_time]
            # ---condition 1: find the nearest TCvital min slp from the analysis
            # find the index/location of 'slp' in fort.10000
            indices = [i for i ,e in enumerate(list_enkf_minSlp) if e == 'slp']
            # assemble a pair of coordinate for each 'slp'
            distances = []
            obs_slp = []
            for it in indices:
                obs_slp.append( float(list_enkf_minSlp[it+9]) )
                lon1 = float(list_enkf_minSlp[it+2])
                lat1 = float(list_enkf_minSlp[it+1])
                distances.append( UD.mercator_distance(lon1, lat1, xa_ms_lon, xa_ms_lat) )
            min_distances = [np.amin(distances) == it for it in distances]
            idx_nearest = min_distances.index(True)
            # ---condition 2: the min slp
            min_obs_slp = [np.amin(obs_slp) == it for it in obs_slp]
            idx_min = min_obs_slp.index(True)
            # ---combine the two conditions
            if idx_min == idx_nearest:
                idx_coor = idx_min
                print('Storm center is choosed with two condtions met!')
            else:
                print('Not sure which obs is the storm center. Go with the min value one!')
                idx_coor = idx_min
            
            # gather this TCvital min slp
            obs_minslp_lat.append( float(list_enkf_minSlp[indices[idx_coor]+1]) )
            obs_minslp_lon.append( float(list_enkf_minSlp[indices[idx_coor]+2]) )
            obs_minslp_value.append( float(list_enkf_minSlp[indices[idx_coor]+9])/100 )

    # Assemble the dictionary
    dict_minSLP['obs_slp'] = np.array(obs_minslp_value)

    return dict_minSLP

# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot accumulated grid scale precipitation per snapshot
# ------------------------------------------------------------------------------------------------------
def read_Precip( Expers,wrf_dirs, DAtime ):

    d_pp_all = {}
    for iExper in Expers:
        idx_exper = Expers.index( iExper )
        # set file names for background and analysis
        mean_dir = [wrf_dirs[idx_exper]+DAtime+'/wrf_enkf_input_d03_mean',]     
        # read the shared variables: lat/lon
        ncdir = nc.Dataset( mean_dir[0] )
        lat = ncdir.variables['XLAT'][0,:,:]
        lon = ncdir.variables['XLONG'][0,:,:]
        # read precipitation
        pp = np.zeros( [len(mean_dir), np.size(lat,0), np.size(lat,1)] )
        for i in range(len(mean_dir)): 
            file_name = mean_dir[i] 
            ncdir = nc.Dataset( file_name )
            pp[i,:,:] = ncdir.variables['RAINNC'][0,:,:]
        d_pp = {'lat':lat,'lon':lon,'Precip':pp}
        
        d_pp_all[iExper] = d_pp
    return d_pp_all


def plot_Precip( Storm, Expers, DAtime, wrf_dirs, plot_dir ):

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
    d_precip = read_Precip( Expers, wrf_dirs, DAtime )

    # Big Assumption: in the same domain!!!
    lat = d_precip[Expers[0]]['lat']
    lon = d_precip[Expers[0]]['lon']
    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )

    # ------------------ Plot -----------------------
    fig, axs=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(12,5), dpi=400)

    # Exper 1: Xb
    min_pp = 0
    max_pp = 20
    bounds = np.linspace(min_pp, max_pp, 6)
    axs.flat[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    axs.flat[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    pp = d_precip[Expers[0]]['Precip'][0,:,:]
    #pp[pp == 0 ] = np.nan # set elements with 0 as np.nan
    mask = pp <= 0.5
    pp_masked = ma.masked_array(pp, mask=mask)
    pp_contourf = axs.flat[0].contourf(lon,lat,pp_masked,cmap='ocean_r',vmin=min_pp,vmax=max_pp,levels=bounds,extend='max',transform=ccrs.PlateCarree())
    # Mark the best track
    if if_btk_exist:
        axs.flat[0].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'red', marker='*',transform=ccrs.PlateCarree())

    # Exper 2: Xb
    axs.flat[1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    axs.flat[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    pp = d_precip[Expers[1]]['Precip'][0,:,:]
    pp[pp == 0 ] = np.nan # set elements with 0 as np.nan
    pp_contourf = axs.flat[1].contourf(lon,lat,pp,cmap='ocean_r',vmin=min_pp,vmax=max_pp,levels=bounds,extend='max',transform=ccrs.PlateCarree())
    # Mark the best track
    if if_btk_exist:
        axs.flat[1].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'red', marker='*',transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.12, 0.1, 0.5, 0.02])
    pp_bar = fig.colorbar(pp_contourf,ax=axs[0:2],orientation="horizontal", cax=caxes)
    pp_bar.ax.tick_params()

    # Exper 1: Xb - Exper 2: Xb
    min_incre = -5
    max_incre = 5
    axs.flat[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    axs.flat[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    pp_incre = d_precip[Expers[1]]['Precip'][:,:] - d_precip[Expers[0]]['Precip'][:,:]
    incre_pp = axs.flat[2].scatter(lon,lat,1.5,c=pp_incre,edgecolors='none', cmap='PiYG', vmin=min_incre, vmax=max_incre,transform=ccrs.PlateCarree())
    if if_btk_exist:
        axs.flat[2].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'black', marker='*',transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.65, 0.1, 0.25, 0.02])
    cb_diff_ticks = np.linspace(min_incre, max_incre, 5, endpoint=True)
    cbar = fig.colorbar(incre_pp, ax=axs[2:], ticks=cb_diff_ticks, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params()

    #subplot title
    axs.flat[0].set_title('IR_WSM6',fontweight='bold')
    axs.flat[1].set_title('IRMW_WSM6',fontweight='bold')
    axs.flat[2].set_title('IRMW-IR', fontweight='bold')

    #title for all
    fig.suptitle(Storm+': '+'('+DAtime+')'+'--Precip (mm)', fontsize=13, fontweight='bold')

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

    savename = plot_dir+DAtime+'_compare_Precip_masked.png'
    plt.savefig( savename, dpi=300 )
    print('Saving the figure: ', savename)
    plt.close()


def Precip( Storm, Expers, DAtimes, big_dir, small_dir ):

    wrf_dirs = []
    for iExper in Expers:
        wrf_dirs.append( big_dir+Storm+'/'+iExper+'/fc/' )
    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        print('Reading WRF background and analysis at ', DAtime)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Expers[0]+'/Vis_analyze/Model/Precip/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
            plot_Precip( Storm, Expers, DAtime, wrf_dirs, plot_dir )
        else:
            plot_Precip( Storm, Expers, DAtime, wrf_dirs, plot_dir )


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = ['CONV',]
    MP = ['WSM6','THO']

    slp_xa = True
    slp_xb = True
    Com_minslp_evo = True

    meanOverEns = False  # mean of H(ens) or H of mean(ens)
    Com_IC_water = False
    Com_Precip = True

    # Time range set up
    start_time_str = '201709030000'
    end_time_str = '201709040000'
    Consecutive_times = True
    # ------------------------------------------------------   

    # Create experiment names
    Expers = []
    for imp in MP:
        for ida in DA:
           Expers.append( UD.generate_one_name( Storm,ida,imp ) )

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

    # Plot the accumulated grid scale precipitation
    if Com_Precip:
        start_time=time.process_time()
        Precip( Storm, Expers, DAtimes, big_dir, small_dir )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')        
 


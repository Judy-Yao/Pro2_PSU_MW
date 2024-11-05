import os # functions for interacting with the operating system
import numpy as np
#import xarray as xr
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
#import metpy.calc as mpcalc
#from metpy.units import units
#import numpy.ma as ma

import Util_data as UD
import Diagnostics as Diag 

# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot psfc per snapshot
# ------------------------------------------------------------------------------------------------------


# Read the background ensemble: min PSFC
def read_mslp_ens( wrf_dir,DAtime ):

    id_ = []
    lon_mslp = []
    lat_mslp = []
    mslp = []

    with open( wrf_dir+ DAtime+'_enkf_input_mslp.txt' ) as f:
        next(f)
        all_lines = f.readlines()

    for line in all_lines:
        split_line = line.split()
        id_.append( split_line[0] )
        lon_mslp.append( float(split_line[1]) )
        lat_mslp.append( float(split_line[2]) )
        mslp.append( float(split_line[3]) )

    d_mslp = {'id':id_,'mslp':mslp,'lat_mslp':lat_mslp,'lon_mslp':lon_mslp}
    return d_mslp


def read_PSFC_mean( wrf_dir ):

    # read the shared variables: lat/lon
    ncdir = nc.Dataset( wrf_dir+'/wrf_enkf_input_d03_mean' )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    psfc = ncdir.variables['PSFC'][0,:,:]/100

    mpsfc =  np.nanmin( psfc )
    idx = np.nanargmin( psfc )
    lat_mpsfc =  lat.flatten()[idx]
    lon_mpsfc = lon.flatten()[idx]

    d_psfc = {'mslp':mpsfc,'lat_mslp':lat_mpsfc,'lon_mslp':lon_mpsfc,'lon':lon,'lat':lat,'psfc':psfc}
    return d_psfc


def read_PSFC_mem( mem ):

    file_name = os.path.basename( mem )
    id_ = file_name.split("_")[-1]

    ncdir = nc.Dataset( mem )
    psfc = ncdir.variables['PSFC'][0,:,:]/100

    d_psfc = {'id':id_,'psfc':psfc}
    return d_psfc


def plot_psfc( Storm, Exper_name, DAtime, wrf_dir, mem, plot_dir, d_xb_ens, ):

    # ------ Read WRFout -------------------
    d_mean = read_PSFC_mean( wrf_dir )
    d_mem = read_PSFC_mem( mem )

    # Locate the TCvitals mslp
    lon_obs = np.array( d_obs[DAtime]['lon'] )
    lat_obs = np.array( d_obs[DAtime]['lat'] )

    lat = d_mean['lat']
    lon = d_mean['lon']
    psfc_mean = d_mean['psfc']
    psfc_mem = d_mem['psfc']

    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )
    #min_psfc = np.min( psfc )
    #max_psfc = np.max( psfc )

    # ------ Plot Figure -------------------
    f, ax=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,2.5), dpi=400)

    # Set the map
    for i in range(3):
        ax.flat[i].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax.flat[i].coastlines(resolution='10m', color='black',linewidth=0.5)

    # background mean
    min_psfc = 995 #970
    max_psfc = 1015 #1015
    bounds = np.linspace(min_psfc, max_psfc, 6)

    ax[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax[0].coastlines (resolution='10m', color='black', linewidth=1)
    #ax[0].contourf(lon,lat,psfc_mean,cmap='inferno',levels=bounds,extend='both',transform=ccrs.PlateCarree())
    ax[0].contourf(lon,lat,psfc_mean,cmap='inferno',vmin=min_psfc,vmax=max_psfc,levels=bounds,extend='both',transform=ccrs.PlateCarree())

    # background member
    ax[1].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax[1].coastlines (resolution='10m', color='black', linewidth=1)
    cf_mem = ax[1].contourf(lon,lat,psfc_mem,cmap='inferno',vmin=min_psfc,vmax=max_psfc,levels=bounds,extend='both',transform=ccrs.PlateCarree())
    #cf_mem = ax[1].contourf(lon,lat,psfc_mem,cmap='inferno',transform=ccrs.PlateCarree())

    # Colorbar
    caxes = f.add_axes([0.12, 0.1, 0.5, 0.02])
    cbar = f.colorbar(cf_mem,ax=ax[0],orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)

    # Mme - mean
    min_diff = -5
    max_diff = 5
    # Create color bar with white at 0
    cmap = plt.get_cmap("bwr")
    #diff_ticks = np.linspace(min_diff, max_diff, 7, endpoint=True)
    #norm = mcolors.TwoSlopeNorm(vmin=min_diff, vcenter=0, vmax=max_diff)

    ax[2].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax[2].coastlines (resolution='10m', color='black', linewidth=1)
    cf_diff = ax[2].scatter(lon,lat,1,psfc_mem-psfc_mean,cmap=cmap,vmin=min_diff,vmax=max_diff,transform=ccrs.PlateCarree())

    # Adding the colorbar
    caxes = f.add_axes([0.65, 0.1, 0.25, 0.02])
    cb_diff_ticks = np.linspace(min_diff, max_diff, 7, endpoint=True)
    #cbar = f.colorbar(cf_diff, ax=ax[1:], ticks=cb_diff_ticks, orientation="horizontal", cax=caxes, extend='both')
    cbar = f.colorbar(cf_diff, ax=ax[1:],orientation="horizontal", cax=caxes, extend='both')
    cbar.ax.tick_params(labelsize=6)

    id_now = d_mem['id']
    idx = d_xb_ens['id'].index( id_now )
    for i in range(3):
        ax[i].scatter(lon_obs, lat_obs, c='green', s=1, marker='*', transform=ccrs.PlateCarree())
        ax[i].scatter(d_xb_ens['lon_mslp'][idx], d_xb_ens['lat_mslp'][idx], c='red', s=1.5, marker='*',  transform=ccrs.PlateCarree())
        ax[i].scatter(d_mean['lon_mslp'], d_mean['lat_mslp'], c='grey', s=1.5, marker='*',  transform=ccrs.PlateCarree())

    # Title
    #subplot title
    #matplotlib.rcParams['mathtext.fontset'] = 'custom'
    #matplotlib.rcParams['mathtext.bf'] = 'STIXGeneral:italic:bold'
    font = {'size':9,}
    ax[0].set_title('Xb: mean', font, )
    ax[1].set_title('Xb: '+d_mem['id'], font, )
    ax[2].set_title('Xb: '+d_mem['id'] +' - mean', font)

    f.suptitle(Storm+': '+Exper_name,fontsize=5, fontweight='bold')

    f.text( 0.05,0.86,'Assimilated obs: '+"{0:.2f}".format((d_obs[DAtime]['obs'][0]/100))+' hPa',fontsize=6,color='green',rotation='horizontal',fontweight='bold')
    f.text( 0.35,0.86,'Mslp in Xb mean: '+"{0:.2f}".format(d_mean['mslp'])+' hPa',fontsize=6,color='grey',rotation='horizontal',fontweight='bold')
    f.text( 0.68,0.86,'Mslp in Xb mem: '+"{0:.2f}".format(d_xb_ens['mslp'][idx])+' hPa',fontsize=6,color='red',rotation='horizontal',fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2))
    for j in range(3):
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
        gl.xlabel_style = {'size': 4}
        gl.ylabel_style = {'size': 6}

    des_path = plot_dir+DAtime+'_'+d_mem['id']+'_psfc.png'
    plt.savefig( des_path, dpi=300 )
    print('Saving the figure: ', des_path)
    plt.close()


def PSFC( Storm, Exper_name, DAtime, big_dir, small_dir ):

    wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
    
    # Read the minimum PSFC for ensemble enkf input
    d_xb_ens = read_mslp_ens( wrf_dir,DAtime )

    files_yb = sorted( glob.glob(wrf_dir + '/wrf_enkf_input_d03_0*') )
    # Loop through each DAtime/analysis
    for file in files_yb:
        print('Reading WRF background member: ', file)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/Ens_PSFC/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
        
        plot_psfc( Storm, Exper_name, DAtime, wrf_dir, file, plot_dir, d_xb_ens, )



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # -------- Configuration -----------------
    Storm = 'IRMA'
    DA = 'CONV'
    MP = 'THO'
    fort_v = ['obs_type','lat','lon','obs']

    # Time range set up
    start_time_str = '201709030000'
    end_time_str = '201709030000'
    Consecutive_times = True

    Plot_PSFC = True

    # -----------------------------------------

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
    Exper_name = UD.generate_one_name( Storm,DA,MP ) 

    # Read assimilated TCvital location
    d_obs = {}
    print('------------ Read obs info from enkf diagnostics fort.10000 --------------')
    for DAtime in DAtimes:
        # Read assimilated obs
        file_Diag = big_dir+Storm+'/'+Exper_name+'/run/'+DAtime+'/enkf/d03/fort.10000'
        d_obs[DAtime] = Diag.Find_min_slp( file_Diag, fort_v )


    # Plot PSFC (surface pressure)
    if Plot_PSFC:
        
        for DAtime in DAtimes:
            start_time=time.process_time()
            PSFC( Storm, Exper_name, DAtime, big_dir, small_dir )
            end_time = time.process_time()
            print ('time needed: ', end_time-start_time, ' seconds')





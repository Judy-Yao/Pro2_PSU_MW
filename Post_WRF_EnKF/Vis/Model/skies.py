#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

from numba import njit, prange
import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
from scipy import interpolate
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from matplotlib import colors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import pickle

import Read_Obspace_IR as ROIR

# ------------------------------------------------------------------------------------------------------
#           Operation: Perform Operations
# ------------------------------------------------------------------------------------------------------

# use boolean(0/1) to indicate if a column of an ensemble member is cloudy or not
@njit(parallel=True)
def boolean_skies( ens_allHydro ):
    res = np.zeros( (ens_allHydro.shape[0],ens_allHydro.shape[2]),  )
    for n in prange( ens_allHydro.shape[2] ):
        for m in range( ens_allHydro.shape[0] ):
            for l in range( ens_allHydro.shape[1] ):
                if ens_allHydro[m,l,n] > 1e-6: # once found a layer meets the threshold, skip the rest of l loop
                    res[m,n] = 1
                    break
                else:
                    res[m,n] = 0
    return res

# Sum up all hydrometeor mixing ratios at a grid point
def allHydro( DAtime,wrf_dir ):
    
    nLevel = 42
    xb_ens_allHydro = np.zeros( shape=[num_ens,nLevel,xmax*ymax] )
    file_xb = sorted( glob.glob(wrf_dir + '/wrf_enkf_input_d03_0*') )
    hydros = ['QCLOUD','QRAIN','QICE','QSNOW','QGRAUP']
    for ifile in file_xb:
        total_hydro = np.zeros( shape=[nLevel,xmax*ymax] )
        ncdir = nc.Dataset( ifile, 'r')
        for ih in hydros:
            var = ncdir.variables[ih][0,:,:,:]
            total_hydro = total_hydro+var.reshape(nLevel,xmax*ymax) 
        idx = file_xb.index( ifile )
        xb_ens_allHydro[idx,:,:] = total_hydro

    # May save the result
    if If_save:
        des_path = wrf_dir+ "xb_d03_3D_ens_allHydro_" + DAtime + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( xb_ens_allHydro, f )
        f.close()
        print('Save '+des_path) 

    return None

# Calculate cloud probability based on ensemble members
def cloud_probability( DAtime,wrf_dir ):

    # Read ensemble's total hydro mixting ratio
    des_path = wrf_dir+ "xb_d03_3D_ens_allHydro_" + DAtime + '.pickle'
    with open( des_path,'rb' ) as f:
        xb_ens_allHydro = pickle.load( f )
    print('Shape of xb_ens_allHydro: '+ str(np.shape(xb_ens_allHydro)))

    # Use 0/1 to indicate if a column of a member is clear or cloudy
    ens_col_sky = boolean_skies( xb_ens_allHydro )
    cloud_p = np.sum( ens_col_sky,0 )/num_ens

    #----------  Plot ---------------
    # Read WRF domain
    wrf_input_mean = wrf_dir+'wrf_enkf_input_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_input_mean )
    ncdir = nc.Dataset( wrf_input_mean, 'r')
    xlat = ncdir.variables['XLAT'][0,:,:]
    xlon = ncdir.variables['XLONG'][0,:,:]

    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, figsize=(6,6), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    # Define probability value threshold
    min_p = 0
    max_p = 1
    bounds = [0.01,0.1,0.3,0.5,0.7,0.9,0.99]#np.linspace(min_p, max_p, 6)

    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax.contourf(xlon, xlat,cloud_p.reshape( xmax,ymax ),\
                cmap='terrain_r',levels=bounds,extend='both',transform=ccrs.PlateCarree())    
    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    color_ticks = [0.1,0.3,0.5,0.7,0.9]
    cbar = fig.colorbar(cs, cax=cbaxes,ticks=color_ticks,fraction=0.046, pad=0.04)
    bounds_str =  [ str(item) for item in color_ticks ]
    cbar.ax.set_xticklabels( bounds_str, rotation=45)
    cbar.ax.tick_params(labelsize=12)

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(1):
        gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')

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

    #title 
    ax.set_title('Cloud Probability',fontsize=12, fontweight='bold')
    fig.suptitle(Storm+': '+Exper_name, fontsize=10, fontweight='bold')

    # Save the figure
    des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/cloud_probability_'+DAtime+'.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)
    plt.close()
    return None




if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'HARVEY'
    Exper_name = 'JerryRun/IR_THO'

    v_interest = [ 'QVAPOR',]
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201708221200'
    end_time_str = '201708221200'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    total_hydro = True
    cloud_P = True
    If_save = True

    # -------------------------------------------------------    

    if not Consecutive_times:
        DAtimes = ['201709160600','201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    
    # Sum up all hydrometeor mixing ratios at a grid point
    if total_hydro:
        start_time=time.process_time()
        print('Sum up mixing ratios of all hydrometeors at a grid point------')
        for DAtime in DAtimes:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            allHydro( DAtime,wrf_dir )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')

    # cloud probability per snapshot
    if cloud_P:
        start_time=time.process_time()
        print('Calculating the cloud probability based on the whole ensemble------')
        for DAtime in DAtimes:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            cloud_probability( DAtime,wrf_dir )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')

































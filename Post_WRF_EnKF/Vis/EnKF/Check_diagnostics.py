#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import pickle
import netCDF4 as nc
import Diagnostics as Diag
import math
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import Util_Vis

def read_wrf_domain( wrf_file ):

    print('Read domain info from: ' + wrf_file)
    ncdir = nc.Dataset(wrf_file, 'r')

    Lat_x = ncdir.variables['XLAT'][0,:,:] #latitude: XLAT(time, y, x)
    Lon_x = ncdir.variables['XLONG'][0,:,:] #longitude: XLONG(time, y, x)

    lat_min = np.min( Lat_x.flatten() )
    lat_max = np.max( Lat_x.flatten() )
    lon_min = np.min( Lon_x.flatten() )
    lon_max = np.max( Lon_x.flatten() )

    d03_list = {'lat_min':lat_min, 'lat_max':lat_max, 'lon_min':lon_min, 'lon_max':lon_max}
    return d03_list



def Plot_IR(Storm, DAtime, Exper, Diag_obs_types):

    # Read WRF domain
    wrf_file = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/JerryRun/MW_THO/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = read_wrf_domain( wrf_file )
    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']
    
    # ------------------ Plot ----------------------
    # Set up figure
    fig = plt.figure( figsize=(6,4), dpi=150 )
    gs = fig.add_gridspec(2,3)

    # Plot Tb values

    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )

    obs_kinds = ['obs','prior_mean','posterior_mean']
    Lat_obs = Diag_obs_types['IR']['lat']
    Lon_obs = Diag_obs_types['IR']['lon']
    for j in range(len(obs_kinds)):
        ax0 = fig.add_subplot( gs[0,j],  projection=ccrs.PlateCarree())
        ax0.set_extent( [lon_min,lon_max,lat_min,lat_max],  crs=ccrs.PlateCarree())
        ax0.coastlines( resolution='10m', color='black',linewidth=0.5 )
        ax0.scatter(Lon_obs,Lat_obs,1.5,c=Diag_obs_types['IR'][obs_kinds[j]],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())

    # Plot spread

    # Define spread threshold
    min_S = 0
    max_S = 35
    spread_kinds = ['prior_spread','posterior_spread']
    for j in range(len(spread_kinds)):
        ax0 = fig.add_subplot( gs[1,j+1],  projection=ccrs.PlateCarree())
        ax0.set_extent( [lon_min,lon_max,lat_min,lat_max],  crs=ccrs.PlateCarree())
        ax0.coastlines( resolution='10m', color='black',linewidth=0.5 )
        cs = ax0.scatter(Lon_obs,Lat_obs,1.5,c=Diag_obs_types['IR'][spread_kinds[j]],edgecolors='none', cmap='jet', vmin=min_S, vmax=max_S,transform=ccrs.PlateCarree())

    #subplot title
    #font = {'size':8,}
    #ax[0].set_title('Yo Assimilated', font, fontweight='bold')
    #ax[1].set_title('H(Xb) Diverged', font, fontweight='bold')
    #ax[2].set_title('H(Xa) Diverged', font, fontweight='bold')
    plt.savefig('test.png', dpi=300)



def Diagnose(big_dir, Storm, Exper, DAtime, v_interest, Check_IR, Check_MW):

    # Print out obs types
    file_diag = big_dir+Storm+'/'+Exper+'/run/'+DAtime+'/enkf/d03/fort.10000' 
    Diag.return_obs_type( file_diag )

    # Only check IR and MW obs
    key_obstype = ['IR', 'MW']
    Diag_obs_types = {key: None for key in key_obstype}
    #if Check_MW:
    #    MW_so = big_dir+Storm+'/'+Exper+'/run/'+DAtime+'/enkf/d03/microwave_'+DAtime+'_so' # Read in MW diagnostics with labelled sensor info for each experiment
    #    Diag_obs_types['MW'] = Diag.label_mw_obs( file_diag, MW_so, v_interest )

        #Plot_MW( Storm, time, Exper, MW_diag_Expers )

    if Check_IR:
        Diag_obs_types['IR'] = Diag.Find_IR( file_diag, v_interest )

        Plot_IR( Storm, DAtime, Exper, Diag_obs_types )


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'

    # Configuration
    Storm = 'HARVEY'
    Exper_name = 'IR+MW-J_DA+J_WRF+J_init-SP-intel19'
    v_interest = ['obs_type','lat','lon','obs','prior_mean','posterior_mean','prior_spread','posterior_spread'] 
    start_time_str = '201708221200'
    end_time_str = '201708221200'
    Consecutive_times = True

    # Identify DA times in the period of interest
    if not Consecutive_times:
        DAtimes = ['201709140000',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Loop through each available fort.10000 
    for DAtime in DAtimes:
        Check_IR = False
        Check_MW = False
        if os.path.exists( big_dir+Storm+'/'+Exper_name+'/run/'+DAtime+'/enkf/d03/radiance_'+DAtime+'_so'):
            Check_IR = True
        if os.path.exists( big_dir+Storm+'/'+Exper_name+'/run/'+DAtime+'/enkf/d03/microwave_'+DAtime+'_so'):
            Check_MW = True
    
        Diagnose(big_dir, Storm, Exper_name, DAtime, v_interest, Check_IR, Check_MW)





















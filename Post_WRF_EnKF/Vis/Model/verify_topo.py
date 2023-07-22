#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from matplotlib import colors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time

import Read_Obspace_IR as ROIR

def Plot_snapshot( wrf_file ):

    # Read files
    ncdir = nc.Dataset( wrf_file, 'r')
    xlon = ncdir.variables['XLONG'][0,:,:]
    xlat = ncdir.variables['XLAT'][0,:,:]
    if plot_geoH:
        PHB = ncdir.variables['PHB'][0,:,:,:]
        PH = ncdir.variables['PH'][0,:,:,:]
        geoHkm = (PHB+PH)/9.8/1000
        geoHkm = geoHkm.reshape( geoHkm.shape[0],-1)
        var_all = geoHkm
        sample = [0,4,8,14,16,18]
    elif plot_pres:
        PB = ncdir.variables['PB'][0,:,:,:]
        P = ncdir.variables['P'][0,:,:,:]
        P_hpa = (PB + P)/100
        P_hpa = P_hpa.reshape( P_hpa.shape[0],-1) 
        var_all = P_hpa
        sample = [0,4,8,14,16,18]
    else:
        pass

    # Read WRF domain
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = ROIR.read_TCvitals(small_dir+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(2, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(9.75,6.5), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    # Define the colorbar
    if plot_geoH:
        max_var = 4 #max(np.amin(abs(Interp_incre)),np.amax(abs(Interp_incre)))
        min_var = 0#0-max_abs
    elif plot_pres:
        max_var = 1000
        min_var = 300
        bounds = np.linspace(min_var, max_var, 8)
    for isub in range(6):
        ax.flat[isub].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.flat[isub].coastlines(resolution='10m', color='black',linewidth=0.5)
        var = var_all[sample[isub],:].reshape( (xlon.shape[0],xlon.shape[1]) )
        cs = ax.flat[isub].contourf(xlon,xlat,var,cmap='terrain_r',vmin=min_var,vmax=max_var,levels=bounds,transform=ccrs.PlateCarree(),extend='both')
        #c = 1.5
        #cs = ax.flat[isub].scatter(xlon,xlat,c,var,cmap='terrain_r',edgecolors='none',vmin=min_var,vmax=max_var,transform=ccrs.PlateCarree(),)

        if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
            ax.flat[isub].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, )
    cbar.set_clim( vmin=min_var, vmax=max_var )
    cbar.ax.tick_params(labelsize=12)

    #subplot title
    font = {'size':15,}
    for isub in range(6):
        ax.flat[isub].set_title( 'ML: '+str(sample[isub]), font, fontweight='bold')

    #title for all
    if plot_geoH:
        fig.suptitle(Storm+':'+Exper_name+' Geo Height (km)', fontsize=10, fontweight='bold')
    elif plot_pres:
        fig.suptitle('noIR  '+Storm+':'+Exper_name+' Pressure (hPa)', fontsize=10, fontweight='bold')
    else:
        pass

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for j in range(6):
        gl = ax.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top = False
        gl.xlabels_bottom = True
        if j==0 or j==3:
            gl.ylabels_left = True
            gl.ylabels_right = False
        else:
            gl.ylabels_left = False
            gl.ylabels_right = False

        if j==3 or j==4 or j==5:
            gl.xlabels_bottom = True
            gl.xlabels_top = False
        else:
            gl.xlabels_bottom = False
            gl.xlabels_top = False

        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10}
        gl.ylabel_style = {'size': 12}

    # Save the figure
    if plot_geoH:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/ML_geoH'+'_'+DAtime+'.png'
    elif plot_pres:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/ML_pressure'+'_'+DAtime+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()




if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    Exper_name = 'IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'

    start_time_str = '201708221200'
    end_time_str = '201708221200'
    Consecutive_times = False

    plot_geoH = False
    plot_pres = True
    # -------------------------------------------------------    
    
    # Identify DA times in the period of interest
    if not Consecutive_times:
        DAtimes = ['201709041600',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]


    # Plot the snapshot
    for DAtime in DAtimes:
        wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
        print('Plotting ',wrf_file)
        Plot_snapshot( wrf_file )









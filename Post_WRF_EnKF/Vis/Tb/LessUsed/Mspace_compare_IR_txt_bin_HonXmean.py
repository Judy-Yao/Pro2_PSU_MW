#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os
import glob
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math
from datetime import datetime, timedelta
import time

import Util_Vis
import Diagnostics as Diag
import Util_data as UD


# Read simulated data on model resolution
def Read_IR_modelRes_HonXmean( DAtime ):
   
    Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime

    # Dimension of the domain
    xmax = 297
    ymax = 297

    # List the Yb and Ya files
    file_yb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem_mean*.bin') )
    file_ya = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_output_mem_mean*.bin') )   

    # Read attributes from a member
    # yb
    yb_control = np.fromfile( file_yb[0],dtype='<f4') # <: little endian; f: float; 4: 4 bytes
    n_ch = len(yb_control)/(xmax*ymax) - 2
    n_ch = int(n_ch)
    if n_ch != len(ch_list):
        print('Error!! # of channels in data is '+str(n_ch))
    yb_data = yb_control[:].reshape(n_ch+2,ymax,xmax)
    lon_x =  yb_data[0,:,:].flatten() 
    lat_x = yb_data[1,:,:].flatten() 
    HonXb = yb_data[2,:,:].flatten() 
    # ya
    ya_control = np.fromfile( file_ya[0],dtype='<f4') # <: little endian; f: float; 4: 4 bytes
    ya_data = ya_control[:].reshape(n_ch+2,ymax,xmax)
    HonXa = ya_data[2,:,:].flatten() 

    lat_x = np.array( lat_x )
    lon_x = np.array( lon_x )
    HonXb = np.array( HonXb )
    HonXa = np.array( HonXa )
    d_IRx = {'Lat_x':lat_x, 'Lon_x':lon_x, 'HonXb':HonXb, 'HonXa':HonXa}
    return d_IRx


def plot_Tb( DAtime, d_obs ):

    # Read simulated IR Tbs
    d_simu = Read_IR_modelRes_HonXmean( DAtime ) 

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(Storm, DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # ------------------ Plot -----------------------
    f, ax=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,2.5), dpi=400)

    # Define the domain
    lat_min = np.amin(d_simu['Lat_x'])
    lat_max = np.amax(d_simu['Lat_x'])
    lon_min = np.amin(d_simu['Lon_x'])
    lon_max = np.amax(d_simu['Lon_x'])
  
    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )
    
    ax[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    ax[0].scatter(d_obs['lon'],d_obs['lat'],1.5,c=d_obs['obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())

    ax[1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    ax[1].scatter(d_simu['Lon_x'], d_simu['Lat_x'],1,c=d_simu['HonXb'],\
                edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax[1].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())

    ax[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax[2].scatter(d_simu['Lon_x'], d_simu['Lat_x'],1,c=d_simu['HonXa'],\
                edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax[2].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())

    # Colorbar
    caxes = f.add_axes([0.4, 0.1, 0.5, 0.02])
    cbar = f.colorbar(cs, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)

    #title for all
    f.suptitle(Storm+': '+Exper_name, fontsize=8, fontweight='bold')

    #subplot title
    font = {'size':8,}
    ax[0].set_title('Yo', font, fontweight='bold')
    ax[1].set_title(r'$H_{Tb}(Xb)$', font, fontweight='bold')
    ax[2].set_title(r'$H_{Tb}(Xa)$', font, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for j in range(3):
        gl = ax[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
       
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
        gl.xlabel_style = {'size': 4}
        gl.ylabel_style = {'size': 6}

    des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IR_HonXmean/'+DAtime+'_'+sensor+'_mspace_HonXmean.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'HARVEY'
    DA = 'IR'
    MP = 'TuneWSM6'

    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201708221200'
    end_time_str = '201708221200'
    Consecutive_times = True

    If_plot = True
    plot_scatter = True
    # -------------------------------------------------------   

    # Create experiment names
    Exper_name =  UD.generate_one_name( Storm,DA,MP )

    if not Consecutive_times:
        IR_times = ['201709041600',]
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
   
    # Plot
    if If_plot:
        for DAtime in IR_times:
            # Read assimilated obs
            file_Diag = big_dir+Storm+'/'+Exper_name+'/run/'+DAtime+'/enkf/d03/fort.10000'
            d_obs = Diag.Find_IR( file_Diag, fort_v )

            print('------------ Plot ----------------------')
            print('DAtime: '+ DAtime)
            plot_Tb( DAtime, d_obs )
            #plot_Tb( Storm, Exper_name, Hx_dir, DAtime, sensor, ch_list, d_obs)
        

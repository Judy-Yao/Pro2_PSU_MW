#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os
import glob
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math
from datetime import datetime, timedelta
import time

import Util_Vis
import Diagnostics as Diag
import Util_data as UD

def RMSE(simu, obs):
    return np.sqrt( ((simu - obs) ** 2).mean() )

def Bias(simu, obs):
    return  np.sum((simu - obs),0)/np.size(obs,0)

def mean_Yo_Hx(simu, obs):
    return  np.sum((obs - simu),0)/np.size(obs,0)

def plot_Tb( Storm, Exper_name, Hxb, Hxa, DAtime, sensor, ch_list, d_obs):

    # Read simulated data
    d_simu = UD.read_simu_IR_single( Hxb, Hxa, ch_list )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir, Storm, DAtime)
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
    cs = ax[1].scatter(d_simu['Lon_x'], d_simu['Lat_x'],1,c=d_simu['Yb_x'],\
                edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax[1].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())

    ax[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax[2].scatter(d_simu['Lon_x'], d_simu['Lat_x'],1,c=d_simu['Ya_x'],\
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

    head_tail = os.path.split( Hxb )
    mem = head_tail[1].replace('TB_GOES_CRTM_input_mem','')
    mem = mem.replace('_d03_2017-09-03_00:00.bin','')
    des_name = small_dir+'/Clean_results/'+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IR_60mem/'+DAtime+'_'+mem+'_mspace.png'
    plt.savefig( des_name, dpi=300)
    plt.close()
    print('Saving the figure: ', des_name)


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    Exper_name = UD.generate_one_name( Storm,'IR','WSM6' )
    Exper_obs =  UD.generate_one_name( Storm,'IR','WSM6' )
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']
    num_ens = 60

    start_time_str = '201709030100'
    end_time_str = '201709030100'
    Consecutive_times = True

    If_plot = True
    If_plot_diff = False
    # -------------------------------------------------------   

    if not Consecutive_times:
        IR_times = ['201709050000','201709050600','201709051200','201709051800','201709060000',]
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
   
    # Plot
    if If_plot:
        for DAtime in IR_times:
            print('DAtime: '+ DAtime)
            # Read assimilated obs
            file_Diag = big_dir+Storm+'/'+Exper_obs+'/run/'+DAtime+'/enkf/d03/fort.10000'
            d_obs = Diag.Find_IR( file_Diag, fort_v )

            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            file_yb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem0*.bin') )
            file_ya = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_output_mem0*.bin') )
            print('------------ Plot ----------------------')
            for i in range(len(file_yb)):
                plot_Tb( Storm, Exper_name, file_yb[i], file_ya[i], DAtime, sensor, ch_list, d_obs)
       
    # Plot
    if If_plot_diff:
        for DAtime in IR_times:
            print('DAtime: '+ DAtime)
            # Read assimilated obs
            file_Diag = big_dir+Storm+'/'+Exper_obs+'/run/'+DAtime+'/enkf/d03/fort.10000'
            d_obs = Diag.Find_IR( file_Diag, fort_v )

            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            file_yb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem0*.bin') )
            file_ya = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_output_mem0*.bin') )
            print('------------ Plot ----------------------')
            for i in range(len(file_yb)):
                plot_Tb_diff( Storm, Exper_name, file_yb[i], file_ya[i], DAtime, sensor, ch_list, d_obs)


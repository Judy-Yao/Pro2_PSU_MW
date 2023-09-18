#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os
import glob
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as patches
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from global_land_mask import globe
import math
from datetime import datetime, timedelta
import time
from scipy.interpolate import interp2d

import Util_data as UD
import Util_Vis
import Diagnostics as Diag

# Plotting domain
def read_wrf_domain( DAtimes ):

    d_allD03 = {}
    for DAtime in DAtimes:
        d_allD03[DAtime] = {}

    for DAtime in DAtimes:
        for iExper in Expers:
            wrf_file = big_dir+Storm+'/'+iExper+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
            print('Read domain info from: ' + wrf_file)
            ncdir = nc.Dataset(wrf_file, 'r')

            Lat_x = ncdir.variables['XLAT'][0,:,:] #latitude: XLAT(time, y, x)
            Lon_x = ncdir.variables['XLONG'][0,:,:] #longitude: XLONG(time, y, x)

            lat_min = np.min( Lat_x.flatten() )
            lat_max = np.max( Lat_x.flatten() )
            lon_min = np.min( Lon_x.flatten() )
            lon_max = np.max( Lon_x.flatten() )
            d_allD03[DAtime][iExper] = {'lat_min':lat_min, 'lat_max':lat_max, 'lon_min':lon_min, 'lon_max':lon_max}

    return d_allD03


def read_allTb( ):

    d_allIR = {}

    for DAtime in DAtimes:
        d_allIR[DAtime] = {}

    for DAtime in DAtimes:
        for iExper in Expers:
            idx_exper = Expers.index( iExper )
            Tb_file = big_dir+Storm+'/'+iExper+'/Obs_Hx/IR/'+DAtime+"/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    
            lat_obs = []
            lon_obs = []
            Yo_obs = []
            meanYb_obs = []
            meanYa_obs = []

            # Read records
            print('Reading ', Tb_file)
            with open(Tb_file) as f:
                next(f)
                all_lines = f.readlines()

            for line in all_lines:
                split_line = line.split()
                lat_obs.append( float(split_line[0]) )
                lon_obs.append( float(split_line[1]) )
                Yo_obs.append( float(split_line[3]) )
                meanYb_obs.append( float(split_line[4]) )
                meanYa_obs.append( float(split_line[5]) )

            lat_obs = np.array( lat_obs )
            lon_obs = np.array( lon_obs )
            Yo_obs = np.array( Yo_obs )
            meanYb_obs = np.array( meanYb_obs )
            meanYa_obs = np.array( meanYa_obs )
            d_Tb =  {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'Yo_obs':Yo_obs, 'meanYb_obs':meanYb_obs, 'meanYa_obs':meanYa_obs}
        
            d_allIR[DAtime][iExper] = d_Tb
        
    return d_allIR

def plot_Tb_withObs( ):

    # Read WRF domain
    d03 = read_wrf_domain( DAtimes )     
    
    # Read IR
    d_allIR = read_allTb(  )

    # ------------------ Plot -----------------------
    #fig, axs=plt.subplots(5, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5,  figsize=(6,10), dpi=400)
    fig, axs=plt.subplots(5, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(6,10), dpi=400)

    c_size = 3
    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )

    for DAtime in DAtimes:
        j = DAtimes.index( DAtime )
        for i in range(5):
            if i == 0:
                axs[i,j].set_extent([d03[DAtime][Expers[0]]['lon_min'],d03[DAtime][Expers[0]]['lon_max'],d03[DAtime][Expers[0]]['lat_min'],d03[DAtime][Expers[0]]['lat_max']], crs=ccrs.PlateCarree())
                axs[i,j].coastlines(resolution='10m', color='black',linewidth=0.5)
                cs = axs[i,j].scatter(d_allIR[DAtime][Expers[0]]['lon_obs'],d_allIR[DAtime][Expers[0]]['lat_obs'],c_size,d_allIR[DAtime][Expers[0]]['Yo_obs'],edgecolors='none',cmap=IRcmap,vmin=min_T,vmax=max_T,transform=ccrs.PlateCarree())
            elif i == 1:
                axs[i,j].set_extent([d03[DAtime][Expers[0]]['lon_min'],d03[DAtime][Expers[0]]['lon_max'],d03[DAtime][Expers[0]]['lat_min'],d03[DAtime][Expers[0]]['lat_max']], crs=ccrs.PlateCarree())
                axs[i,j].coastlines(resolution='10m', color='black',linewidth=0.5)
                axs[i,j].scatter(d_allIR[DAtime][Expers[0]]['lon_obs'],d_allIR[DAtime][Expers[0]]['lat_obs'],c_size,d_allIR[DAtime][Expers[0]]['meanYb_obs'],edgecolors='none',cmap=IRcmap,vmin=min_T,vmax=max_T,transform=ccrs.PlateCarree())
            elif i == 2:
                axs[i,j].set_extent([d03[DAtime][Expers[0]]['lon_min'],d03[DAtime][Expers[0]]['lon_max'],d03[DAtime][Expers[0]]['lat_min'],d03[DAtime][Expers[0]]['lat_max']], crs=ccrs.PlateCarree())
                axs[i,j].coastlines(resolution='10m', color='black',linewidth=0.5)
                axs[i,j].scatter(d_allIR[DAtime][Expers[0]]['lon_obs'],d_allIR[DAtime][Expers[0]]['lat_obs'],c_size,d_allIR[DAtime][Expers[0]]['meanYa_obs'],edgecolors='none',cmap=IRcmap,vmin=min_T,vmax=max_T,transform=ccrs.PlateCarree())
            elif i == 3:
                axs[i,j].set_extent([d03[DAtime][Expers[1]]['lon_min'],d03[DAtime][Expers[1]]['lon_max'],d03[DAtime][Expers[1]]['lat_min'],d03[DAtime][Expers[1]]['lat_max']], crs=ccrs.PlateCarree())
                axs[i,j].coastlines(resolution='10m', color='black',linewidth=0.5)
                axs[i,j].scatter(d_allIR[DAtime][Expers[1]]['lon_obs'],d_allIR[DAtime][Expers[1]]['lat_obs'],c_size,d_allIR[DAtime][Expers[1]]['meanYb_obs'],edgecolors='none',cmap=IRcmap,vmin=min_T,vmax=max_T,transform=ccrs.PlateCarree())
            elif i == 4:
                axs[i,j].set_extent([d03[DAtime][Expers[1]]['lon_min'],d03[DAtime][Expers[1]]['lon_max'],d03[DAtime][Expers[1]]['lat_min'],d03[DAtime][Expers[1]]['lat_max']], crs=ccrs.PlateCarree()) 
                axs[i,j].coastlines(resolution='10m', color='black',linewidth=0.5)
                axs[i,j].scatter(d_allIR[DAtime][Expers[1]]['lon_obs'],d_allIR[DAtime][Expers[1]]['lat_obs'],c_size,d_allIR[DAtime][Expers[1]]['meanYa_obs'],edgecolors='none',cmap=IRcmap,vmin=min_T,vmax=max_T,transform=ccrs.PlateCarree())


    # Colorbar
    caxes = fig.add_axes([0.13, 0.06, 0.75, 0.02])
    tb_range = np.linspace(min_T,max_T,8)
    cbar = fig.colorbar(cs,ticks=tb_range,orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=13)

    #subplot title
    font = {'size':12,}
    for DAtime in DAtimes:
        idx_t = DAtimes.index( DAtime )
        axs[0,idx_t].set_title(DAtime, font, fontweight='bold')

    #title for all
    fig.suptitle(Storm, fontsize=8, fontweight='bold')

    # ylabel
    axs[0,0].set_ylabel('Obs',fontsize=8)
    axs[1,0].set_ylabel('THO_H(Xb)',fontsize=12, fontweight='bold')    
    axs[2,0].set_ylabel('THO_H(Xa)',fontsize=12, fontweight='bold')
    axs[3,0].set_ylabel('WSM6_H(Xb)',fontsize=12, fontweight='bold')
    axs[4,0].set_ylabel('WSM6_H(Xa)',fontsize=12, fontweight='bold')

    # Axis labels
    for j in range(3):
        for i in range(5):
            gl = axs[i,j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
            if i == 1 or i == 2 or i == 3 :
                lon_min = d03[DAtime][Expers[0]]['lon_min']
                lon_max = d03[DAtime][Expers[0]]['lon_max']
                lat_min = d03[DAtime][Expers[0]]['lat_min']
                lat_max = d03[DAtime][Expers[0]]['lat_max']
            else:
                lon_min = d03[DAtime][Expers[1]]['lon_min']
                lon_max = d03[DAtime][Expers[1]]['lon_max']
                lat_min = d03[DAtime][Expers[1]]['lat_min']
                lat_max = d03[DAtime][Expers[1]]['lat_max']

            lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
            lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

            if j==0:
                gl.ylabels_left = True
                gl.ylabels_right = False
            else:
                gl.ylabels_left = False
                gl.ylabels_right = False
    
            if i == 4:
                gl.xlabels_top = False
                gl.xlabels_bottom = True    
            else:
                gl.xlabels_top = False
                gl.xlabels_bottom = False

            gl.ylocator = mticker.FixedLocator(lat_ticks)
            gl.xlocator = mticker.FixedLocator(lon_ticks)
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 8}
            gl.ylabel_style = {'size': 10}

    des_name = 'IR_withObs.png'
    #des_name = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Expers[0]+'/Vis_analyze/Tb/IR/Obspace/'+DAtime+'_'+sensor+'_Obspace_expers.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)


def plot_Tb_withoutObs( ):

    # Read WRF domain
    d03 = read_wrf_domain( DAtimes )

    # Read IR
    d_allIR = read_allTb(  )

    # ------------------ Plot -----------------------
    fig, axs=plt.subplots(4, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(6,8), dpi=400)

    c_size = 3
    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )

    for DAtime in DAtimes:
        j = DAtimes.index( DAtime )
        for i in range(4):
            if i == 0:
                axs[i,j].set_extent([d03[DAtime][Expers[0]]['lon_min'],d03[DAtime][Expers[0]]['lon_max'],d03[DAtime][Expers[0]]['lat_min'],d03[DAtime][Expers[0]]['lat_max']], crs=ccrs.PlateCarree())
                axs[i,j].coastlines(resolution='10m', color='black',linewidth=0.5)
                cs = axs[i,j].scatter(d_allIR[DAtime][Expers[0]]['lon_obs'],d_allIR[DAtime][Expers[0]]['lat_obs'],c_size,d_allIR[DAtime][Expers[0]]['meanYb_obs'],edgecolors='none',cmap=IRcmap,vmin=min_T,vmax=max_T,transform=ccrs.PlateCarree())
            elif i == 1:
                axs[i,j].set_extent([d03[DAtime][Expers[0]]['lon_min'],d03[DAtime][Expers[0]]['lon_max'],d03[DAtime][Expers[0]]['lat_min'],d03[DAtime][Expers[0]]['lat_max']], crs=ccrs.PlateCarree())
                axs[i,j].coastlines(resolution='10m', color='black',linewidth=0.5)
                axs[i,j].scatter(d_allIR[DAtime][Expers[0]]['lon_obs'],d_allIR[DAtime][Expers[0]]['lat_obs'],c_size,d_allIR[DAtime][Expers[0]]['meanYa_obs'],edgecolors='none',cmap=IRcmap,vmin=min_T,vmax=max_T,transform=ccrs.PlateCarree())
            elif i == 2:
                axs[i,j].set_extent([d03[DAtime][Expers[1]]['lon_min'],d03[DAtime][Expers[1]]['lon_max'],d03[DAtime][Expers[1]]['lat_min'],d03[DAtime][Expers[1]]['lat_max']], crs=ccrs.PlateCarree())
                axs[i,j].coastlines(resolution='10m', color='black',linewidth=0.5)
                axs[i,j].scatter(d_allIR[DAtime][Expers[1]]['lon_obs'],d_allIR[DAtime][Expers[1]]['lat_obs'],c_size,d_allIR[DAtime][Expers[1]]['meanYb_obs'],edgecolors='none',cmap=IRcmap,vmin=min_T,vmax=max_T,transform=ccrs.PlateCarree())
            elif i == 3:
                axs[i,j].set_extent([d03[DAtime][Expers[1]]['lon_min'],d03[DAtime][Expers[1]]['lon_max'],d03[DAtime][Expers[1]]['lat_min'],d03[DAtime][Expers[1]]['lat_max']], crs=ccrs.PlateCarree())
                axs[i,j].coastlines(resolution='10m', color='black',linewidth=0.5)
                axs[i,j].scatter(d_allIR[DAtime][Expers[1]]['lon_obs'],d_allIR[DAtime][Expers[1]]['lat_obs'],c_size,d_allIR[DAtime][Expers[1]]['meanYa_obs'],edgecolors='none',cmap=IRcmap,vmin=min_T,vmax=max_T,transform=ccrs.PlateCarree())


    # Colorbar
    caxes = fig.add_axes([0.13, 0.06, 0.75, 0.02])
    tb_range = np.linspace(min_T,max_T,8)
    cbar = fig.colorbar(cs,ticks=tb_range,orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=13)
    
  


    #subplot title
    font = {'size':12,}
    for DAtime in DAtimes:
        idx_t = DAtimes.index( DAtime )
        axs[0,idx_t].set_title(DAtime, font, fontweight='bold')

    #title for all
    fig.suptitle(Storm, fontsize=8, fontweight='bold')

    # Axis labels
    for j in range(3):
        for i in range(4):
            gl = axs[i,j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
            if i == 1 or i == 2:
                lon_min = d03[DAtime][Expers[0]]['lon_min']
                lon_max = d03[DAtime][Expers[0]]['lon_max']
                lat_min = d03[DAtime][Expers[0]]['lat_min']
                lat_max = d03[DAtime][Expers[0]]['lat_max']
            else:
                lon_min = d03[DAtime][Expers[1]]['lon_min']
                lon_max = d03[DAtime][Expers[1]]['lon_max']
                lat_min = d03[DAtime][Expers[1]]['lat_min']
                lat_max = d03[DAtime][Expers[1]]['lat_max']

            lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
            lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

            if j==0:
                gl.ylabels_left = True
                gl.ylabels_right = False
            else:
                gl.ylabels_left = False
                gl.ylabels_right = False

            if i == 3:
                gl.xlabels_top = False
                gl.xlabels_bottom = True
            else:
                gl.xlabels_top = False
                gl.xlabels_bottom = False

            gl.ylocator = mticker.FixedLocator(lat_ticks)
            gl.xlocator = mticker.FixedLocator(lon_ticks)
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 8}
            gl.ylabel_style = {'size': 10}

    des_name = 'IR_withoutObs.png'
    #des_name = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Expers[0]+'/Vis_analyze/Tb/IR/Obspace/'+DAtime+'_'+sensor+'_Obspace_expers.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    MP = ['WSM6','THO']
    DA = 'IR'
    
    sensor = 'abi_gr'
    ch_list = ['8',]

    # Time range set up
    start_time_str = '201709030500'
    end_time_str = '201709030700'
    Consecutive_times = True

    # -------------------------------------------------------    

    # Create experiment names
    Expers = []
    for imp in MP:
        Expers.append( UD.generate_one_name( Storm,DA,imp ))

    if not Consecutive_times:
        DAtimes = ['201709050000','201709050600','201709051200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    plot_Tb_withoutObs( )

















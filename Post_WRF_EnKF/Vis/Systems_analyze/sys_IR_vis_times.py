#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os
import glob
import numpy as np
import netCDF4 as nc
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as patches
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from global_land_mask import globe
import math
from datetime import datetime, timedelta
import time

import Util_data as UD
import Util_Vis

# Read variables at obs resolution/location for one experiment 
def read_Tbs_obsRes_oneExper(istorm,imp,ida,Exper_names,DAtimes,sensor):

    Hx_dir = big_dir+istorm+'/'+Exper_names[istorm][imp][ida]+'/Obs_Hx/IR/'
    dict_allTb = {}

    for DAtime in DAtimes[istorm]:

        Tb_file = Hx_dir+DAtime+'/mean_obs_res_d03_' + DAtime + '_' +  sensor + '.txt'
        lat_obs = []
        lon_obs = []
        Yo_obs = []
        # meanYb_obs = []
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
            #meanYb_obs.append( float(split_line[4]) )
            meanYa_obs.append( float(split_line[5]) )

        lat_obs = np.array( lat_obs )
        lon_obs = np.array( lon_obs )
        #ch_obs = np.array( ch_obs )
        Yo_obs = np.array( Yo_obs )
        #meanYb_obs = np.array( meanYb_obs )
        meanYa_obs = np.array( meanYa_obs )
    
        dict_allTb[DAtime] = {'lat_obs':lat_obs, 'lon_obs':lon_obs,'Yo_obs':Yo_obs, 'meanYa_obs':meanYa_obs}

    return dict_allTb

def plot_twoTimes_IR( Exper_Tb ):

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(6, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0.1}, linewidth=0.5, figsize=(9,15), dpi=400)

    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )

    for istorm in Storms:
        j = Storms.index( istorm )
        for i in range(6):
            # domain range
            if i <= 2:
                obs_range = Exper_Tb[istorm]['WSM6']['conv'][DAtimes[istorm][0]]
            else:
                obs_range = Exper_Tb[istorm]['WSM6']['conv'][DAtimes[istorm][1]]
            lat_min = np.min( obs_range['lat_obs'] )
            lat_max = np.max( obs_range['lat_obs'] )
            lon_min = np.min( obs_range['lon_obs'] )
            lon_max = np.max( obs_range['lon_obs'] )
            ax[i,j].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
            ax[i,j].coastlines(resolution='10m', color='black',linewidth=0.5)
            # row 0 - row 2
            conv_tb = Exper_Tb[istorm]['WSM6']['conv'][DAtimes[istorm][0]]
            ir_tb = Exper_Tb[istorm]['WSM6']['IR'][DAtimes[istorm][0]]
            ax[0,j].scatter(conv_tb['lon_obs'],conv_tb['lat_obs'],2.5,c=conv_tb['Yo_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
            ax[1,j].scatter(conv_tb['lon_obs'],conv_tb['lat_obs'],2.5,c=conv_tb['meanYa_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
            ax[2,j].scatter(ir_tb['lon_obs'],ir_tb['lat_obs'],2.5,c=ir_tb['meanYa_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
            # row 3 - row 5
            conv_tb = Exper_Tb[istorm]['WSM6']['conv'][DAtimes[istorm][1]]
            ir_tb = Exper_Tb[istorm]['WSM6']['IR'][DAtimes[istorm][1]]
            ax[3,j].scatter(conv_tb['lon_obs'],conv_tb['lat_obs'],2.5,c=conv_tb['Yo_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
            ax[4,j].scatter(conv_tb['lon_obs'],conv_tb['lat_obs'],2.5,c=conv_tb['meanYa_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
            cs = ax[5,j].scatter(ir_tb['lon_obs'],ir_tb['lat_obs'],2.5,c=ir_tb['meanYa_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
            # grid lines
            lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
            lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2)) 
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
            gl.xlabels_top = False
            gl.xlabels_bottom = True
            gl.ylabels_left = True
            gl.ylabels_right = False
            gl.ylocator = mticker.FixedLocator(lat_ticks)
            gl.xlocator = mticker.FixedLocator(lon_ticks)
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 10.5}
            gl.ylabel_style = {'size': 10}

    # Colorbar
    cbaxes = fig.add_axes([0.90, 0.1, 0.03, 0.78])
    cbar = fig.colorbar(cs,cax=cbaxes,fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=18)

    # subplot title
    for istorm in Storms:
        j = Storms.index( istorm )
        ax[0,j].set_title(istorm, fontsize='20', fontweight='bold')

    des_name = plot_dir+'twoTimes.png'
    plt.savefig(des_name,dpi=300)
    print('Saving the figure: ',des_name)


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['JOSE','MARIA','IRMA']#['HARVEY','IRMA','JOSE','MARIA']
    DA = ['conv','IR']
    MP = ['WSM6',]
    sensor = 'abi_gr'

    DAtimes = {'IRMA':['201709030000','201709040000'],'JOSE':['201709050000','201709060000'],'MARIA':['201709160000','201709170000']}
    if_plot = True

    #-----------------------------------

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = {}
        for imp in MP:
            Exper_names[istorm][imp] = {}
            for ida in DA:
                Exper_names[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # Read obs, Hxb, and Hxa of all files
    Exper_Tb = {}
    for istorm in Storms:
        Exper_Tb[istorm] = {}
        for imp in MP:
            Exper_Tb[istorm][imp] = {}
            for ida in DA:
                iExper = Exper_names[istorm][imp][ida]
                if iExper is not None:
                    Exper_Tb[istorm][imp][ida] = read_Tbs_obsRes_oneExper(istorm,imp,ida,Exper_names,DAtimes,sensor)
                else:
                    Exper_Tb[istorm][imp][ida] = None

    # Plot Tbs
    if if_plot:
        # Create plot dir
        plot_dir = small_dir+'/SYSTEMS/Vis_analyze/Tb/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
        # plot 
        plot_twoTimes_IR( Exper_Tb )










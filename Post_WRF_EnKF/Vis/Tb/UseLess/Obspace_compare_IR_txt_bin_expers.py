#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os
import glob
import numpy as np
import Util_Vis
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from global_land_mask import globe
import math
from datetime import datetime, timedelta
import matlab.engine
import time


# Plotting domain
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

# Storm center produced by ATCF of the NWP
def read_TCvitals(tc_file, DAtime):

    with open(tc_file) as tmp:
        tc_all = tmp.readlines()

    tc_lat = []
    tc_lon = []
    for line in tc_all:
        line_split = line.split()
        tc_time = line_split[3]+line_split[4]

        if tc_time == DAtime:
            print('Time from TCvitals:', tc_time)
            # Read latitude
            if 'N' in line_split[5]:
                tc_lat.append(float(line_split[5].replace('N',''))/10)
            else:
                tc_lat.append( 0-float(line_split[5].replace('S',''))/10)
            # Read longitude
            if 'W' in line_split[6]:
                tc_lon.append(0-float(line_split[6].replace('W',''))/10)
            else:
                tc_lon.append(float(line_split[6].replace('E',''))/10)

            break

    return tc_lon, tc_lat



def read_allTb( Expers,wrf_dirs, DAtime, sensor ):

    d_IR_all = {}
    for iExper in Expers:
        idx_exper = Expers.index( iExper )
        Tb_file = wrf_dirs[idx_exper]+"/mean_obs_res_d03" + DAtime + '_' +  sensor + '.txt' 

        lat_obs = []
        lon_obs = []
        ch_obs = []
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
            ch_obs.append( int(split_line[2]) )
            Yo_obs.append( float(split_line[3]) )
            meanYb_obs.append( float(split_line[4]) )
            meanYa_obs.append( float(split_line[5]) )
    
        lat_obs = np.array( lat_obs )
        lon_obs = np.array( lon_obs )
        ch_obs = np.array( ch_obs )
        Yo_obs = np.array( Yo_obs )
        meanYb_obs = np.array( meanYb_obs )
        meanYa_obs = np.array( meanYa_obs )
        print('Number of NaN in meanYa_obs', sum(np.isnan(meanYa_obs)))
        d_Tb =  {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'ch_obs':ch_obs, 'Yo_obs':Yo_obs, 'meanYb_obs':meanYb_obs, 'meanYa_obs':meanYa_obs}
        d_IR_all[iExper] = d_Tb

    return d_IR_all


def plot_Tb(Storm, Expers, DAtime, sensor, ch_list, big_dir ):

    # Read WRF domain
    wrf_file = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Expers[0]+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = read_wrf_domain( wrf_file )

    # Read Tbs of obs, Hxb, Hxa
    wrf_dirs = []
    for iExper in Expers:
        wrf_dirs.append( big_dir+Storm+'/'+iExper+'/Obs_Hx/IR/'+DAtime )
    d_all = read_allTb( Expers,wrf_dirs, DAtime, sensor )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = read_TCvitals('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # ------------------ Plot -----------------------
    fig, axs=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,2.5), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )

    # Exper 1
    axs.flat[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    axs.flat[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    axs.flat[0].scatter(d_all[Expers[0]]['lon_obs'],d_all[Expers[0]]['lat_obs'],1.5,c=d_all[Expers[0]]['meanYb_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        axs.flat[0].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())

    # Exper 2
    axs.flat[1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    axs.flat[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    xb_Tb = axs.flat[1].scatter(d_all[Expers[1]]['lon_obs'], d_all[Expers[1]]['lat_obs'],1.5,c=d_all[Expers[1]]['meanYb_obs'],\
                edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        axs.flat[1].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.12, 0.1, 0.45, 0.02])
    xwv_bar = fig.colorbar(xb_Tb,ax=axs[0:2],orientation="horizontal", cax=caxes)
    xwv_bar.ax.tick_params()

    # Exper 1 - Exper 2
    min_incre = -20
    max_incre = 20
    axs.flat[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    axs.flat[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = axs.flat[2].scatter(d_all[Expers[0]]['lon_obs'], d_all[Expers[0]]['lat_obs'],1.5,c=d_all[Expers[0]]['meanYb_obs']-d_all[Expers[1]]['meanYb_obs'],\
                edgecolors='none', cmap='RdBu_r', vmin=min_incre, vmax=max_incre, transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        axs.flat[2].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.65, 0.1, 0.25, 0.02])
    cb_diff_ticks = np.linspace(min_incre, max_incre, 5, endpoint=True)
    cbar = fig.colorbar(cs, ax=axs[2:], ticks=cb_diff_ticks, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params()

    #subplot title
    font = {'size':8,}
    axs.flat[0].set_title('Yo', font, fontweight='bold')
    axs.flat[1].set_title('H(Xb)', font, fontweight='bold')
    axs.flat[2].set_title('H(Xa)', font, fontweight='bold')

    #title for all
    fig.suptitle(Storm+': '+DAtime, fontsize=8, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

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
        gl.xlabel_style = {'size': 4}
        gl.ylabel_style = {'size': 6}

    des_name = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Expers[0]+'/Vis_analyze/Tb/IR/Obspace/'+DAtime+'_'+sensor+'_Obspace_expers.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)






if __name__ == '__main__':

    Storm = 'JOSE'
    Expers = ['J_DA+J_WRF+J_init-SP-intel17-THO-30hr-hroi900','IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900',]
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    sensor = 'abi_gr'
    ch_list = ['8',]
    Plot_IR = True

    # Time range set up
    start_time_str = '201709160000'
    end_time_str = '201709160000'
    Consecutive_times = False

    if not Consecutive_times:
        DAtimes = ['201709050000','201709050600','201709051200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Plot precipitable water (integral column of water vapor)
    if Plot_IR:
        start_time=time.process_time()
        for DAtime in DAtimes:
            plot_Tb( Storm, Expers, DAtime,  sensor, ch_list, big_dir)
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')



import os
import glob
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as patches
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math
from datetime import datetime, timedelta
import time

import Util_data as UD
import Util_Vis


# Read observed IR obs
def readIR_obs( Tb_file ):

    lat_obs = []
    lon_obs = []
    ch_obs = []
    Yo_obs = []

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

    #if np.size(lat_obs) != dict_ss_len[sensor]:
    #    raise ValueError('The length of post-processed file is not equal to the pre-processed file!')

    lat_obs = np.array( lat_obs )
    lon_obs = np.array( lon_obs )
    ch_obs = np.array( ch_obs )
    Yo_obs = np.array( Yo_obs )

    lat_min = np.min( lat_obs )
    lat_max = np.max( lat_obs )
    lon_min = np.min( lon_obs)
    lon_max = np.max( lon_obs )

    d_IRobs = {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'ch_obs':ch_obs, 'Yo_obs':Yo_obs,'lat_min':lat_min,'lat_max':lat_max,'lon_min':lon_min,'lon_max':lon_max}
    return d_IRobs

# Read simulated IR at obs resolution/location
def read_simuIR_obsLoc( Tb_file ):

    lat_obs = []
    lon_obs = []
    #ch_obs = []
    #Yo_obs = []
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
        #ch_obs.append( int(split_line[2]) )
        #Yo_obs.append( float(split_line[3]) )
        meanYb_obs.append( float(split_line[4]) )
        meanYa_obs.append( float(split_line[5]) )

    lat_obs = np.array( lat_obs )
    lon_obs = np.array( lon_obs )
    #ch_obs = np.array( ch_obs )
    #Yo_obs = np.array( Yo_obs )
    meanYb_obs = np.array( meanYb_obs )
    meanYa_obs = np.array( meanYa_obs )
    print('Number of NaN in meanYa_obs', sum(np.isnan(meanYa_obs)))

    d_simu = {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'meanYb_obs':meanYb_obs, 'meanYa_obs':meanYa_obs}
    return d_simu

#Layout:
#Columns: storms
#Rows: obs,WSM6_CONV,WSM6_IR,WSM6_IRMW,THO_CONV,THO_IR,THO_IRMW
def plot_snapshot():

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200) # standard: 6.5,8.5
    grids = fig.add_gridspec(ncols=4,nrows=7,hspace=0.0,wspace=0.0,top=0.93,left=0.12,)
    ax = {}
    ax['obs'] = {}
    for ida in DA:
        ax['WSM6_'+ida] = {}
        ax['THO_'+ida] = {}
    # obs
    for ist in Storms:
        ax['obs'][ist] = fig.add_subplot( grids[0,Storms.index(ist)],projection=ccrs.PlateCarree())
    # WSM6
    for ida in DA:
        for ist in Storms:
            ax['WSM6_'+ida][ist] = fig.add_subplot( grids[1+DA.index(ida),Storms.index(ist)],projection=ccrs.PlateCarree())
    # THO
    for ida in DA:
        for ist in Storms:
            ax['THO_'+ida][ist] = fig.add_subplot( grids[4+DA.index(ida),Storms.index(ist)],projection=ccrs.PlateCarree() )

    # Customization
    min_tb = 185
    max_tb = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )

    # Set map
    row = ['obs',]
    for ida in DA:
        row.append('WSM6_'+ida)
    for ida in DA:
        row.append('THO_'+ida)
    for ir in row:
        for ist in Storms:
            ax[ir][ist].set_extent([IR_obs[ist]['lon_min'],IR_obs[ist]['lon_max'],IR_obs[ist]['lat_min'],IR_obs[ist]['lat_max']], crs=ccrs.PlateCarree())
            ax[ir][ist].coastlines(resolution='10m', color='black',linewidth=0.5)

    # plot obs
    for ist in Storms:
        obs_s = ax['obs'][ist].scatter(IR_obs[ist]['lon_obs'],IR_obs[ist]['lat_obs'],1.5,c=IR_obs[ist]['Yo_obs'],edgecolors='none', cmap=IRcmap, vmin=min_tb, vmax=max_tb,transform=ccrs.PlateCarree())
    # plot analyzed Tb for WSM6
    for ida in DA:
        for ist in Storms:
            ax['WSM6_'+ida][ist].scatter(IR_simu[ist]['WSM6'][ida]['lon_obs'],IR_simu[ist]['WSM6'][ida]['lat_obs'],1.5,c=IR_simu[ist]['WSM6'][ida]['meanYa_obs'],edgecolors='none', cmap=IRcmap, vmin=min_tb, vmax=max_tb,transform=ccrs.PlateCarree())
    # plot analyzed Tb for THO
    for ida in DA:
        for ist in Storms:
            ax['THO_'+ida][ist].scatter(IR_simu[ist]['THO'][ida]['lon_obs'],IR_simu[ist]['THO'][ida]['lat_obs'],1.5,c=IR_simu[ist]['THO'][ida]['meanYa_obs'],edgecolors='none', cmap=IRcmap, vmin=min_tb, vmax=max_tb,transform=ccrs.PlateCarree())

    # colorbar
    cbar_ax = fig.add_axes([0.10, 0.06, 0.8, 0.015])
    cbar = fig.colorbar(obs_s, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('GOES-16 Ch8 Tb (K)')

    # labels
    row = ['obs',]
    for ida in DA:
        row.append('WSM6_'+ida)
    for ida in DA:
        row.append('THO_'+ida)

    # Add storm information
    for ist in Storms:
        if Storms.index(ist) == 0:
            fig.text(0.21,0.95,ist, fontsize=12, ha='center', va='center')
        elif Storms.index(ist) == 1:
            fig.text(0.42,0.95,ist, fontsize=12, ha='center', va='center')
        elif Storms.index(ist) == 2:
            fig.text(0.61,0.95,ist, fontsize=12, ha='center', va='center')
        else:
            fig.text(0.8,0.95,ist, fontsize=12, ha='center', va='center')
    # Add obs and experiment
    top_h = [0.88 - it for it in [0,0.12,0.24,0.35,0.48,0.58,0.71]]
    for ir in row:
        if ir == 'obs':
            ir_n = 'OBS'
        else:
            ir_n = ir
        fig.text(0.07,top_h[row.index(ir)],ir_n,fontsize=8,ha='center',va='center',rotation='vertical')

    # axis
    for ist in Storms:
        lon_ticks = list(range(math.ceil(IR_obs[ist]['lon_min'])-2, math.ceil(IR_obs[ist]['lon_max'])+2,2))
        lat_ticks = list(range(math.ceil(IR_obs[ist]['lat_min'])-2, math.ceil(IR_obs[ist]['lat_max'])+2,2))
        for ir in row:
            gl = ax[ir][ist].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5,alpha=0.7,color='gray',linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.left_labels = True
            if ir == row[-1]:
                gl.bottom_labels = True
            else:
                gl.bottom_labels = False
                
            gl.ylocator = mticker.FixedLocator(lat_ticks)
            gl.xlocator = mticker.FixedLocator(lon_ticks)
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 4}
            gl.ylabel_style = {'size': 5}


    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/snapshot_IR_Cycle1.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY','JOSE','MARIA','IRMA']
    MP = ['WSM6','THO']
    DA = ['CONV','IR','IR+MW']
    # sensor
    sensor = 'abi_gr'
    ch_list = ['8',]
    # time
    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    t_incre = 1 # 1 hour
    #------------------------------------
    wrf_dir = big_dir

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = {}
        for imp in MP:
            Exper_names[istorm][imp] = {}
            for ida in DA:
                Exper_names[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # DA time of interest
    DAtime = {}
    for ist in Storms:
        tmp_time = datetime.strptime(start_time_str[ist],"%Y%m%d%H%M") + timedelta(hours=t_incre)
        DAtime[ist] = tmp_time.strftime("%Y%m%d%H%M") 

    # Read observed IR radiances
    IR_obs = {}
    for ist in Storms:
        Hx_file = big_dir+ist+'/'+Exper_names[ist]['WSM6']['IR']+'/Obs_Hx/IR/'+DAtime[ist]+'/mean_obs_res_d03_'+DAtime[ist]+'_'+sensor+'.txt'
        IR_obs[ist] = readIR_obs( Hx_file )  

    # Read simulated IR radiances
    IR_simu = {}
    for ist in Storms:
        IR_simu[ist] = {}
        for imp in MP:
            IR_simu[ist][imp] = {}
            for ida in DA:
                Hx_file = big_dir+ist+'/'+Exper_names[ist][imp][ida]+'/Obs_Hx/IR/'+DAtime[ist]+'/mean_obs_res_d03_'+DAtime[ist]+'_'+sensor+'.txt'
                IR_simu[ist][imp][ida] = read_simuIR_obsLoc( Hx_file ) 

    plot_snapshot()























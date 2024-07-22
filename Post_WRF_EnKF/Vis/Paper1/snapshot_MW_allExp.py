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
import pandas as pd

import Util_data as UD
import Util_Vis

# Read observed MW obs
def readMW_obs( ist ):
    
    # Read the file into a pandas DataFrame
    Hx_dir = big_dir+ist+'/'+Exper_names[ist][MP[0]]['IR+MW']+'/Obs_Hx/MW/'+DAtime[ist]+'/'
    file_name = 'mean_obs_res_d03_'+DAtime[ist]+'.tb.'+d_highf[ist]['sensor']+'.crtm.conv.txt'
    # Read the file into a pandas DataFrame
    df = pd.read_csv( Hx_dir+file_name,delim_whitespace=True,skiprows=1,header=None,dtype=str )
    ch_obs = np.array( [int(it) for it in df[2].values] )
    lat_obs = np.array( [float(it) for it in df[0].values] )
    lon_obs = np.array( [float(it) for it in df[1].values] )
    mw_obs = np.array( [float(it) for it in df[3].values] )
    # Only take rows that are associated with the specified channel
    lat_obs = lat_obs[ch_obs == d_highf[ist]['ch']]
    lon_obs = lon_obs[ch_obs == d_highf[ist]['ch']]
    mw_obs = mw_obs[ch_obs == d_highf[ist]['ch']]

    lat_min = np.min( lat_obs )
    lat_max = np.max( lat_obs )
    lon_min = np.min( lon_obs)
    lon_max = np.max( lon_obs )

    d_MWobs = {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'mw_obs':mw_obs,'lat_min':lat_min,'lat_max':lat_max,'lon_min':lon_min,'lon_max':lon_max}
    return d_MWobs

# Read simulated IR at obs resolution/location
def read_simuMW_obsLoc( ist,imp,ida ):

    # Read the file into a pandas DataFrame
    Hx_dir = big_dir+ist+'/'+Exper_names[ist][imp][ida]+'/Obs_Hx/MW/'+DAtime[ist]+'/'
    file_name = 'mean_obs_res_d03_'+DAtime[ist]+'.tb.'+d_highf[ist]['sensor']+'.crtm.conv.txt'
    # Read the file into a pandas DataFrame
    df = pd.read_csv( Hx_dir+file_name,delim_whitespace=True,skiprows=1,header=None,dtype=str )
    ch_obs = np.array( [int(it) for it in df[2].values] )
    lat_obs = np.array( [float(it) for it in df[0].values] )
    lon_obs = np.array( [float(it) for it in df[1].values] )
    mw_simu = np.array( [float(it) for it in df[5].values] )
    # Only take rows that are associated with the specified channel
    lat_obs = lat_obs[ch_obs == d_highf[ist]['ch']]
    lon_obs = lon_obs[ch_obs == d_highf[ist]['ch']]
    mw_simu = mw_simu[ch_obs == d_highf[ist]['ch']]

    d_MWsimu = {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'meanYa_obs':mw_simu}
    return d_MWsimu

#Layout:
#Columns: storms
#Rows: obs,WSM6_CONV,WSM6_IR,WSM6_IRMW,THO_CONV,THO_IR,THO_IRMW
def plot_snapshot():

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=300) # standard: 6.5,8.5
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
    max_tb=300
    min_tb=80
    min_Jet=150
    MWJet = Util_Vis.newJet(300,80,150)

    # Set map
    row = ['obs',]
    for ida in DA:
        row.append('WSM6_'+ida)
    for ida in DA:
        row.append('THO_'+ida)
    for ir in row:
        for ist in Storms:
            ax[ir][ist].set_extent([MW_obs[ist]['lon_min'],MW_obs[ist]['lon_max'],MW_obs[ist]['lat_min'],MW_obs[ist]['lat_max']], crs=ccrs.PlateCarree())
            ax[ir][ist].coastlines(resolution='10m', color='black',linewidth=0.5)

    # plot obs
    for ist in Storms:
        obs_s = ax['obs'][ist].scatter(MW_obs[ist]['lon_obs'],MW_obs[ist]['lat_obs'],1.5,c=MW_obs[ist]['mw_obs'],edgecolors='none', cmap=MWJet, vmin=min_tb, vmax=max_tb,transform=ccrs.PlateCarree())
    # plot analyzed Tb for WSM6
    for ida in DA:
        for ist in Storms:
            ax['WSM6_'+ida][ist].scatter(MW_simu[ist]['WSM6'][ida]['lon_obs'],MW_simu[ist]['WSM6'][ida]['lat_obs'],1.5,c=MW_simu[ist]['WSM6'][ida]['meanYa_obs'],edgecolors='none', cmap=MWJet, vmin=min_tb, vmax=max_tb,transform=ccrs.PlateCarree())
    # plot analyzed Tb for THO
    for ida in DA:
        for ist in Storms:
            ax['THO_'+ida][ist].scatter(MW_simu[ist]['THO'][ida]['lon_obs'],MW_simu[ist]['THO'][ida]['lat_obs'],1.5,c=MW_simu[ist]['THO'][ida]['meanYa_obs'],edgecolors='none', cmap=MWJet, vmin=min_tb, vmax=max_tb,transform=ccrs.PlateCarree())

    # colorbar
    cbar_ax = fig.add_axes([0.10, 0.06, 0.78, 0.015])
    cbar = fig.colorbar(obs_s, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('GPM MW Tb (K)')

    # labels
    row = ['obs',]
    for ida in DA:
        row.append('WSM6_'+ida)
    for ida in DA:
        row.append('THO_'+ida)

    # Add storm information
    for ist in Storms:
        if Storms.index(ist) == 0:
            fig.text(0.21,0.98,ist, fontsize=12, ha='center', va='center')
            fig.text(0.21,0.96,DAtime[ist], fontsize=10, ha='center', va='center')
            fig.text(0.21,0.94,d_highf[ist]['sensor'], fontsize=10, ha='center', va='center')
        elif Storms.index(ist) == 1:
            fig.text(0.42,0.98,ist, fontsize=12, ha='center', va='center')
            fig.text(0.42,0.96,DAtime[ist], fontsize=10, ha='center', va='center')
            fig.text(0.42,0.94,d_highf[ist]['sensor'], fontsize=10, ha='center', va='center')
        elif Storms.index(ist) == 2:
            fig.text(0.61,0.98,ist, fontsize=12, ha='center', va='center')
            fig.text(0.61,0.96,DAtime[ist], fontsize=10, ha='center', va='center')
            fig.text(0.61,0.94,d_highf[ist]['sensor'], fontsize=10, ha='center', va='center')
        else:
            fig.text(0.8,0.98,ist, fontsize=12, ha='center', va='center')
            fig.text(0.8,0.96,DAtime[ist], fontsize=10, ha='center', va='center')
            fig.text(0.8,0.94,d_highf[ist]['sensor'], fontsize=10, ha='center', va='center')
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
        lon_ticks = list(range(math.ceil(MW_obs[ist]['lon_min'])-2, math.ceil(MW_obs[ist]['lon_max'])+2,2))
        lat_ticks = list(range(math.ceil(MW_obs[ist]['lat_min'])-2, math.ceil(MW_obs[ist]['lat_max'])+2,2))
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
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/snapshot_MW.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    MP = ['WSM6','THO']
    DA = ['CONV','IR','IR+MW']
    # sensor and channel
    d_highf = {'HARVEY':{'sensor':'ssmis_f18','ch':9},
                'IRMA':{'sensor':'mhs_metop-b','ch':5},
                'JOSE':{'sensor':'ssmis_f16','ch':9},
                'MARIA':{'sensor':'saphir_meghat','ch':5}}
    # time
    DAtime = {'HARVEY':'201708221300','IRMA':'201709030100','JOSE':'201709050600','MARIA':'201709160200'}

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

    # Read observed MW radiances
    MW_obs = {}
    for ist in Storms:
        MW_obs[ist] = readMW_obs( ist )


    # Read simulated IR radiances
    MW_simu = {}
    for ist in Storms:
        MW_simu[ist] = {}
        for imp in MP:
            MW_simu[ist][imp] = {}
            for ida in DA:
                MW_simu[ist][imp][ida] = read_simuMW_obsLoc( ist,imp,ida )

    plot_snapshot()












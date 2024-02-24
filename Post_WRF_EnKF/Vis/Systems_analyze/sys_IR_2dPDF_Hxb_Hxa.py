import os
import glob
import numba
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import math
from datetime import datetime, timedelta
import time
from fast_histogram import histogram2d as hist2d
import time
import matplotlib.colors as mcolors

import Util_data as UD
import Util_Vis
import Diagnostics as Diag

# Read variables at obs resolution/location for one experiment 
def read_Tbs_obsRes_oneExper(istorm,imp,ida,Exper_names,d_times,sensor):

    Hx_dir = big_dir+istorm+'/'+Exper_names[istorm]+'/Obs_Hx/IR/'
    dict_allTb = {}

    Yo_all = []
    meanYb_all = []
    meanYa_all = []

    for DAtime in d_times[istorm]:

        Tb_file = Hx_dir+DAtime+'/mean_obs_res_d03_' + DAtime + '_' +  sensor + '.txt'
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
            Yo_obs.append( float(split_line[3]) )
            meanYb_obs.append( float(split_line[4]) )
            meanYa_obs.append( float(split_line[5]) )

        # add up records
        Yo_all.extend( Yo_obs )
        meanYb_all.extend( meanYb_obs )
        meanYa_all.extend( meanYa_obs )

    dict_allTb['All_times'] = {'Yo_obs':Yo_all, 'meanYb_obs':meanYb_all, 'meanYa_obs':meanYa_all}

    return dict_allTb

def combine_storms_allTimes( Exper_Tb ):

    allStorm_tb = {}
    all_Yo = []
    all_meanYb = []
    all_meanYa = []

    for istorm in Storms:
        tb_case = Exper_Tb[istorm]['All_times']
        all_Yo.extend( tb_case['Yo_obs'] )
        all_meanYb.extend( tb_case['meanYb_obs'] )
        all_meanYa.extend( tb_case['meanYa_obs'] )

    all_Yo = np.array( all_Yo )
    all_meanYb = np.array( all_meanYb )
    all_meanYa = np.array( all_meanYa )
    
    allStorm_tb = {'Yo_obs':all_Yo, 'meanYb_obs':all_meanYb, 'meanYa_obs':all_meanYa}

    return allStorm_tb

# Plot 2D histogram of simulated IR and Obs
def Plot_2Dhist( d_hcount ):

    # Set up figure
    fig,axs = plt.subplots(1, 2, figsize=(10,6), dpi=300 )

    # Verification using scatter function
    #axs.flat[0].scatter(Storms_Tb['Yo_obs'],Storms_Tb['meanYb_obs'],s=2)
    #axs.flat[1].scatter(Storms_Tb['Yo_obs'],Storms_Tb['meanYa_obs'],s=2)

    # Customize the colormap
    color_intervals = [0,0.5,1,1.5,2,2.5,3,3.5,4.0]
    exist_cmap = plt.cm.jet.reversed()
    colors = exist_cmap(np.linspace(0,1,len(color_intervals)-1))
    new_map = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals)-1)

    # Plot 2d histogram
    axs.flat[0].imshow(np.transpose(d_hcount['meanYb_Yo']),cmap=new_map,aspect='equal',
                            vmin=0,vmax=4,extent=None)
    show = axs.flat[1].imshow(np.transpose(d_hcount['meanYa_Yo']),cmap=new_map,aspect='equal',
                            vmin=0,vmax=4,extent=None)
                            
    # Add color bar below the plot
    caxes = fig.add_axes([0.12, 0.88, 0.78, 0.03])
    color_bar = fig.colorbar(show,cax=caxes,orientation='horizontal')#ticks=bounds)
    color_bar.ax.xaxis.set_major_locator(plt.FixedLocator([0,1,2,3,4]))
    #color_bar.ax.set_xticklabels( bounds_str, rotation=45)
    color_bar.ax.tick_params(labelsize=15)
    #color_bar.ax.set_xlabel('Domain-mean Increment',fontsize=15)

    # Plot lines
    x = np.linspace(0,number_bins) 
    for j in range(2):
        axs.flat[j].plot(x, x, color='grey', linestyle='dashed', linewidth=2)
        axs.flat[j].plot(x, x-10*scale_factor, color='grey', linestyle='dashed', linewidth=2)
        axs.flat[j].plot(x, x+10*scale_factor, color='grey', linestyle='dashed', linewidth=2)
        axs.flat[j].plot(x, x-20*scale_factor, color='grey', linestyle='dashed', linewidth=2)
        axs.flat[j].plot(x, x+20*scale_factor, color='grey', linestyle='dashed', linewidth=2)

    # Plot a patch
    cloud_th = 205
    cloud_th_point = (cloud_th-min_Tb_rg)*scale_factor
    for j in range(2):
        axs.flat[j].fill_betweenx(x, 0, cloud_th_point,alpha=0.2,facecolor='grey')

    # ticks and labels
    for j in range(2):
        tick = list(range(0,number_bins+1,10*scale_factor))
        tick_label = list(range(0,max_Tb_rg-min_Tb_rg+1,10))
        axs.flat[j].set_xticks( tick )
        axs.flat[j].set_yticks( tick )
        axs.flat[j].set_xticklabels( [str(it+min_Tb_rg) for it in tick_label],fontsize=15 )
        axs.flat[j].set_yticklabels( [str(it+min_Tb_rg) for it in tick_label],fontsize=15 )
        axs.flat[j].set_xlim(xmin=0,xmax=number_bins)
        axs.flat[j].set_ylim(ymin=0,ymax=number_bins)
        axs.flat[j].set_xlabel('Observation (K)',fontsize=15)
    axs.flat[0].set_ylabel('Background (K)',fontsize=15)
    axs.flat[1].yaxis.set_label_position("right")
    axs.flat[1].set_ylabel('Analysis (K)',fontsize=15)

    # set title
    fig.suptitle( 'Density scatterplots of IR Tbs over '+str(len(dict_times[Storms[0]]))+' Cycles', fontsize=15, fontweight='bold')

    # Save
    des_name = small_dir+'SYSTEMS/Vis_analyze/Tb/IR_2dPDF_'+str(len(dict_times[Storms[0]]))+'cycles.png'
    plt.savefig( des_name )
    print( 'Saving the figure: ', des_name )
    plt.close()

    return None



if __name__ == '__main__':

    big_dir = '/scratch_S2/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['IRMA','JOSE','MARIA','HARVEY']#['HARVEY','IRMA','JOSE','MARIA']
    DA = 'IR'
    MP = 'WSM6'
    sensor = 'abi_gr'

    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708241200','IRMA':'201709050000','JOSE':'201709070000','MARIA':'201709180000'}
    Consecutive_times = True

    # Range
    min_Tb_rg = 190
    max_Tb_rg = 260

    # number of bins
    scale_factor = 2 # 1: 1k per bin; 2: 0.5k per bin
    number_bins = (max_Tb_rg-min_Tb_rg)*scale_factor

    If_plot = True
    #------------------------------------

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = UD.generate_one_name( istorm,DA,MP )

    # Identify DA times in the period of interest
    dict_times = UD.generate_times( Storms, start_time_str, end_time_str )

    # Read obs, Hxb, and Hxa of all files
    Exper_Tb = {}
    for istorm in Storms:
        Exper_Tb[istorm] = read_Tbs_obsRes_oneExper(istorm,MP,DA,Exper_names,dict_times,sensor)

    # Combine data for all storms
    Storms_Tb = combine_storms_allTimes( Exper_Tb )

    # Make bins
    dcount = {}
    meanYb_Yo = hist2d(Storms_Tb['Yo_obs'],Storms_Tb['meanYb_obs'],range=[[min_Tb_rg,max_Tb_rg],[min_Tb_rg,max_Tb_rg]],bins=number_bins)
    meanYa_Yo = hist2d(Storms_Tb['Yo_obs'],Storms_Tb['meanYa_obs'],range=[[min_Tb_rg,max_Tb_rg],[min_Tb_rg,max_Tb_rg]],bins=number_bins)
    # Compute log10 of non-zero elements, return NAN where elements are ZERO
    with np.errstate(divide='ignore', invalid='ignore'):
        dcount['meanYb_Yo'] = np.where(meanYb_Yo != 0, np.log10(meanYb_Yo), np.nan)
        dcount['meanYa_Yo'] = np.where(meanYa_Yo != 0, np.log10(meanYa_Yo), np.nan)

    # Plot
    if If_plot:
        Plot_2Dhist( dcount )#Storms_Tb )

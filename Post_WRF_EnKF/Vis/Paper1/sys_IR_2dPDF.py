import os
import glob
import numba
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib import ticker
import math
from datetime import datetime, timedelta
import time
from fast_histogram import histogram2d as hist2d
import time
import matplotlib.colors as mcolors
from scipy import interpolate
import matplotlib.path as mpath
import matplotlib.patches as mpatches

import Util_data as UD
import Util_Vis
import IR_2dPDF as IRpdf

# Generate time series
def generate_times( Storms, start_time_str, end_time_str, interval ):

    dict_times = {}
    for ist in Storms:
        time_diff = datetime.strptime(end_time_str[ist],"%Y%m%d%H%M") - datetime.strptime(start_time_str[ist],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[ist],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+interval, interval))]
        dict_times[ist] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times

# Read storm HPI info from pre-calculated files
def Read_storm_HPI( mfile ):

    d_model = {}
    print('Reading ', mfile)
    with open(mfile) as f:
        next(f)
        all_lines = f.readlines()

    for line in all_lines:
        if 'Success' in line:
            break
        split_line = line.split()
        itime = split_line[0] 
        lon_sc = float(split_line[1])
        lat_sc = float(split_line[2]) 
        d_model[itime] = (lon_sc,lat_sc)

    return d_model

# Read Tbs and subset them if necessary
def Read_Subset_Tbs( istorm,exp_name,DAtimes ):

    # Calculate the modelled SC if needed
    if Subset_data:
        # Find the storm center
        if model_center:
            if istorm == 'HARVEY':
                st_time = '201708221200'
                end_time = '201708231200'
            elif istorm == 'IRMA':
                st_time = '201709030000'
                end_time = '201709040000'
            elif istorm == 'JOSE':
                st_time = '201709050000'
                end_time = '201709060000'
            elif istorm == 'MARIA':
                st_time = '201709160000'
                end_time = '201709170000'

            mfile = small_dir+istorm+'/'+exp_name+'/Data_analyze/'+'HPI_wrf_enkf_input_d03_mean.'+st_time+'_'+end_time+'.txt'
            d_model = Read_storm_HPI( mfile )
        else: # best track
            # Read the best track position every 6 hours
            d_btk6Hrs = UD.read_bestrack(Storm)
            # Read the hourly best track position 
            d_btkHr = UD.interpolate_locations( DAtimes, d_btk6Hrs) # time,lat,lon

    # Read IR Tbs
    Hx_dir = big_dir+istorm+'/'+exp_name+'/Obs_Hx/IR/'
    
    d_allTb = {}
    Yo_all = []
    meanYb_all = []
    meanYa_all = []

    for DAtime in DAtimes:

        t_idx = DAtimes.index( DAtime )
        Tb_file = Hx_dir+DAtime+'/mean_obs_res_d03_' + DAtime + '_' +  sensor + '.txt'
        Yo_obs = []
        meanYb_obs = []
        meanYa_obs = []

        # Read records for file at DAtime
        print('Reading ', Tb_file)
        with open(Tb_file) as f:
            next(f)
            all_lines = f.readlines()

        for line in all_lines:
            split_line = line.split()
            if Subset_data:
                # determine if this obs is qualified
                lat_obs = float(split_line[0])
                lon_obs = float(split_line[1])
                sc_lon = d_model[DAtime][0]
                sc_lat = d_model[DAtime][1]
                dist = UD.mercator_distance(lon_obs, lat_obs, sc_lon, sc_lat)
                if dist <= max_radius:
                    Yo_obs.append( float(split_line[3]) )
                    meanYb_obs.append( float(split_line[4]) )
                    meanYa_obs.append( float(split_line[5]) )
                else:
                    continue
            else:
                Yo_obs.append( float(split_line[3]) )
                meanYb_obs.append( float(split_line[4]) )
                meanYa_obs.append( float(split_line[5]) )

        # Add up records
        Yo_all.extend( Yo_obs )
        meanYb_all.extend( meanYb_obs )
        meanYa_all.extend( meanYa_obs )

    d_allTb['All_times'] = {'Yo_obs':Yo_all, 'meanYb_obs':meanYb_all, 'meanYa_obs':meanYa_all}

    return d_allTb


# Plot BMO V.S. either background or Obs
def plt_BMO( d_hcount ):
    # Set up figure
    fig = plt.figure( figsize=(6.5,8),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=3,top=0.95,bottom=0.08,right=0.90,hspace=0.05)
    ax = {}
    row = ['bmo_Yo_conv','bmo_Yo_ir','amo_Yo_ir']
    for ir in range(len(row)):
        ax[row[ir]] = {}
        for imp in MP:
            ist_grids = outer_grids[ir].subgridspec(1, 2, wspace=0.05)
            ic = MP.index(imp)
            ax[row[ir]][imp] = fig.add_subplot( ist_grids[ic] )

    # Customize the colormap
    if calc_log:
        if Subset_data:
            color_intervals = [0,2.5]
            new_map = plt.cm.cividis.reversed()
        else:
            color_intervals = [0, 1, 2, 3, 4]  # Define your intervals
            exist_cmap = plt.cm.cividis.reversed()
            colors = exist_cmap(np.linspace(0,1,len(color_intervals)-1))
            new_map = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals)-1)
    else:
        if Subset_data:
            color_intervals = [0,150]
            # define the threshold values (boundaries) and corresponding colors
            #thresholds = [1,30,60,90,120,150]  # These are the threshold boundaries
            #colors = ['#dcbeff','#aaffc3','#f032e6','#4363d8','#911eb4','#f58231','#e6194B','#800000','#000000']
            # create a colormap using the list of colors
            #new_map = mcolors.ListedColormap(colors)
            #new_map = plt.cm.magma.reversed()
            # create a normalizer using BoundaryNorm
            #norm = mcolors.BoundaryNorm(thresholds, new_map.N)
            #color_intervals = [0,30,60,90,120,150]
            color_intervals = [0,15,30,45,60,75]
            #color_intervals = [0,500,1000,1500,2000,2500,3000,3500,4000]  # Define your intervals
            exist_cmap = plt.cm.inferno.reversed()
            colors = exist_cmap(np.linspace(0,1,len(color_intervals)-1))
            new_map = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals)-1)

        else:
            color_intervals = [0,10000]
            new_map = plt.cm.cividis.reversed()
            # define the threshold values (boundaries) and corresponding colors
            #thresholds = [10,50,250,500,1000,2250,5000,7500,10000,12500]  # These are the threshold boundaries
            #colors = ['#dcbeff','#aaffc3','#f032e6','#4363d8','#911eb4','#f58231','#e6194B','#800000','#000000']
            # create a colormap using the list of colors
            #new_map = mcolors.ListedColormap(colors)
            #new_map = plt.cm.Set1.reversed()
            # create a normalizer using BoundaryNorm
            #norm = mcolors.BoundaryNorm(thresholds, new_map.N)
            #color_intervals = [0,500,1000,1500,2000,2500,3000,3500,4000]  # Define your intervals
            #exist_cmap = plt.cm.cividis.reversed()
            #colors = exist_cmap(np.linspace(0,1,len(color_intervals)-1))
            #new_map = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals)-1)
            #new_map = plt.cm.cividis.reversed()


    # Plot 2d histogram
    for ir in row:
        if ir == 'bmo_Yo_conv':
            ida = 'CONV'
            iss = 'bmo_Yo'
        elif ir == 'bmo_Yo_ir':
            ida = 'IR'
            iss = 'bmo_Yo'
        elif ir == 'amo_Yo_ir':
            ida = 'IR'
            iss = 'amo_Yo'

        for imp in MP:
            show = ax[ir][imp].imshow( d_hcount[imp][ida][iss],cmap=new_map,aspect='equal', #norm=norm)
                            vmin=color_intervals[0],vmax=color_intervals[-1])
            #show = ax[ir][imp].imshow(np.transpose(d_hcount[imp][ida][iss]),cmap=new_map,aspect='equal',
            #                vmin=color_intervals[0],vmax=color_intervals[-1])

    # Add color bar below the plot
    caxes = fig.add_axes([0.9, 0.07, 0.02, 0.88])
    color_bar = fig.colorbar(show,cax=caxes,orientation='vertical',extend='max',ticks=color_intervals)
    #color_bar.ax.yaxis.set_major_locator(ticker.FixedLocator([0,1,2,3,4]))
    #color_bar.ax.yaxis.set_major_formatter(ticker.FixedFormatter(['$10^0$', '$10^1$', '$10^2$', '$10^3$','$10^4$']))
    #color_bar.ax.tick_params(labelsize=12)

    # ticks and labels
    for ir in row:
        for imp in MP:
            tick = list(range(0,number_bins+1,10*scale_factor))
            ytick_label = list(range(0,max_Tb_rg-min_Tb_rg+1,10))
            ytick_labels = [min_Tb_rg+it for it in ytick_label]
            xtick_label = list(range(0,max_Tbdiff_rg-min_Tbdiff_rg+1,10))
            xtick_labels = [min_Tbdiff_rg+it for it in xtick_label]
            # plot line where y=y0 is
            f_xinterp = interpolate.interp1d( xtick_labels, tick)
            ax[ir][imp].axvline( x = f_xinterp(0), color='grey', alpha=1, linestyle='dashed', linewidth=2)
            ax[ir][imp].axvline( x = f_xinterp(15), color='grey', alpha=0.6, linestyle='dashed', linewidth=2)
            ax[ir][imp].axvline( x = f_xinterp(-15), color='grey', alpha=0.6, linestyle='dashed', linewidth=2)
            # labels
            #ax[ir][imp].set_xticks( tick )
            ax[ir][imp].set_xticks( tick )
            ax[ir][imp].set_ylim(ymin=0,ymax=number_bins)
            ax[ir][imp].set_xlim(xmin=0,xmax=number_bins)
            # y ticklabel 
            if MP.index( imp ) == 0:
                ax[ir][imp].set_ylabel('Observation (K)',fontsize=12)
                ax[ir][imp].set_yticks( tick )
                ax[ir][imp].set_yticklabels( [str(it) for it in ytick_labels],fontsize=10 )
            else:
                ax[ir][imp].set_yticks( tick )
                ax[ir][imp].set_yticklabels( [''] * len(tick) )
            # x ticklabel
            if ir == 'bmo_Yo_conv':
                ax[ir][imp].text(0.95, 0.95, f"BMO", transform=ax[ir][imp].transAxes, 
                    fontsize=10, verticalalignment='top', horizontalalignment='right')
            elif ir == 'bmo_Yo_ir':
                ax[ir][imp].text(0.95, 0.95, f"BMO", transform=ax[ir][imp].transAxes,       
                    fontsize=10, verticalalignment='top', horizontalalignment='right')
            elif ir == 'amo_Yo_ir':
                ax[ir][imp].text(0.95, 0.95, f"AMO", transform=ax[ir][imp].transAxes,       
                    fontsize=10, verticalalignment='top', horizontalalignment='right')
            # x tick lable
            if ir == row[-1]:
                ax[ir][imp].set_xticks( tick )
                ax[ir][imp].set_xticklabels( [str(it) for it in xtick_labels],fontsize=10 )
                ax[row[-1]][imp].set_xlabel('Model minus Obs (K)',fontsize=10)
            else:
                ax[ir][imp].set_xticks( tick )
                ax[ir][imp].set_xticklabels( [''] * len(tick) )
    # text
    fig.text(0.03,0.82,'CONV',rotation='vertical',fontsize=12, ha='center', va='center')
    fig.text(0.03,0.53,'IR',rotation='vertical',fontsize=12, ha='center', va='center')
    fig.text(0.03,0.23,'IR',rotation='vertical',fontsize=12, ha='center', va='center')

    fig.text(0.32,0.96,MP[0], fontsize=12, ha='center', va='center')
    fig.text(0.73,0.96,MP[1], fontsize=12, ha='center', va='center')

    # Save
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/IR_2dPDF_BMO_'+str(len(d_hrs[Storms[0]]))+'cycles.png'
    plt.savefig( des_name )
    print( 'Saving the figure: ', des_name )
    plt.close()


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/Clean_results/'

    #--------Configuration------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    DA = ['CONV','IR']
    MP = ['WSM6','THO'] 
    sensor = 'abi_gr'

    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708231200','IRMA':'201709040000','JOSE':'201709060000','MARIA':'201709170000'}
    Consecutive_times = True

    # action
    Subset_data = False
    model_center = True
    max_radius = 200 #km

    # Range
    min_Tb_rg = 190
    max_Tb_rg = 260
    min_Tbdiff_rg = -35
    max_Tbdiff_rg = 35
    
    bin_BMO = True # background minue observation
    calc_log = False

    # number of bins
    scale_factor = 1 # 1: 1k per bin; 2: 0.5k per bin
    number_bins = (max_Tb_rg-min_Tb_rg)*scale_factor

    # Plot
    plot_obs_BMO = True

    #------------------------------------

    # Create experiment names
    Exper_names = {}
    for ist in Storms:
        Exper_names[ist] = {}
        for imp in MP:
            Exper_names[ist][imp] = {}
            for ida in DA:
                Exper_names[ist][imp][ida] = UD.generate_one_name( ist,ida,imp )

    # Identify DA times in the period of interest
    d_hrs = generate_times( Storms, start_time_str, end_time_str, 1 )

    # Read obs, Hxb, and Hxa of all files
    # Might need to subset data
    d_Tbs = {}
    for imp in MP:
        d_Tbs[imp] = {}
        for ida in DA:
            d_Tbs[imp][ida] = {}
            for ist in Storms:
                d_Tbs[imp][ida][ist] = Read_Subset_Tbs( ist,Exper_names[ist][imp][ida],d_hrs[ist] )


    # Combine data for all storms
    combined_Tb = {}
    for imp in MP:
        combined_Tb[imp] = {}
        for ida in DA:
            combined_Tb[imp][ida] = IRpdf.combine_storms_allTimes( Storms,d_Tbs[imp][ida] )

    # Make bins
    dcount = {}

    if bin_BMO:
        for imp in MP:
            dcount[imp] = {}
            for ida in DA:
                dcount[imp][ida] = {}
                Storms_Tb = combined_Tb[imp][ida]
                #!!!!!!! BMO V.S. Obs
                bmo_Yo = hist2d(Storms_Tb['Yo_obs'],Storms_Tb['meanYb_obs']-Storms_Tb['Yo_obs'],range=[[min_Tb_rg,max_Tb_rg],[min_Tbdiff_rg,max_Tbdiff_rg]],bins=number_bins)
                amo_Yo = hist2d(Storms_Tb['Yo_obs'],Storms_Tb['meanYa_obs']-Storms_Tb['Yo_obs'],range=[[min_Tb_rg,max_Tb_rg],[min_Tbdiff_rg,max_Tbdiff_rg]],bins=number_bins)
                print(imp,ida)
                print(np.sum( bmo_Yo ), np.sum(amo_Yo))
                # if Compute log10
                if calc_log:
                    # Compute log10 of non-zero elements, return NAN where elements are ZERO
                    with np.errstate(divide='ignore', invalid='ignore'):
                        dcount[imp][ida]['bmo_Yo'] = np.where(bmo_Yo != 0, np.log10(bmo_Yo), np.nan)
                        dcount[imp][ida]['amo_Yo'] = np.where(amo_Yo != 0, np.log10(amo_Yo), np.nan)
                else:
                    bmo_Yo[bmo_Yo == 0] = np.nan
                    amo_Yo[amo_Yo == 0] = np.nan
                    dcount[imp][ida]['bmo_Yo'] = bmo_Yo
                    dcount[imp][ida]['amo_Yo'] = amo_Yo

    else:
        for imp in MP:
            dcount[imp] = {}
            for ida in DA:
                dcount[imp][ida] = {}
                Storms_Tb = combined_Tb[imp][ida]
                #!!!!!!! Model V.S. Obs
                meanYb_Yo = hist2d(Storms_Tb['Yo_obs'],Storms_Tb['meanYb_obs'],range=[[min_Tb_rg,max_Tb_rg],[min_Tb_rg,max_Tb_rg]],bins=number_bins)
                meanYa_Yo = hist2d(Storms_Tb['Yo_obs'],Storms_Tb['meanYa_obs'],range=[[min_Tb_rg,max_Tb_rg],[min_Tb_rg,max_Tb_rg]],bins=number_bins)
                # Compute log10 of non-zero elements, return NAN where elements are ZERO
                with np.errstate(divide='ignore', invalid='ignore'):
                    dcount['meanYb_Yo'] = np.where(meanYb_Yo != 0, np.log10(meanYb_Yo), np.nan)
                    dcount['meanYa_Yo'] = np.where(meanYa_Yo != 0, np.log10(meanYa_Yo), np.nan)

    # Plot
    if plot_obs_BMO:
        plt_BMO( dcount )




























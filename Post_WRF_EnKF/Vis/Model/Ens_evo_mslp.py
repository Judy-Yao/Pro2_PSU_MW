import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import math
import matplotlib
from scipy import interpolate
matplotlib.use("agg")
import matplotlib.ticker as ticker
from matplotlib import pyplot as plt
from matplotlib import colors
import pandas as pd

import EnKF_Vmax as Vx #Vmax
import EnKF_minSlp as SC
import Util_data as UD

# Generate time series
def generate_times( start_time_str, end_time_str, interval ):

    time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
    time_diff_hour = time_diff.total_seconds() / 3600
    time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+interval, interval))]
    dict_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times

# Interpolate from 6-hour value to every hour value
def Interpolate_hourly( d_6hrs,d_obs6Hr,d_hrs):

    # group pairs of synoptic times 
    syno_ts = d_6hrs
    syno_pair = [[] for i in range(2)]
    for i in range(len(syno_ts)-1):
        syno_pair[0].append( syno_ts[i] )
        syno_pair[1].append( syno_ts[i+1] )

    # Linearly interpolate the Vmax between two synoptic times
    d_ObsHr = []

    for tt in d_hrs: #target time

        if any( hh in tt[8:10] for hh in ['00','06','12','18']):
            d_ObsHr.append( d_obs6Hr[tt]  )
        else:
            lg = [ syno_pair[0][ip] < tt < syno_pair[1][ip] for ip in range(len(syno_pair[0]))]
            idx = int( np.where(lg)[0] )

            # Calculate time differences
            time_diff_total = (datetime.strptime(syno_pair[1][idx],"%Y%m%d%H%M") - datetime.strptime(syno_pair[0][idx],"%Y%m%d%H%M")).total_seconds()
            time_diff_target = (datetime.strptime(tt,"%Y%m%d%H%M") - datetime.strptime(syno_pair[0][idx],"%Y%m%d%H%M")).total_seconds()

            # Calculate interpolation factor
            interp_factor = time_diff_target / time_diff_total

            # Perform linear interpolation for latitude and longitude
            d_ObsHr.append( d_obs6Hr[syno_pair[0][idx]] + interp_factor * (d_obs6Hr[syno_pair[1][idx]] - d_obs6Hr[syno_pair[0][idx]]) )

    return d_ObsHr


# Read Ens HPI info from pre-calculated files
# e.g., DAtime/DAtime_enkf_output_mslp.txt 
def read_mslp_ens( wrf_file  ):

    d_mslp = {}

    with open( wrf_file ) as f:
        next(f)
        all_lines = f.readlines()

    for line in all_lines:
        split_line = line.split()
        id_ = split_line[0]
        d_mslp[id_] = float(split_line[3])

    return d_mslp

# Read  HPI info from pre-calculated files
# e.g., HPI_wrf_enkf_input_d03_.201709030000_201709040000.txt 
def Read_mslp_EnsMean( mfile ):

    d_mean = {}
    print('Reading ', mfile)
    with open(mfile) as f:
        next(f)
        all_lines = f.readlines()

    for line in all_lines:
        if 'Success' in line:
            break
        split_line = line.split()
        time = split_line[0]
        lon = float(split_line[1])
        lat = float(split_line[2])
        mslp = = float(split_line[3])
        d_mean[time] = (lon,lat,mslp)

    return d_mean


def plot_slp_timeseries( small_dir, Storm, Expers, DAtimes, Evo_slp ):

    # Set up figure
    fig = plt.figure( figsize=(10,9), dpi=300 )
    ax = plt.subplot(1,1,1)
    dates = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]

    dates_zip = list( chain.from_iterable( zip(dates,dates)) )
    len_seg = len(dates_zip)
    # Customize labels
    labels = {}

    for iid in id_: # linestyle=(0, (5, 1))
        # zip data
        t_zip = list( chain.from_iterable( zip(lead_t,lead_t)) )
        len_seg = len( t_zip )
        slp_zip = list( chain.from_iterable( zip(Exper_minSLP[istorm][imp][ida]['xb_slp'],Exper_minSLP[istorm][imp][ida]['xa_slp']) ) )

        if ida == 'CONV':
            lines = '-'
        elif ida == 'IR':
            lines = '--'
        else:
            lines = (0, (1, 1))

        for i in range(1,len_seg):
        # specify which segment uses which line style
            if i % 2 == 0:
                line = lines #'-' #'--'
            else:
                line = lines # '-'
        # customize the line
            if i == 1:
                ax.plot(t_zip[i-1:i+1],slp_zip[i-1:i+1],linestyle=line,linewidth='6',color=colorset[istorm])
            elif i == 2:
                ax.plot(t_zip[i-1:i+1],slp_zip[i-1:i+1],linestyle=line,linewidth='4',color=colorset[istorm],label=ida )#istorm)
            else:
                ax.plot(t_zip[i-1:i+1],slp_zip[i-1:i+1],linestyle=line,linewidth='4',color=colorset[istorm])

    # Plot 
    if if_:
        pass



    fig.suptitle( supt, fontsize=17, fontweight='bold')

    leg = plt.legend(loc='lower left',fontsize=22)
    # Set X/Y labels
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45,labelsize=15)
    ax.set_ylabel('Minimum Sea Level Pressure (hPa)',fontsize=24)
    ax.set_ylim(920,1020)  #( 900,1000 ) #(970,1020)

    #ax.set_title( 'mim SLP (hPa)',fontweight="bold",fontsize='15' )

    # Save the figure
    save_des = small_dir+Storm+'/'+Expers[1]+'/Vis_analyze/Model/minslp_'+DAtimes[0]+'_'+DAtimes[-1]+'_WSM6_THO.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()



def plot_bins_with_t( count, ix, mean=None, btk=None):

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)
    
    # Set up coordinate
    x_rg = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]
    x_axis_edges = range(len(x_rg)+1)
    # Replace -inf and inf with np.nan
    y_rg = [np.nan if x == float('inf') or x == -float('inf') else x for x in bins]
    # Interpolate the bin centers (y axis)
    y_axis_rg = range(len(y_rg))
    f_yinterp = interpolate.interp1d( y_rg, y_axis_rg)

    # Plot occurence frequency
    # percentage
    vmin = 0
    vmax = 80
    cs = ax.pcolormesh( x_axis_edges,y_axis_rg[::-1], np.flipud(count), cmap='gist_heat_r', vmin=vmin,vmax=vmax)
    #ax.scatter(0,f_yinterp( 940 )  ): test
    #cs = ax.imshow(count,aspect='auto', cmap='gist_heat_r',origin='upper',vmin=vmin,vmax=vmax)

    # colorbar
    cbar = fig.colorbar( cs,extend='max',fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel('Percentage (%)',fontsize=12)
    # x_axis_centers
    x_axis_centers = [(x_axis_edges[i] + x_axis_edges[i + 1]) / 2 for i in range(len(x_axis_edges) - 1)]
    formatted_date = [it.strftime("%m/%d %H UTC") for it in x_rg]
    ax.xaxis.set_major_locator(ticker.FixedLocator( x_axis_centers[::6] ))
    ax.set_xticklabels( formatted_date[::6] )
    ax.xaxis.set_minor_locator(ticker.FixedLocator( x_axis_centers ))
    ax.tick_params(axis='x', labelrotation=15, labelsize=12)
    # yticks
    yticks = y_axis_rg  #[i for i in y_axis_rg if not np.isnan(y_rg[i])]
    ytick_labels = [str(i) for i in y_rg if not np.isnan(i) ]
    ytick_labels = ['-Inf'] + ytick_labels + ['Inf']
    ax.set_yticks( yticks )
    ax.set_yticklabels( ytick_labels,fontsize=15 )
    if not at_obs_res:
        ax.set_ylabel('MSLP (hPa)',fontsize=15)
    else:
        ax.set_ylabel('SLP (hPa)',fontsize=15)

    # Add line for ensemble mean
    if if_btk:
        btk_axis_y = []
        for it in btk:
            if it <= 930:
                btk_axis_y.append( 0 )
            elif it >= 1000:
                btk_axis_y.append( 9 )
            else:
                btk_axis_y.append( f_yinterp( it ) )
        ax.plot(x_axis_centers, btk_axis_y, color='black', linewidth=2, marker='o', label='Best Track')
        ax.legend(loc='upper right',fontsize=12)

    # Add line for ensemble mean
    if if_mean:
        mslp_mean = list( mean.values() )
        mean_axis_y = []
        for it in mslp_mean:
            if it <= 930:
                mean_axis_y.append( 0 )
            elif it >= 1000:
                mean_axis_y.append( 9 )
            else:
                mean_axis_y.append( f_yinterp( it ) )
        ax.plot(x_axis_centers, mean_axis_y, color='blue', linewidth=2, marker='o', label='Ens Mean')
        ax.legend(loc='upper right',fontsize=12)

    # titles
    if not at_obs_res:
        if ix == 'xb':
            title_name = 'Occurrence Frequency for MSLP in Ens EnKF input (Xb)'
        else:
            title_name = 'Occurrence Frequency for MSLP in Ens EnKF output (Xa)'
    else:
        if ix == 'xb':
            title_name = 'Occurrence Frequency for SLP at Assimilated MSLP Loc in Ens EnKF input (Xb)'
        else:
            title_name = 'Occurrence Frequency for SLP at Assimilated MSLP Loc in Ens EnKF output (Xa)'
    ax.set_title( title_name,fontweight="bold",fontsize=11 )
    fig.suptitle(Storm+': '+Exper_name, fontsize=10, fontweight='bold')

    # Save the figure
    if not at_obs_res:
        if ix == 'xb':
            save_des = plot_dir+'Freq_mslp_ens_input_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
        else:
            save_des = plot_dir+'Freq_mslp_ens_output_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    else:
        if ix == 'xb':
            save_des = plot_dir+'Freq_slp_obsLoc_ens_input_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
        else:
            save_des = plot_dir+'Freq_slp_obsLoc_ens_output_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'

    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR-onlyTCV'
    MP = 'WSM6'
    
    Num_ens = 60
    # Time range set up
    start_time_str = '201709030000'
    end_time_str = '201709040000'
    Consecutive_times = True

    # at obs location
    at_obs_res = False

    # Actions
    if_mean = True # read ensemble mean or not
    if_tcvital = False
    if_btk = True

    Bin_w_time = True
    if Bin_w_time:
        # Define the bin edges and labels
        bins = [-float('inf'),930,940,950,960,970,980,990,1000,float('inf')] # hPa
        bin_labels = ['<930','930-940','940-950','950-960','960-970','970-980','980-990','990-1000','>1000']
    
    # find Ens mslp: identify_mslp_ens in ../EnKF/explain_horizontal_onepoint_corr.py 
    # find  mslp: Identify_EnKF_StormCenter.py in ../Paper1/Identify_EnKF_StormCenter.py 
    # ------------------------------------------------------   

    DAtimes = generate_times( start_time_str, end_time_str, 1 )
    t_6hrs = generate_times(  start_time_str, end_time_str, 6 )

    # Experiment name
    Exper_name = UD.generate_one_name( Storm,DA,MP )

    # Read best-track
    if if_btk:
        btk6Hrs = SC.MSLP_btk( small_dir,Storm,t_6hrs )
        btkHr = Interpolate_hourly( t_6hrs,btk6Hrs['mslp'],DAtimes)

    # Read mslp of the whole ensemble
    xb_ens = {}
    xa_ens = {}
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
        if not at_obs_res:
            xb_ens[DAtime] = read_mslp_ens( wrf_dir+ DAtime+'_enkf_input_mslp.txt' )
            xa_ens[DAtime] = read_mslp_ens( wrf_dir+ DAtime+'_enkf_output_mslp.txt'  )
        else:
            xb_ens[DAtime] = read_mslp_ens( wrf_dir+ DAtime+'_enkf_input_slp_obsLoc.txt' )
            xa_ens[DAtime] = read_mslp_ens( wrf_dir+ DAtime+'_enkf_output_slp_obsLoc.txt'  )

    # Read mslp of the ensemble mean 
    if if_mean:
        saved_dir = small_dir+Storm+'/'+Exper_name+'/Data_analyze/'
        xb_mean = Read_mslp_mean( saved_dir+'HPI_wrf_enkf_input_d03_mean.'+start_time_str+'_'+end_time_str+'.txt' )
        xa_mean = Read_mslp_mean( saved_dir+'HPI_wrf_enkf_output_d03_mean.'+start_time_str+'_'+end_time_str+'.txt' )
    
    # Bin the mslp for the ensemble with time
    if Bin_w_time:
        d_hcount = {}
        # Data process
        xb_hcount = np.zeros( (len(bin_labels),len(DAtimes)) )
        xa_hcount = np.zeros( (len(bin_labels),len(DAtimes)) )
        d_hcount['xb'] = xb_hcount
        d_hcount['xa'] = xa_hcount
        for DAtime in DAtimes:
            it = DAtimes.index( DAtime )
            # FOR Xb
            # create a DataFrame
            df_xb = pd.DataFrame(list(xb_ens[DAtime].values()), columns=['values'])
            # bin the data and count the occurrences
            df_xb['binned'] = pd.cut(df_xb['values'], bins=bins, labels=bin_labels)
            tmp_count = df_xb['binned'].value_counts().sort_index()
            d_hcount['xb'][:,it] = tmp_count.values
            # FOR Xa
            # create a DataFrame
            df_xa = pd.DataFrame(list(xa_ens[DAtime].values()), columns=['values'])
            # bin the data and count the occurrences
            df_xa['binned'] = pd.cut(df_xa['values'], bins=bins, labels=bin_labels)
            tmp_count = df_xa['binned'].value_counts().sort_index()
            d_hcount['xa'][:,it] = tmp_count.values

        # Plot
        for ix in ['xb','xa']:
            plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/Ens_MSLP_time/'
            plotdir_exists = os.path.exists( plot_dir )
            if plotdir_exists == False:
                os.mkdir(plot_dir)
           
            if ix == 'xb':
                plot_bins_with_t( d_hcount[ ix ]/Num_ens*100, ix, mean=xb_mean, btk=btkHr)
            else:
                plot_bins_with_t( d_hcount[ ix ]/Num_ens*100, ix, mean=xa_mean, btk=btkHr)
        

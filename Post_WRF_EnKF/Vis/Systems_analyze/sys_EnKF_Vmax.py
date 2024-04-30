import os
import glob
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib
import math
from datetime import datetime, timedelta
import time
import netCDF4 as nc
from wrf import getvar
import subprocess
from itertools import chain

import Util_data as UD

# Calculate mean absolute error
def mean_absolute_error(obs,model):
    return np.mean(np.abs(np.subtract(obs, model)))

# Generate time series
def generate_times( Storms, start_time_str, end_time_str, interval ):

    dict_times = {}
    for istorm in Storms:
        time_diff = datetime.strptime(end_time_str[istorm],"%Y%m%d%H%M") - datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+interval, interval))]
        dict_times[istorm] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times

def find_Vmax( wrfout ):

    ncdir = nc.Dataset( wrfout )
    U10 = ncdir.variables['U10'][0,:,:]
    V10 = ncdir.variables['V10'][0,:,:]
    windspeed = (U10 ** 2 + V10 ** 2) ** 0.5
        
    return np.amax(windspeed)

# Obtain max 10-meter wind from WRF output
def model_Vmax( big_dir,istorm,iExper,DAtimes,Vmax_xa,Vmax_xb=None ):

    dict_Vmax = {}

    if Vmax_xa:
        xa_vmax_value = []

    if Vmax_xb:
        xb_vmax_value = []

    for DAtime in DAtimes:
        if Vmax_xa:
            xa_dir = big_dir+istorm+'/'+iExper+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
            print('Reading the EnKF posterior mean from ', xa_dir)
            xa_vmax_value.append( find_Vmax( xa_dir ) )

        if Vmax_xb:
            xb_dir = big_dir+istorm+'/'+iExper+'/fc/'+DAtime+'/wrf_enkf_input_d03_mean'
            print('Reading the EnKF prior mean from ', xb_dir)
            xb_vmax_value.append( find_Vmax( xb_dir ) )

    if Vmax_xa:
        dict_Vmax['xa_vmax'] = xa_vmax_value
    if Vmax_xb:
        dict_Vmax['xb_vmax'] = xb_vmax_value

    return dict_Vmax

# Obtain sustained max wind speed from best-track file (m/s)
def Vmax_btk( small_dir,Storm,DAt_6hrs ):

    btk_file = os.listdir(small_dir+Storm+'/Post_Storm_btk/')
    with open (small_dir+Storm+'/Post_Storm_btk/'+btk_file[0]) as f:
        btk_all = f.readlines()

    # btk times
    btk_times = []
    for line in btk_all:
        split_line = line.split()
        btk_times.append(split_line[2].replace(',','') + '00')

    # Filter out repetitive times
    idx_lines = []
    for itx in DAt_6hrs:
        lg = [itx == itobs for itobs in btk_times]
        idx_lines.append( np.where(lg)[0][0] ) # only choose the first occurrence 

    # Vmax
    d_obs6Hr = {}
    for i in idx_lines:
        line = btk_all[i]
        split_line = line.split()
        btk_time = split_line[2].replace(',','') + '00'
        d_obs6Hr[btk_time] = float(split_line[8].replace(',',''))*0.51444 # knots to m/s
    
    return d_obs6Hr


# Obtain sustained max wind speed from TCvitals (m/s)
def Vmax_TCvitals( small_dir,Storm,DAt_6hrs ):

    tc_file = small_dir+Storm+'/TCvitals/'+Storm+'_tcvitals'
    with open(tc_file) as tmp:
        tc_all = tmp.readlines()

    # TCvital times (may have repetitive values)
    tc_times = []
    for line in tc_all:
        line_split = line.split()
        tc_times.append( line_split[3]+line_split[4] )

    # Filter out repetitive times
    idx_lines = []
    for itx in DAt_6hrs:
        lg = [itx == itobs for itobs in tc_times]
        idx_lines.append( np.where(lg)[0][0] ) # only choose the first occurrence
    
    d_obs6Hr = {}
    for i in idx_lines:
        line = tc_all[i]
        line_split = line.split()
        tc_time = line_split[3]+line_split[4] 
        d_obs6Hr[tc_time] = float(line_split[11])/10 

    return d_obs6Hr


# Interpolate from 6-hour value to every hour value
def Interpolate_hourly( Storm,d_6hrs,d_obs6Hr,d_hrs):

    # group pairs of synoptic times 
    syno_ts = d_6hrs[Storm]
    syno_pair = [[] for i in range(2)]
    for i in range(len(syno_ts)-1):
        syno_pair[0].append( syno_ts[i] )
        syno_pair[1].append( syno_ts[i+1] )

    # Linearly interpolate the Vmax between two synoptic times
    d_ObsHr = []

    for tt in d_hrs[Storm]: #target time

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


def plot_sys_Vmax():

    # Obtain hourly TC-vital max wind speed
    d_obs = {}
    for istorm in Storms:
        d_obs6Hrs = Vmax_TCvitals( istorm,d_6hrs[istorm] )
        d_obs[istorm] = Interpolate_hourly( istorm,d_6hrs,d_obs6Hrs,d_hrs)

    # Obtain model simulated min vmax
    Exper_Vmax = {}
    for istorm in Storms:
        Exper_Vmax[istorm] = {}
        for imp in MP:
            Exper_Vmax[istorm][imp] = {}
            for ida in DA:
                iExper = Exper_names[istorm][imp][ida]
                if iExper is None:
                    Exper_Vmax[istorm][imp][ida] = None
                else:
                    Exper_Vmax[istorm][imp][ida] = model_Vmax( istorm,Exper_names[istorm][imp][ida],d_hrs[istorm] )

    # Obtain mean absolute error
    MAE = {}
    for istorm in Storms:
        MAE[istorm] = {}
        for imp in MP:
            MAE[istorm][imp] = {}
            for ida in DA:
                MAE[istorm][imp][ida] = {}
                if Vmax_xb:
                    MAE[istorm][imp][ida]['xb'] = mean_absolute_error(d_obs[istorm],Exper_Vmax[istorm][imp][ida]['xb_vmax'])
                if Vmax_xa:
                    if iExper is None:
                        MAE[istorm][imp][ida]['xa'] = None
                    else:
                        MAE[istorm][imp][ida]['xa'] = mean_absolute_error(d_obs[istorm],Exper_Vmax[istorm][imp][ida]['xa_vmax'])

    # Set up figure
    fig = plt.figure( figsize=(18,9), dpi=300 )
    ax = plt.subplot(1,1,1)

    # Set colors
    colorset = {'HARVEY':'#FF13EC','JOSE':'#0D13F6','MARIA':'#FFA500','IRMA':'#DB2824'}
    lead_t = list(range(0, cycles, 1))

    # Plot obs
    for istorm in Storms:
        ax.plot(lead_t, d_obs[istorm],color='black',linestyle='-',linewidth=3)#,label=istorm)

    # plot Vmax_xa
    if Vmax_xa and not Vmax_xb:
        for istorm in Storms:
            for imp in MP:
                for ida in DA: # linestyle=(0, (5, 1))
                    if ida == 'CONV':
                        ax.plot(lead_t, Exper_Vmax[istorm][imp][ida]['xa_vmax'],color=colorset[istorm],linestyle='-',linewidth=4,label=istorm) # linestyle=(0, (5, 1))
                    elif ida == 'IR':
                        ax.plot(lead_t, Exper_Vmax[istorm][imp][ida]['xa_vmax'],color=colorset[istorm],linestyle='--',linewidth=4)
                    else:
                        ax.plot(lead_t, Exper_Vmax[istorm][imp][ida]['xa_vmax'],color=colorset[istorm],linestyle=(0, (1, 1)),linewidth=4)
                        #ax.scatter(lead_t, Exper_Vmax[istorm][imp][ida]['xa_vmax'],s=5,color=colorset[istorm],marker='*',markersize=12)


    # plot saw-tooth lines
    if Vmax_xa and Vmax_xb:
        for istorm in Storms:
            for imp in MP:
                for ida in DA: # linestyle=(0, (5, 1))
                    # zip data
                    t_zip = list( chain.from_iterable( zip(lead_t,lead_t)) )
                    len_seg = len( t_zip )
                    Vmax_zip = list( chain.from_iterable( zip(Exper_Vmax[istorm][imp][ida]['xb_vmax'],Exper_Vmax[istorm][imp][ida]['xa_vmax']) ) )

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
                            ax.plot(t_zip[i-1:i+1],Vmax_zip[i-1:i+1],linestyle=line,linewidth='6',color=colorset[istorm])
                        elif i == 2:
                            ax.plot(t_zip[i-1:i+1],Vmax_zip[i-1:i+1],linestyle=line,linewidth='4',color=colorset[istorm],label=ida )#istorm)
                        else:
                            ax.plot(t_zip[i-1:i+1],Vmax_zip[i-1:i+1],linestyle=line,linewidth='4',color=colorset[istorm])

    # labels and attributes
    leg = ax.legend(loc='lower left',fontsize=25)
    # Set X/Y labels
    ax.set_xticks( lead_t[::4] )
    ax.set_xlim([0,cycles-1])
    ax.set_xlabel('EnKF Cycle',fontsize=28)
    ax.grid(True,linewidth=1, color='gray', alpha=0.5, linestyle='-')
    ax.set_ylabel('Maximum 10-m Wind (m/s)',fontsize=28)
    ax.set_ylim( 0,70 )
    ax.tick_params(axis='both',labelsize=28)

    # Titles
    supt = 'mean absolute error of Xa\n'
    for ida in DA:
        for imp in MP:
            supt = supt+ida+'_'+imp+': '
            storm_str = ''
            for istorm in Storms:
                pass
                storm_str = storm_str+istorm+' '+'%.2f' %MAE[istorm][imp][ida]['xa']+'; '
            supt = supt+storm_str+'\n'
    fig.suptitle( supt, fontsize=15, fontweight='bold')

    #if Vmax_xa and not Vmax_xa:
    #    fig.suptitle('TCvtial V.S. Analysis', fontsize=24, fontweight='bold')

    # Save the figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Model/HARVEY_Vmax_THO_CONV_IR_MW.png'
    plt.savefig( des_name )
    print( 'Saving the figure: ', des_name )
    plt.close()


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY']
    DA = ['IR','IR+MW']
    MP = ['THO',] #

    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708231200','IRMA':'201709040000','JOSE':'201709060000','MARIA':'201709170000'}

    Vmax_xa = True
    Vmax_xb = False
    cycles = 25
    Plot_vmax_compare = True

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = {}
        for imp in MP:
            Exper_names[istorm][imp] = {}
            for ida in DA:
                Exper_names[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # Identify DA times in the period of interest
    d_hrs = generate_times( Storms, start_time_str, end_time_str, 1 )
    d_6hrs = generate_times( Storms, start_time_str, end_time_str, 6 )

    # Compare analysis mean with TCvital
    if Plot_vmax_compare:
        start_time=time.process_time()
        plot_sys_Vmax()
        end_time = time.process_time()
        print('time needed: ', end_time-start_time, ' seconds')






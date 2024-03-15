#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7
import os
import glob
import numba
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
import scipy as sp
import scipy.ndimage
from itertools import chain

import Util_data as UD

matplotlib.rcParams['font.size'] = 20

# Calculate mean absolute error
def mean_absolute_error(obs,model):
    return np.mean(np.abs(np.subtract(obs, model)))

# Generate time series
def generate_times( Storms, start_time_str, end_time_str ):

    dict_times = {}
    for istorm in Storms:
        time_diff = datetime.strptime(end_time_str[istorm],"%Y%m%d%H%M") - datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        dict_times[istorm] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times

# Obtain assimilated min slp from fort.10000
def assimilated_obs( storm,iExper,DAtimes ):

    # Read model min slp
    d_model = model_minSLP( storm,iExper,DAtimes )

    obs_minslp_lat = []
    obs_minslp_lon = []
    obs_minslp_value = []

    for DAtime in DAtimes:
        # collect assimilated min slp obs from TCvital record in fort.10000
        diag_enkf = big_dir+storm+'/'+iExper+'/run/'+DAtime+'/enkf/d03/fort.10000'
        print('Reading the EnKF diagnostics from ', diag_enkf)
        enkf_minSlp = subprocess.run(["grep","slp",diag_enkf],stdout=subprocess.PIPE,text=True)
        list_enkf_minSlp = enkf_minSlp.stdout.split()
        # condition on the number of assimilated min slp
        if list_enkf_minSlp.count('slp') == 0 :
            raise ValueError('No min slp is assimilated!')
        elif list_enkf_minSlp.count('slp') == 1 :
            obs_minslp_lat.append( float(list_enkf_minSlp[1]) )
            obs_minslp_lon.append( float(list_enkf_minSlp[2]) )
            obs_minslp_value.append( float(list_enkf_minSlp[9])/100 )
        else : # at least two min slp records are assimilated
            print('At least two min slp obs are assimilated!')
            # find the index of the current time
            idx_time = DAtimes.index( DAtime )
            # assemble the diagnosed min slp from an analysis
            xa_ms_lat = d_model['xa_lat'][idx_time]
            xa_ms_lon = d_model['xa_lon'][idx_time]
            # ---condition 1: find the nearest TCvital min slp from the analysis
            # find the index/location of 'slp' in fort.10000
            indices = [i for i ,e in enumerate(list_enkf_minSlp) if e == 'slp']
            # assemble a pair of coordinate for each 'slp'
            distances = []
            obs_slp = []
            for it in indices:
                obs_slp.append( float(list_enkf_minSlp[it+9]) )
                lon1 = float(list_enkf_minSlp[it+2])
                lat1 = float(list_enkf_minSlp[it+1])
                distances.append( UD.mercator_distance(lon1, lat1, xa_ms_lon, xa_ms_lat) )
            min_distances = [np.amin(distances) == it for it in distances]
            idx_nearest = min_distances.index(True)
            # ---condition 2: the min slp
            min_obs_slp = [np.amin(obs_slp) == it for it in obs_slp]
            idx_min = min_obs_slp.index(True)
            # ---combine the two conditions
            if idx_min == idx_nearest:
                idx_coor = idx_min
                print('Storm center is choosed with two condtions met!')
            else:
                print('Not sure which obs is the storm center. Go with the min value one!')
                idx_coor = idx_min
            # gather this TCvital min slp
            obs_minslp_lat.append( float(list_enkf_minSlp[indices[idx_coor]+1]) )
            obs_minslp_lon.append( float(list_enkf_minSlp[indices[idx_coor]+2]) )
            obs_minslp_value.append( float(list_enkf_minSlp[indices[idx_coor]+9])/100 )

    # Assemble the dictionary
    return obs_minslp_value


def find_minSLP( wrfout ):
    
    ncdir = nc.Dataset( wrfout )
    slp = getvar(ncdir, 'slp')
    min_slp = np.amin( slp )
    slp_smooth = sp.ndimage.gaussian_filter(slp, [11,11])
    idx = np.nanargmin( slp_smooth )
    lat_minslp = ncdir.variables['XLAT'][:].flatten()[idx]
    lon_minslp = ncdir.variables['XLONG'][:].flatten()[idx]

    min_slp = [lat_minslp,lon_minslp,min_slp.values]
    return min_slp


# Obtain min slp from WRF output
def model_minSLP( istorm,iExper,DAtimes ):

    dict_minSLP = {}

    if slp_xa:
        xa_minslp_lat = [] 
        xa_minslp_lon = [] 
        xa_minslp_value = [] 

    if slp_xb:
        xb_minslp_lat = [] 
        xb_minslp_lon = [] 
        xb_minslp_value = [] 

    for DAtime in DAtimes:
        if slp_xa:
            xa_dir = big_dir+istorm+'/'+iExper+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
            print('Reading the EnKF posterior mean from ', xa_dir)
            list_xa = find_minSLP( xa_dir )
            xa_minslp_lat.append( list_xa[0] )
            xa_minslp_lon.append( list_xa[1] )
            xa_minslp_value.append( list_xa[2] )

        if slp_xb:
            xb_dir = big_dir+istorm+'/'+iExper+'/fc/'+DAtime+'/wrf_enkf_input_d03_mean'
            print('Reading the EnKF prior mean from ', xb_dir)
            list_xb = find_minSLP( xb_dir )
            xb_minslp_lat.append( list_xb[0] )
            xb_minslp_lon.append( list_xb[1] )
            xb_minslp_value.append( list_xb[2] )

    if slp_xa:
        dict_minSLP['xa_slp'] = xa_minslp_value
        dict_minSLP['xa_lat'] = xa_minslp_lat
        dict_minSLP['xa_lon'] = xa_minslp_lon
    if slp_xb:
        dict_minSLP['xb_slp'] = xb_minslp_value
        dict_minSLP['xb_lat'] = xb_minslp_lat
        dict_minSLP['xb_lon'] = xb_minslp_lon

    return dict_minSLP


def plot_sys_minslp():

    # Obtain assimilated TC-vital min slp
    d_tcvital = {}
    for istorm in Storms:
        d_tcvital[istorm] = assimilated_obs( istorm, Exper_names[istorm]['WSM6']['IR'], dict_times[istorm]) 

    # Obtain model simulated min slp
    Exper_minSLP = {}
    for istorm in Storms:
        Exper_minSLP[istorm] = {}
        for imp in MP:
            Exper_minSLP[istorm][imp] = {}
            for ida in DA:
                iExper = Exper_names[istorm][imp][ida]
                if iExper is None:
                    Exper_minSLP[istorm][imp][ida] = None
                else:
                    Exper_minSLP[istorm][imp][ida] = model_minSLP( istorm,Exper_names[istorm][imp][ida],dict_times[istorm] )

    # Obtain mean absolute error
    MAE = {}
    for istorm in Storms:
        MAE[istorm] = {}
        for imp in MP:
            MAE[istorm][imp] = {}
            for ida in DA:
                MAE[istorm][imp][ida] = {}
                if slp_xb:
                    MAE[istorm][imp][ida]['xb'] = mean_absolute_error(d_tcvital[istorm],Exper_minSLP[istorm][imp][ida]['xb_slp']) 
                if slp_xa:
                    MAE[istorm][imp][ida]['xa'] = mean_absolute_error(d_tcvital[istorm],Exper_minSLP[istorm][imp][ida]['xa_slp'])


    # Set up figure
    fig = plt.figure( figsize=(18,9), dpi=300 )
    ax = plt.subplot(1,1,1) 

    # Set colors
    colorset = {'HARVEY':'#FF13EC','JOSE':'#0D13F6','MARIA':'#FFA500','IRMA':'#DB2824'}
    lead_t = list(range(0, cycles, 1))
   
    # Plot obs
    for istorm in Storms:
        ax.plot(lead_t, d_tcvital[istorm],color='black',linestyle='-',linewidth=3)#,label=istorm)

    # plot slp_xa
    if slp_xa and not slp_xb:
        for istorm in Storms:
            for imp in MP:
                for ida in DA: # linestyle=(0, (5, 1))
                    if imp == 'WSM6':
                        ax.plot(lead_t, Exper_minSLP[istorm][imp][ida]['xa_slp'],color=colorset[istorm],linestyle='-',linewidth=4,label=istorm) # linestyle=(0, (5, 1))
                    elif imp == 'THO':
                        ax.plot(lead_t, Exper_minSLP[istorm][imp][ida]['xa_slp'],color=colorset[istorm],linestyle='--',linewidth=4,markersize=12)
                    else:
                        ax.plot(lead_t, Exper_minSLP[istorm][imp][ida]['xa_slp'],color=colorset[istorm],marker='*',linewidth=4,markersize=12)
                        #ax.scatter(lead_t, Exper_minSLP[istorm][imp][ida]['xa_slp'],s=5,color=colorset[istorm],marker='*')
                        
    # plot saw-tooth lines
    if slp_xa and  slp_xb:
        for istorm in Storms:
            for imp in MP:
                for ida in DA: # linestyle=(0, (5, 1))
                    # zip data
                    t_zip = list( chain.from_iterable( zip(lead_t,lead_t)) )
                    len_seg = len( t_zip )
                    slp_zip = list( chain.from_iterable( zip(Exper_minSLP[istorm][imp][ida]['xb_slp'],Exper_minSLP[istorm][imp][ida]['xa_slp']) ) )

                    for i in range(1,len_seg):
                    # specify which segment uses which line style
                        if i % 2 == 0:
                            line = '-' #'--'
                        else:
                            line = '-'
                    # customize the line
                        if i == 1:
                            ax.plot(t_zip[i-1:i+1],slp_zip[i-1:i+1],linestyle=line,linewidth='6',color=colorset[istorm])
                        elif i == 2:
                            ax.plot(t_zip[i-1:i+1],slp_zip[i-1:i+1],linestyle=line,linewidth='4',color=colorset[istorm],label=istorm)
                        else:
                            ax.plot(t_zip[i-1:i+1],slp_zip[i-1:i+1],linestyle=line,linewidth='4',color=colorset[istorm])


    # labels and attributes
    leg = ax.legend(loc='lower left',fontsize=25)
    # Set X/Y labels
    ax.set_xticks( lead_t[::4] )
    ax.set_xlim([0,cycles-1])
    ax.set_xlabel('EnKF Cycle',fontsize=28)
    ax.grid(True,linewidth=1, color='gray', alpha=0.5, linestyle='-')
    ax.set_ylabel('Minimum Sea Level Pressure (hPa)',fontsize=28)
    ax.set_ylim(970,1020) #( 920,1020 )
    ax.tick_params(axis='both',labelsize=28)

    # Titles
    supt = 'mean absolute error of Xa\n'
    for ida in DA:
        for imp in MP:
            supt = supt+ida+'_'+imp+': '
            storm_str = ''
            for istorm in Storms:
                storm_str = storm_str+istorm+' '+'%.2f' %MAE[istorm][imp][ida]['xa']+'; '
            supt = supt+storm_str+'\n'
    fig.suptitle( supt, fontsize=17, fontweight='bold')
 
    #if slp_xa and not slp_xb:
    #    fig.suptitle('TCvtial V.S. Analysis', fontsize=24, fontweight='bold')

    # Save the figure
    des_name = 'Harvey.png'#small_dir+'SYSTEMS/Vis_analyze/Model/sys_minslp_IR_WSM6_IR_TuneWSM6.png'
    plt.savefig( des_name )
    print( 'Saving the figure: ', des_name )
    plt.close()


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY',]#'JOSE','MARIA','IRMA']
    DA = ['IR',]
    MP = ['WSM6','THO','TuneWSM6'] #

    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708241200','IRMA':'201709050000','JOSE':'201709070000','MARIA':'201709180000'}

    slp_xa = True
    slp_xb = False
    Plot_minslp_compare = True
    cycles = 49

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = {}
        for imp in MP:
            Exper_names[istorm][imp] = {}
            for ida in DA:
                Exper_names[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # Identify DA times in the period of interest
    dict_times = generate_times( Storms, start_time_str, end_time_str )

    # Compare analysis mean with TCvital
    if Plot_minslp_compare:
        start_time=time.process_time()
        plot_sys_minslp()
        end_time = time.process_time()
        print('time needed: ', end_time-start_time, ' seconds')




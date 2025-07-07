#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7
import os,fnmatch
import glob
import numba
import numpy as np
#import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib
import math
from datetime import datetime, timedelta
import time
import netCDF4 as nc
import subprocess
import scipy as sp
import scipy.ndimage
from itertools import chain
import numpy.ma as ma

import Util_data as UD
import EnKF_Vmax as Vx 

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

def Read_HPI_file( file_hpi ):

    with open(file_hpi, 'r') as file:
         # Skip the first 21 lines and Read all lines into a list
         lines = file.readlines()[21:]
    hpi_line = lines[0].split()
    hpi_lat = float(hpi_line[4])
    hpi_lon = float(hpi_line[5])
    return (hpi_lon,hpi_lat)


# Obtain assimilated min slp from fort.10000
def assimilated_obs(big_dir,small_dir,Storm,iExper,DAtimes,slp_xa,slp_xb=False  ):

    # Read model min slp
    #d_model = model_minSLP( big_dir,Storm,iExper,DAtimes,slp_xa )

    # Determine the min slp and its location
    d_obs = {}

    obs_minslp_lat = []
    obs_minslp_lon = []
    obs_minslp_value = []

    for DAtime in DAtimes:
        # collect assimilated min slp obs from TCvital record in fort.10000
        diag_enkf = big_dir+Storm+'/'+iExper+'/run/'+DAtime+'/enkf/d03/fort.10000'
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
            dt = datetime.strptime(DAtime, "%Y%m%d%H%M")
            DAtime_wrf = dt.strftime("%Y-%m-%d_%H:%M")
            file_hpi = small_dir+'/Obs_input_EnKF/'+Storm+'/HPI/HPI_obs_gts_'+DAtime_wrf+':00.3DVAR'
            print('Reading '+file_hpi+' for more information......')
            hpi_loc = Read_HPI_file( file_hpi )
            # find the index/location of 'slp' in fort.10000
            indices = [i for i ,e in enumerate(list_enkf_minSlp) if e == 'slp']
            # assemble a pair of coordinate for each 'slp'
            for it in indices:
                lon = float(list_enkf_minSlp[it+2])
                lat = float(list_enkf_minSlp[it+1])
                loc_candis = (lon,lat)
                if loc_candis == hpi_loc:
                    obs_minslp_lat.append( float(list_enkf_minSlp[it+1]) )
                    obs_minslp_lon.append( float(list_enkf_minSlp[it+2]) )
                    obs_minslp_value.append( float(list_enkf_minSlp[it+9])/100 )

            # find the index of the current time
#            idx_time = DAtimes.index( DAtime )
#            # assemble the diagnosed min slp from an analysis
#            xa_ms_lat = d_model['xa_lat'][idx_time]
#            xa_ms_lon = d_model['xa_lon'][idx_time]
#            # ---condition 1: find the nearest TCvital min slp from the analysis
#            # find the index/location of 'slp' in fort.10000
#            indices = [i for i ,e in enumerate(list_enkf_minSlp) if e == 'slp']
#            # assemble a pair of coordinate for each 'slp'
#            distances = []
#            obs_slp = []
#            for it in indices:
#                obs_slp.append( float(list_enkf_minSlp[it+9]) )
#                lon1 = float(list_enkf_minSlp[it+2])
#                lat1 = float(list_enkf_minSlp[it+1])
#                distances.append( UD.mercator_distance(lon1, lat1, xa_ms_lon, xa_ms_lat) )
#            min_distances = [np.amin(distances) == it for it in distances]
#            idx_nearest = min_distances.index(True)
#            # ---condition 2: the min slp
#            min_obs_slp = [np.amin(obs_slp) == it for it in obs_slp]
#            idx_min = min_obs_slp.index(True)
#            # ---combine the two conditions
#            if idx_min == idx_nearest:
#                idx_coor = idx_min
#                print('Storm center is choosed with two condtions met!')
#            else:
#                print('Not sure which obs is the storm center. Go with the min value one!')
#                idx_coor = idx_min
#            # gather this TCvital min slp
#            obs_minslp_lat.append( float(list_enkf_minSlp[indices[idx_coor]+1]) )
#            obs_minslp_lon.append( float(list_enkf_minSlp[indices[idx_coor]+2]) )
#            obs_minslp_value.append( float(list_enkf_minSlp[indices[idx_coor]+9])/100 )

    d_obs['lon_obs'] = obs_minslp_lon
    d_obs['lat_obs'] = obs_minslp_lat
    d_obs['minslp_obs'] = obs_minslp_value

    # Assemble the dictionary
    return d_obs

# Obtain minimum sea level pressure from best-track file (m/s)
def MSLP_btk( small_dir,Storm,DAt_6hrs ):

    Best_track_path = sorted(fnmatch.filter(os.listdir(small_dir+Storm+'/TC_Guidance/'),'bal*'))
    Best_track_file = small_dir+Storm+'/TC_Guidance/'+Best_track_path[0]
    with open(Best_track_file, 'r') as f:
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

    # MSLP and location 
    d_obs6Hr = {}
    d_obs6Hr['lat'] = {}
    d_obs6Hr['lon'] = {}
    d_obs6Hr['mslp'] = {}
    for i in idx_lines:
        line = btk_all[i]
        split_line = line.split()
        btk_time = split_line[2].replace(',','') + '00'
        # Read latitude
        lat_line = split_line[6].replace(',','')
        if 'N' in lat_line:
            d_obs6Hr['lat'][btk_time] = float(lat_line.replace('N',''))/10
        else:
            d_obs6Hr['lat'][btk_time] = 0-float(lat_line.replace('S',''))/10
        # Read longitute
        lon_line = split_line[7].replace(',','')
        if 'W' in lon_line:
            d_obs6Hr['lon'][btk_time] = 0-float(lon_line.replace('W',''))/10
        else:
            d_obs6Hr['lon'][btk_time] = float(lon_line.replace('E',''))/10
        # Read min sea level pressure
        d_obs6Hr['mslp'][btk_time] = float(split_line[9].replace(',','')) # mb

    return d_obs6Hr


def find_minSLP( Storm, wrfout, DAtime ):

    ncdir = nc.Dataset( wrfout )
    # geographical 
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # sea level pressure
    slp = UD.compute_slp( ncdir,NX=297,NY=297 )
    #slp = getvar(ncdir, 'slp')
    # original SLP
    slp_values = slp #slp.values
    slp_values[slp_values > 1030] = np.nan
    # smoothed SLP
    slp_smt_values = sp.ndimage.gaussian_filter(slp, [1,1]) #[11,11]
    slp_smt_values[slp_smt_values > 1030] = np.nan
    # simulated storm center
    if Storm == 'HARVEY': # Harvey is special with its location near land!

        # eyeball where the storm is
        if DAtime <= '201708221600':
            lat_mask = lat <= 18
            lon_mask = lon <= -91.5
            mask = lat_mask | lon_mask
        elif (DAtime >= '201708221700') & (DAtime <= '201708222300'):
            lat_mask = lat <= 18
            lon_mask =  (lon >= -88) | (lon <= -91.5)
            mask = lat_mask | lon_mask
        else:
            mask = lat <= 18
        slp_masked = ma.masked_array(slp_values, mask=mask)
        minslp = np.nanmin( slp_masked ) 
        
        slp_smooth_masked = ma.masked_array(slp_smt_values, mask=mask)
        idx = np.nanargmin( slp_smooth_masked )
        lat_minslp =  lat.flatten()[idx] 
        lon_minslp = lon.flatten()[idx] 
    else:
        minslp = np.nanmin( slp_values ) 
        idx = np.nanargmin( slp_smt_values )
        lat_minslp =  lat.flatten()[idx] 
        lon_minslp = lon.flatten()[idx]
    
    return [lat_minslp,lon_minslp,minslp]


# Obtain min slp from WRF output
def model_minSLP( big_dir,istorm,iExper,DAtimes,slp_xa,slp_xb=False ):

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
            list_xa = find_minSLP( istorm,xa_dir,DAtime )
            xa_minslp_lat.append( list_xa[0] )
            xa_minslp_lon.append( list_xa[1] )
            xa_minslp_value.append( list_xa[2] )

        if slp_xb:
            xb_dir = big_dir+istorm+'/'+iExper+'/fc/'+DAtime+'/wrf_enkf_input_d03_mean'
            print('Reading the EnKF prior mean from ', xb_dir)
            list_xb = find_minSLP( istorm,xb_dir,DAtime )
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


# Read pre-calculated HPI data
def read_precompute_HPI( save_dir ):

    model_hpi = {}
    # file paths
    file_xb = save_dir+'HPI_wrf_enkf_input_d03_mean.'+start_time_str[ist]+'_'+end_time_str[ist]+'.txt'
    file_xa = save_dir+'HPI_wrf_enkf_output_d03_mean.'+start_time_str[ist]+'_'+end_time_str[ist]+'.txt'
    files = [file_xb,file_xa]

    for file in files:
        mslp_lon = []
        mslp_lat = []
        mslp = []
        vmax = []
        with open(file,'r') as f:
            all_lines = f.readlines()[1:]
        for line in all_lines:
            if 'Success' in line:
                continue
            split_lines = line.split()
            mslp_lon.append( float(split_lines[1]) )
            mslp_lat.append( float(split_lines[2]) )
            mslp.append( float(split_lines[3]) )
            vmax.append( float(split_lines[6]) )
        # fill in the dictionary
        if 'input' in file:
            model_hpi['lon'] = mslp_lon
            model_hpi['lat'] = mslp_lat
            model_hpi['xb_mslp'] = mslp
            model_hpi['xb_vmax'] = vmax
        elif 'output' in file:
            model_hpi['lon'] = mslp_lon
            model_hpi['lat'] = mslp_lat
            model_hpi['xa_mslp'] = mslp
            model_hpi['xa_vmax'] = vmax

    return model_hpi

# ---------------------------------------------------------
#       Min SLP Only
# ---------------------------------------------------------

def plot_sys_minslp( d_model ):

    # Obtain mean absolute error
    MAE = {}
    for istorm in Storms:
        MAE[istorm] = {}
        for imp in MP:
            MAE[istorm][imp] = {}
            for ida in DA:
                MAE[istorm][imp][ida] = {}
                if slp_xb:
                    MAE[istorm][imp][ida]['xb'] = mean_absolute_error(d_tcvital[istorm]['minslp_obs'],d_model[istorm][imp][ida]['xb_slp']) 
                if slp_xa:
                    if Exper_names[ist][imp][ida] is None:
                         MAE[istorm][imp][ida]['xa'] = None
                    else:
                        MAE[istorm][imp][ida]['xa'] = mean_absolute_error(d_tcvital[istorm]['minslp_obs'],d_model[istorm][imp][ida]['xa_slp'])


    # Set up figure
    fig = plt.figure( figsize=(18,9), dpi=300 )
    ax = plt.subplot(1,1,1) 

    # Set colors
    #colorset = {'HARVEY':'#FF13EC','JOSE':'#0D13F6','MARIA':'#FFA500','IRMA':'#DB2824'}
    colorset = {'WSM6':'#FF13EC','THO':'#FFA500'}
    lead_t = list(range(0, cycles, 1))
   
    # Plot obs
    for istorm in Storms:
        ax.plot(lead_t, d_tcvital[istorm]['minslp_obs'],color='black',linestyle='-',linewidth=3)#,label=istorm)

    # plot slp_xa
    if slp_xa and not slp_xb:
        for istorm in Storms:
            for imp in MP:
                for ida in DA: # linestyle=(0, (5, 1))
                    if ida == 'CONV':
                        ax.plot(lead_t, d_model[istorm][imp][ida]['xa_slp'],color=colorset[imp],linestyle='-',linewidth=4,label=imp+':'+ida) # linestyle=(0, (5, 1))
                    elif ida == 'IR':
                        ax.plot(lead_t, d_model[istorm][imp][ida]['xa_slp'],color=colorset[imp],linestyle='--',linewidth=4,label=imp+':'+ida)
                    else:
                        ax.plot(lead_t, d_model[istorm][imp][ida]['xa_slp'],color=colorset[imp],linestyle=(0, (1, 1)),linewidth=4,label=imp+':'+ida)
                        #ax.scatter(lead_t, d_model[istorm][imp][ida]['xa_slp'],s=5,color=colorset[istorm],marker='*',markersize=12)
                        
    # plot saw-tooth lines
    if slp_xa and  slp_xb:
        for istorm in Storms:
            for imp in MP:
                for ida in DA: # linestyle=(0, (5, 1))
                    # zip data
                    t_zip = list( chain.from_iterable( zip(lead_t,lead_t)) )
                    len_seg = len( t_zip )
                    slp_zip = list( chain.from_iterable( zip(d_model[istorm][imp][ida]['xb_slp'],d_model[istorm][imp][ida]['xa_slp']) ) )

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
                            ax.plot(t_zip[i-1:i+1],slp_zip[i-1:i+1],linestyle=line,linewidth='6',color=colorset[imp])
                        elif i == 2:
                            ax.plot(t_zip[i-1:i+1],slp_zip[i-1:i+1],linestyle=line,linewidth='4',color=colorset[imp],label=ida )#istorm)
                        else:
                            ax.plot(t_zip[i-1:i+1],slp_zip[i-1:i+1],linestyle=line,linewidth='4',color=colorset[imp])


    # labels and attributes
    leg = ax.legend(loc='lower left',fontsize=25)
    # Set X/Y labels
    ax.set_xticks( lead_t[::4] )
    ax.set_xlim([0,cycles-1])
    ax.set_xlabel('EnKF Cycle',fontsize=28)
    ax.grid(True,linewidth=1, color='gray', alpha=0.5, linestyle='-')
    ax.set_ylabel('Minimum Sea Level Pressure (hPa)',fontsize=28)
    ax.set_ylim( 1000,1015 )
    #ax.set_ylim( 920,1000 )
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
 
    #if slp_xa and not slp_xb:
    #    fig.suptitle('TCvtial V.S. Analysis', fontsize=24, fontweight='bold')

    # Save the figure
    save = Storms[0]+'_'+MP[0]+'_'+DA[0]+'_'+DA[1]+'_'+DA[2]+'.png'
    des_name = small_dir+'Clean_results/SYSTEMS/Vis_analyze/Model/'+save
    plt.savefig( des_name )
    print( 'Saving the figure: ', des_name )
    plt.close()

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['JOSE']
    DA = ['CONV','IR','MW']
    MP = ['THO',] #

    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708231200','IRMA':'201709040000','JOSE':'201709060000','MARIA':'201709170000'}

    # observation
    if_tcvital = True
    if_btk = False

    slp_xa = True
    slp_xb = True
    Plot_minslp_compare = True
    cycles = 25
    multiple_segment = False

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

    # Read assimilated TC-vital min slp and its location 
    if if_tcvital:
        d_tcvital = {}
        for istorm in Storms:
            d_tcvital[istorm] = assimilated_obs( big_dir,small_dir,istorm,Exper_names[istorm][MP[0]]['IR'],dict_times[istorm],slp_xa)

    # Read best-track storm min slp
    if if_btk:
        d_btk = {}
        for istorm in Storms:
            d_btk6Hrs = MSLP_btk( small_dir,istorm,d_6hrs[istorm] )
            d_btkHr[istorm][iv] = Vx.Interpolate_hourly( istorm,d_6hrs,d_btk6Hrs['mslp'],d_hrs)

    # Read simulations
    d_model = {}
    for ist in Storms:
        d_model[ist] = {}
        for imp in MP:
            d_model[ist][imp] = {}
            for ida in DA:
                d_model[ist][imp][ida] = {}
                if not os.path.exists( big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida] ):
                    print('Path does NOT exist'+big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida])
                    d_model[ist][imp][ida] = None
                else:
                    # Read the precomputed data
                    precompute_dir = small_dir+'Clean_results/'+ist+'/'+Exper_names[ist][imp][ida]+'/Data_analyze/'
                    model_hpi = read_precompute_HPI( precompute_dir )
                    if slp_xb:
                        d_model[ist][imp][ida]['xb_slp'] = model_hpi['xb_mslp']
                        #print( d_model[ist][imp][ida]['xb_mslp']  )
                        #print( d_model['JOSE']['THO']['CONV']['xb_slp'] )
                    if slp_xa:
                        d_model[ist][imp][ida]['xa_slp'] = model_hpi['xa_mslp']



    # Compare analysis mean with TCvital
    if Plot_minslp_compare:
        print('Comparing min slp...')
        start_time=time.process_time()
        plot_sys_minslp( d_model )
        end_time = time.process_time()
        print('time needed: ', end_time-start_time, ' seconds')

    #if Plot_track_compare:
    #    print('Comparing track...')
    #    start_time=time.process_time()
    #    plot_track_compare()
    #    end_time = time.process_time()
    #    print('time needed: ', end_time-start_time, ' seconds')


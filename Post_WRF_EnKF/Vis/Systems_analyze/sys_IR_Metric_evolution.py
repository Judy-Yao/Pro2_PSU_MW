#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import subprocess
from itertools import chain

import Util_data as UD
import  Read_Obspace_IR as IR_obs

# setting font size
plt.rcParams.update({'font.size': 15})

# Generate time series
def generate_times( Storms, start_time_str, end_time_str ):

    dict_times = {}
    for istorm in Storms:
        time_diff = datetime.strptime(end_time_str[istorm],"%Y%m%d%H%M") - datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        dict_times[istorm] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times

def RMSE(simu, obs):
    return np.sqrt(np.nanmean((simu - obs) ** 2))

def Bias(simu, obs):
    return  np.nansum((simu - obs),0)/np.sum(~np.isnan(obs), 0)   # ignore nan values


def IR_metric( ist,ida,imp ):

    print(ist,ida,imp,Expers[ist][imp][ida])

    xb_bias = []
    xa_bias = []
    xb_rmse = []
    xa_rmse = []

    for DAtime in DAtimes[ist]:
        Tb_file = big_dir+ist+'/'+Expers[ist][imp][ida]+'/Obs_Hx/IR/'+DAtime+"/mean_obs_res_d03_"+DAtime+'_'+ sensor+'.txt'
        if os.path.isfile( Tb_file ):
            d_all = IR_obs.read_Tb_obsRes(Tb_file, sensor )
        else:
            raise ValueError(Tb_file+' does not exist!')

        # limitation
        if limit:
            condi = d_all['Yo_obs'] <= 220
            idx_lm = np.where( condi )[0]
        else:
            idx_lm = range(len(d_all['Yo_obs'])) 
        
        # Bias
        xb_bias.append( Bias(d_all['meanYb_obs'][idx_lm], d_all['Yo_obs'][idx_lm] ))
        xa_bias.append( Bias(d_all['meanYa_obs'][idx_lm], d_all['Yo_obs'][idx_lm] ))

        # RMSE
        xb_rmse.append( RMSE(d_all['meanYb_obs'][idx_lm], d_all['Yo_obs'][idx_lm] ))
        xa_rmse.append( RMSE(d_all['meanYa_obs'][idx_lm], d_all['Yo_obs'][idx_lm] ))

    # Assemble the dictionary
    d_metric = {'xb_rmse': xb_rmse,'xa_rmse':xa_rmse,'xb_bias':xb_bias,'xa_bias':xa_bias }
    return d_metric

# Sum up the domain-mean bias over all cycles across all storms
# Calculate the mean
def mean_bias( d_bias ):
    
    # Collect biases
    collect_bias = []

    for ist in Storms:
        for imp in MP:
            for ida in DA:
                collect_bias.append( d_bias[ist][imp][ida]['xb'] )

    # Calculate the number of cycles
    num_cycle = 0
    for ist in Storms:
        num_cycle += len(DAtimes[ist])

    # Average 
    return collect_bias, num_cycle

if __name__ == '__main__':


    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    DA = ['IR']
    MP = ['THO',]

    sensor = 'abi_gr'
    ch_list = ['8',]

    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708231200','IRMA':'201709040000','JOSE':'201709060000','MARIA':'201709170000'}
    Consecutive_times = True

    # limitations
    limit = False

    if_bias = True
    if_rmse = False
    # ------------------------------------------------------   

    # Create experiment names
    Expers = {}
    for ist in Storms:
        Expers[ist] = {}
        for imp in MP:
            Expers[ist][imp] = {}
            for ida in DA:
                Expers[ist][imp][ida] = UD.generate_one_name( ist,ida,imp )

    # Identify DA times in the period of interest
    DAtimes = generate_times( Storms, start_time_str, end_time_str )

    # Calculate metrics
    d_bias = {}
    d_rmse = {}
    for ist in Storms:
        d_bias[ist] = {}
        d_rmse[ist] = {}
        for imp in MP:
            d_bias[ist][imp] = {}
            d_rmse[ist][imp] = {}
            for ida in DA:
                d_bias[ist][imp][ida] = {}
                d_rmse[ist][imp][ida] = {}
                # calculate metrics
                d_metric = IR_metric( ist,ida,imp )
                d_bias[ist][imp][ida]['xb'] = d_metric['xb_bias']
                d_bias[ist][imp][ida]['xa'] = d_metric['xa_bias']
                d_rmse[ist][imp][ida]['xb'] = d_metric['xb_rmse']
                d_rmse[ist][imp][ida]['xa'] = d_metric['xa_rmse']

    # Calculate the average bias
    collect_bias, num_cycle = mean_bias( d_bias )
    mean_bias = np.sum(collect_bias)/num_cycle
    print('For '+str(num_cycle)+' cycles:')
    print('the averaged value for the domain-mean biases is '+str(mean_bias))


    # Plot the metrics
    #if if_bias:
    #    Plot_bias( d_bias ) 
    #if if_rmse:
    #    Plot_rmse( d_rmse )


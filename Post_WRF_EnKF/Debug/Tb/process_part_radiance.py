#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import time

from Track_xbxa import read_HPI_model
import Util_data as UD
import Read_Obspace_IR as ROIR
import Diagnostics as Diag


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    model_resolution = 3000 #m

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR'
    MP = 'WSM6'

    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709030700'
    end_time_str = '201709030700'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    deep_slp_incre = True
    incre_slp_th = 20
    plot_circle = True
    radius_th = 200 # km

    to_obs_res = True

    outputIR = True
    If_save = True
    # -------------------------------------------------------    

    # Create experiment names
    #Exper_name = UD.generate_one_name( Storm,DA,MP )
    Exper_name = 'IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'

    if not Consecutive_times:
        IR_times = []
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Find the value of minimum slp increment
    if deep_slp_incre:
        # ----- Read min slp from model-----------------
        HPI_models = {}
        DAtimes_dir = [big_dir+Storm+'/'+Exper_name+'/fc/'+it for it in IR_times]
        file_kinds = ['wrf_enkf_input_d03_mean','wrf_enkf_output_d03_mean']
        for ifk in file_kinds:
            idx = file_kinds.index( ifk )
            HPI_models[ifk] = read_HPI_model( Storm, Exper_name, ifk, DAtimes_dir )
        incre_slp = np.array(HPI_models['wrf_enkf_output_d03_mean']['min_slp']) - np.array(HPI_models['wrf_enkf_input_d03_mean']['min_slp'])

    # output IR obs of interest
    if outputIR and deep_slp_incre:
        for DAtime in IR_times:
            idx_t = IR_times.index( DAtime )
            if abs(incre_slp[idx_t]) > incre_slp_th:
                print('At '+DAtime)
                # Find IR of interest
                Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
                ct_lon = HPI_models['wrf_enkf_output_d03_mean']['lon'][idx_t]
                ct_lat = HPI_models['wrf_enkf_output_d03_mean']['lat'][idx_t]
                Tb_file = Hx_dir + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
                d_all = ROIR.read_Tb_obsRes(Tb_file, sensor )
                idx_hxb = UD.find_circle_area_ungrid( ct_lon, ct_lat, d_all['lon_obs'], d_all['lat_obs'], radius_th )
                # Find innovation (yo-Hxb) > 0
                inno = d_all['Yo_obs']-d_all['meanYb_obs']
                idx_2 = np.where(inno<0)[0]
                idx_3 = np.where(d_all['lat_obs']<=18)[0]
                idx_inter1 = np.intersect1d(idx_hxb,idx_2)
                idx_inter = np.intersect1d(idx_inter1,idx_3)
                # Name a file
                file_name = 'part_radiance_'+DAtime+'_so'
                # Define the format specifier
                formatSpec = '{:>12s}{:>12s}{:>12d}{:>12.3f}{:>12.3f}{:>12.3f}{:>12d}{:>12d}{:>12.3f}{:>12.3f}\n'
                # Write these IR obs to a file
                with open(file_name,'w') as f:
                    # Write the record to the file serially
                    for it in idx_inter:
                        data = [DAtime,"abi_gr",8,d_all['lat_obs'][it],d_all['lon_obs'][it],d_all['Yo_obs'][it],0,200,3,35000]
                        # Write the record using the format specifier
                        formatted_record = formatSpec.format(*data)
                        f.write(formatted_record)





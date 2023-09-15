#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os
import glob
import numba
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import math
from datetime import datetime, timedelta
import time

import Util_data as UD
import Util_Vis
#import matlab.engine
import Diagnostics as Diag

# numba histogram
#@numba.njit
#def hist1d(v, b, r):
#    return np.histogram(v, b, r)[0]
    

# Read variables at obs resolution/location
def read_Tbs_obsRes(Tb_files, sensor ):

    dict_allTb = {}

    if bin_Tbdiff:
        diff_ob_all = []
        diff_oa_all = []
    else:
        Yo_all = []
        meanYb_all = []
        meanYa_all = []
 
    for DAtime in IR_times:

        Tb_file = Hx_dir+DAtime+'/mean_obs_res_d03_' + DAtime + '_' +  sensor + '.txt'
        if bin_Tbdiff:
            diff_ob = []
            diff_oa = [] 
        else:
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
            if bin_Tbdiff:
                diff_ob.append( float(split_line[4])-float(split_line[3]) ) # Hxb - Yo 
                diff_oa.append( float(split_line[5])-float(split_line[3]) ) # Hxa - Yo
            else:
                Yo_obs.append( float(split_line[3]) )
                meanYb_obs.append( float(split_line[4]) )
                meanYa_obs.append( float(split_line[5]) )

        # add up records
        if bin_Tbdiff:
            diff_ob_all.extend( diff_ob )
            diff_oa_all.extend( diff_oa )
        else:
            Yo_all.extend( Yo_obs )
            meanYb_all.extend( meanYb_obs )
            meanYa_all.extend( meanYa_obs )       

        # convert to array
        if bin_Tbdiff:
            diff_ob = np.array( diff_ob )
            diff_oa = np.array( diff_oa )
            dict_allTb[DAtime] = {'diff_ob':diff_ob, 'diff_oa':diff_oa}
        else:
            Yo_obs = np.array( Yo_obs )
            meanYb_obs = np.array( meanYb_obs )
            meanYa_obs = np.array( meanYa_obs ) 
            dict_allTb[DAtime] = {'Yo_obs':Yo_obs, 'meanYb_obs':meanYb_obs, 'meanYa_obs':meanYa_obs}

    if bin_Tbdiff:
        diff_ob_all = np.array( diff_ob_all )
        diff_oa_all = np.array( diff_oa_all )
        dict_allTb['All'] = {'diff_ob':diff_ob_all, 'diff_oa':diff_oa_all}
    else:
        Yo_all = np.array( Yo_all )
        meanYb_all = np.array( meanYb_all )
        meanYa_all = np.array( meanYa_all )
        dict_allTb['All'] = {'Yo_obs':Yo_all, 'meanYb_obs':meanYb_all, 'meanYa_obs':meanYa_all}

    return dict_allTb

# Plot histograms
def Plot_IR_hist( d_hcount ):


    fig,axs = plt.subplots(1,1, figsize=(8,8), dpi=300 )

    warm_color = ["#c23728","#e14b31","#de6e56","#e1a692","#786028","#a57c1b","#d2980d","#ffb400","#503f3f","#6d4b4b","#a86464","#e27c7c"] #redish
    cold_color = ["#115f9a", "#1984c5", "#22a7f0", "#48b5c4", "#48446e", "#5e569b", "#776bcd", "#9080ff","#3c4e4b", "#466964", "#599e94", "#6cd4c5"] #blueish

    x_axis = (range_bins[:-1]+range_bins[1:])/2
    idx = 0

    if bin_Tbdiff:
        for outkey in d_hcount:
            if outkey == 'All':
                for key in d_hcount[outkey]:
                    if key == 'diff_ob':
                        axs.plot(x_axis,d_hcount[outkey][key],color='red',linestyle='-',linewidth='4',label='H(Xb)-Obs')
                    elif key == 'diff_oa':
                        axs.plot(x_axis,d_hcount[outkey][key],color='blue',linestyle='-',linewidth='4',label='H(Xa)-Obs')
            else:
                for key in d_hcount[outkey]:
                    if key == 'diff_ob':
                        axs.plot(x_axis,d_hcount[outkey][key],color="#de6e56",linestyle='-',linewidth='2')
                        #axs.plot(x_axis,d_hcount[outkey][key],color=warm_color[idx],linestyle='-',linewidth='2')
                    elif key == 'diff_oa':
                        axs.plot(x_axis,d_hcount[outkey][key],color="#22a7f0",linestyle='-',linewidth='2')
                        #axs.plot(x_axis,d_hcount[outkey][key],color=cold_color[idx],linestyle='-',linewidth='2')
                idx = idx+1
    else:
        for outkey in d_hcount:
            if outkey == 'All':
                for key in d_hcount[outkey]:
                    if key == 'Yo_obs': 
                        axs.plot(x_axis,d_hcount[outkey][key],color='black',linestyle='-',linewidth='4',label='Obs:all cycles')
                    elif key == 'meanYb_obs':
                        axs.plot(x_axis,d_hcount[outkey][key],color='red',linestyle='-',linewidth='4',label='H(Xb):all cycles')
                    elif key == 'meanYa_obs':
                        axs.plot(x_axis,d_hcount[outkey][key],color='blue',linestyle='-',linewidth='4',label='H(Xa):all cycles')
            else:
                for key in d_hcount[outkey]:
                    if key == 'Yo_obs':
                        pass
                        #axs.plot(x_axis,d_hcount[outkey][key],color='black',linestyle='--',linewidth='0.5') 
                    elif key == 'meanYb_obs':
                        axs.plot(x_axis,d_hcount[outkey][key],color="#de6e56",linestyle='-',linewidth='2')
                        #axs.plot(x_axis,d_hcount[outkey][key],color=warm_color[idx],linestyle='-',linewidth='2')
                    elif key == 'meanYa_obs':
                        axs.plot(x_axis,d_hcount[outkey][key],color="#22a7f0",linestyle='-',linewidth='2')
                        #axs.plot(x_axis,d_hcount[outkey][key],color=cold_color[idx],linestyle='-',linewidth='2')
                    else:
                        pass
                idx = idx+1




    axs.legend(frameon=True,loc='upper right',fontsize='14')
    axs.grid(True,linestyle='--',alpha=0.5)
    axs.set_title( DA+'_'+MP,fontweight="bold",fontsize='15' )
    
    axs.tick_params(axis='x', labelsize=15)
    axs.tick_params(axis='y', labelsize=15)
    if bin_Tbdiff:
        axs.set_ylim(ymin=0,ymax=0.25)
        axs.set_xlim(xmin=min_Tbdiff_rg ,xmax=max_Tbdiff_rg )
    else:
        axs.set_ylim(ymin=0,ymax=0.12)
        axs.set_xlim(xmin=min_Tb_rg,xmax=max_Tb_rg)
    axs.set_xlabel('Brightness Temperature (K)',fontsize=15)
    axs.set_ylabel('Density',fontsize=15)
    
    if bin_Tbdiff:
        fig.suptitle( Storm+': PDF of Bias of IR Tb over '+str(len(IR_times))+' Cycles', fontsize=15, fontweight='bold')
        des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IR_PDF_'+IR_times[0]+'_'+IR_times[-1]+'_bias.png'
    else:
        fig.suptitle( Storm+': PDF of IR Tbs over '+str(len(IR_times))+' Cycles', fontsize=15, fontweight='bold')
        des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IR_PDF_'+IR_times[0]+'_'+IR_times[-1]+'.png'
    
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name+'!' )



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR'
    MP = 'THO'

    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709030000'
    end_time_str = '201709031100'
    Consecutive_times = True

    number_bins = 50

    min_Tb_rg = 190
    max_Tb_rg = 260
    min_Tbdiff_rg = -25
    max_Tbdiff_rg = 25

    Make_bins = True
    bin_Tbdiff = True
    
    If_plot = True
    # ------------------------------------------------------   

    # Create experiment names
    Exper_name =  UD.generate_one_name( Storm,DA,MP )

    if not Consecutive_times:
        IR_times = ['201709030100',]
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]


    # Make bins and count the number in bins
    if Make_bins:
        # Read obs, Hxb, and Hxa of all files 
        Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'
        d_allTb =  read_Tbs_obsRes(Hx_dir, sensor )
        # Make bins
        print('Counting number per bins...')
        start_time=time.process_time()
        d_hcount = {}
        for key in d_allTb: # initialize an empty nested dictionary
            d_hcount[key] = {} 

        for outer_key in d_allTb:
            print(outer_key)
            for inner_key in d_allTb[outer_key]:
                print(inner_key)
                if bin_Tbdiff:
                    hist,range_bins = np.histogram(d_allTb[outer_key][inner_key], number_bins, (min_Tbdiff_rg,max_Tbdiff_rg),density=True )
                else:
                    hist,range_bins = np.histogram(d_allTb[outer_key][inner_key], number_bins, (min_Tb_rg,max_Tb_rg),density=True )
                d_hcount[outer_key][inner_key] = hist
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds') 
    
        if If_plot:
             Plot_IR_hist( d_hcount )




    








































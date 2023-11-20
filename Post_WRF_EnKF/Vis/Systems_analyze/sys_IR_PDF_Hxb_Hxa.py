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


# Generate time series
def generate_times( Storms, start_time_str, end_time_str ):

    dict_times = {}
    for istorm in Storms:
        time_diff = datetime.strptime(end_time_str[istorm],"%Y%m%d%H%M") - datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        dict_times[istorm] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times


# Read variables at obs resolution/location for one experiment 
def read_Tbs_obsRes_oneExper(istorm,imp,ida,Exper_names,d_times,sensor):

    Hx_dir = big_dir+istorm+'/'+Exper_names[istorm][imp][ida]+'/Obs_Hx/IR/'
    dict_allTb = {}

    if bin_Tbdiff:
        diff_ob_all = []
        diff_oa_all = []
    else:
        Yo_all = []
        meanYb_all = []
        meanYa_all = []

    for DAtime in d_times[istorm]:

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
                diff_oa.append( float(split_line[5])-float(split_line[3])) # Hxa - Yo
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

    if bin_Tbdiff:
        #diff_ob_all = np.array( diff_ob_all )
        #diff_oa_all = np.array( diff_oa_all )
        #print(np.shape(diff_oa_all))
        dict_allTb['All_times'] = {'diff_ob':diff_ob_all, 'diff_oa':diff_oa_all}
    else:
        #Yo_all = np.array( Yo_all )
        #meanYb_all = np.array( meanYb_all )
        #meanYa_all = np.array( meanYa_all )
        dict_allTb['All_times'] = {'Yo_obs':Yo_all, 'meanYb_obs':meanYb_all, 'meanYa_obs':meanYa_all}

    return dict_allTb

def combine_storms_allTimes( Exper_Tb ):

    allStorm_tb = {}
    if bin_Tbdiff:
        all_diff_ob = []
        all_diff_oa = []
    else:
        all_Yo = []
        all_meanYb = []
        all_meanYa = []

    for istorm in Storms:
        tb_case = Exper_Tb[istorm][imp][ida]['All_times']
        if bin_Tbdiff:
            all_diff_ob.extend( tb_case['diff_ob'] )
            all_diff_oa.extend( tb_case['diff_oa'] )
        else:
            all_Yo.extend( tb_case['Yo_obs'] )
            all_meanYb.extend( tb_case['meanYb_obs'] )
            all_meanYa.extend( tb_case['meanYa_obs'] )

    if bin_Tbdiff:
        all_diff_ob = np.array( all_diff_ob )
        all_diff_oa = np.array( all_diff_oa )
        allStorm_tb['All_times'] = {'diff_ob':all_diff_ob, 'diff_oa':all_diff_oa}
    else:
        all_Yo = np.array( all_Yo )
        all_meanYb = np.array( all_meanYb )
        all_meanYa = np.array( all_meanYa )
        allStorm_tb['All_times'] = {'Yo_obs':all_Yo, 'meanYb_obs':all_meanYb, 'meanYa_obs':all_meanYa} 

    return allStorm_tb

def csm_color(imp,ida):
    if ida == 'conv' and imp == 'WSM6':
        return '#FFED13'
    if ida == 'IR' and imp == 'WSM6':
        return '#DB2824'
    elif ida == 'IR+MW' and imp == 'WSM6':
        return '#FF13EC'
    elif ida == 'conv' and imp == 'THO':
        return '#B1D866'
    elif ida == 'IR' and imp == 'THO':
        return '#0D13F6'
    elif ida == 'IR+MW' and imp == 'THO':
        return '#684BA0'


# Plot the PDF
def Plot_hist_IRsum( d_hcount ):

    fig,axs = plt.subplots(1,1, figsize=(10.3,8), dpi=300 )

    # customize for color
    line_style = {}
    line_style['diff_ob'] = '-'
    line_style['diff_oa'] = (0, (5, 1)) # densely dashed

    x_axis = (range_bins[:-1]+range_bins[1:])/2
    idx = 0

    if bin_Tbdiff:
        for imp in MP:
            for ida in DA:
                for key in d_hcount[imp][ida]:
                    if 'b' in key:
                        hx=r'$\mathbf{\overline{H(Xb)}}$'
                    else:
                        hx=r'$\mathbf{\overline{H(Xa)}}$'
                    axs.plot(x_axis,d_hcount[imp][ida][key],color=csm_color(imp,ida),linestyle=line_style[key],linewidth='6',label=ida+'_'+imp+':'+hx) 
    else:
        for imp in MP:
            for ida in DA:
                for key in d_hcount[imp][ida]:
                    axs.plot(x_axis,d_hcount[imp][ida][key],color=csm_color(imp,ida),linestyle=line_style[key],linewidth='6',label=ida+'_'+imp)

    # Shrink current axis by 20%
    box = axs.get_position()
    axs.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    axs.legend(loc='center left', frameon=True,bbox_to_anchor=(0.7, 0.5),fontsize='20')
    axs.grid(True,linestyle='--',alpha=0.5)
    axs.set_title( 'IR Tb Bias ('+r'$\mathbf{\overline{H(X)}}$'+'- Obs) of All Storms',fontweight="bold",fontsize='18' )

    axs.tick_params(axis='x', labelsize=24)
    axs.tick_params(axis='y', labelsize=20)
    if bin_Tbdiff:
        axs.set_ylim(ymin=0,ymax=0.30)
        axs.set_xlim(xmin=min_Tbdiff_rg ,xmax=max_Tbdiff_rg )
    else:
        axs.set_ylim(ymin=0,ymax=0.05)
        axs.set_xlim(xmin=min_Tb_rg,xmax=max_Tb_rg)
    axs.set_xlabel('Brightness Temperature (K)',fontsize=24)
    axs.set_ylabel('Density',fontsize=24)

    if bin_Tbdiff:
        fig.suptitle( 'PDF over '+str(len(dict_times[Storms[0]]))+' Cycles', fontsize=15, fontweight='bold')
        des_name = small_dir+'SYSTEMS/Vis_analyze/Tb/IR_PDF_'+str(len(dict_times[Storms[0]]))+'cycles_bias_multi_IR_WSM6_THO.png'
    else:
        fig.suptitle( 'PDF over '+str(len(IR_times))+' Cycles', fontsize=15, fontweight='bold')
        des_name = small_dir+'SYSTEMS/Vis_analyze/Tb/IR_PDF_'+str(len(IR_times))+'cycles_multi.png'

    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name+'!' )



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['IRMA','JOSE','MARIA']#['HARVEY','IRMA','JOSE','MARIA']
    DA = ['IR',]
    MP = ['THO','WSM6'] 
    sensor = 'abi_gr'

    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708231200','IRMA':'201709040000','JOSE':'201709060000','MARIA':'201709170000'}
    Consecutive_times = True

    number_bins = 50

    min_Tb_rg = 190
    max_Tb_rg = 260
    min_Tbdiff_rg = -25
    max_Tbdiff_rg = 25

    Make_bins = True
    bin_Tbdiff = True
 
    If_plot = True
    #------------------------------------

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

    # Number of kinds of experiments
    num_kinds = len(DA)*len(MP)

    # Read obs, Hxb, and Hxa of all files
    Exper_Tb = {}
    for istorm in Storms:
        Exper_Tb[istorm] = {}
        for imp in MP:
            Exper_Tb[istorm][imp] = {}
            for ida in DA:
                iExper = Exper_names[istorm][imp][ida]
                if iExper is not None:
                    Exper_Tb[istorm][imp][ida] = read_Tbs_obsRes_oneExper(istorm,imp,ida,Exper_names,dict_times,sensor) 
                else:
                    Exper_Tb[istorm][imp][ida] = None
    
    # Combine data for all storms
    Storms_Tb = {}
    for imp in MP:
        Storms_Tb[imp] = {}
        for ida in DA:
            Storms_Tb[imp][ida] = combine_storms_allTimes( Exper_Tb ) 

    # Make bins based on experiment types
    start_time=time.process_time()
    d_hcount = {}
    for imp in MP:
        d_hcount[imp] = {}
        for ida in DA:
            d_hcount[imp][ida] = {}
            for inner_key in Storms_Tb[imp][ida]['All_times']:
                if bin_Tbdiff:
                    hist,range_bins = np.histogram(Storms_Tb[imp][ida]['All_times'][inner_key], number_bins, (min_Tbdiff_rg,max_Tbdiff_rg),density=True )
                else:
                    hist,range_bins = np.histogram(Storms_Tb[imp][ida]['All_times'][inner_key], number_bins, (min_Tb_rg,max_Tb_rg),density=True )
                d_hcount[imp][ida][inner_key] = hist
    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    if If_plot:
        Plot_hist_IRsum( d_hcount )















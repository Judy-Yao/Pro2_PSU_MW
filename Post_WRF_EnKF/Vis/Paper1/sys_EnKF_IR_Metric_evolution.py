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
import Obspace_compare_IR_txt_bin as IR_obs

# setting font size
#plt.rcParams.update({'font.size': 15})

# Generate time series
def generate_times( Storms, start_time_str, end_time_str, interval ):

    dict_times = {}
    for istorm in Storms:
        time_diff = datetime.strptime(end_time_str[istorm],"%Y%m%d%H%M") - datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+interval, interval))]
        dict_times[istorm] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times


def RMSE(simu, obs):
    return np.sqrt(np.nanmean((simu - obs) ** 2))

def Bias(simu, obs):
    return  np.nansum((simu - obs),0)/np.sum(~np.isnan(obs), 0)   # ignore nan values


# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot the evolution of IR metrics
# ------------------------------------------------------------------------------------------------------
# Bias over time
def Plot_rmse( d_rmse,ida ):

    # ------ Plot Figure -------------------
    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=300)
    grid = fig.add_gridspec(ncols=1,nrows=4,top=0.95,right=0.95,)
    ax = {}
    for ist in Storms:
        ax[ist] = {}
        ir = Storms.index( ist )
        # gridspec inside gridspec
        ax[ist] = fig.add_subplot( grid[ir] )

    # Customize labels
    labels = {}
    for imp in MP:
        if imp == 'THO':
            labels[imp] = {'xa':'THO_DA','xb':'THO_FC'}
        elif imp == 'WSM6':
            labels[imp] = {'xa':'WSM6_DA','xb':'WSM6_FC'}

    # Customize colors
    colors = {'WSM6':'red','THO':'blue'}

    # Customize vmin and vmax
    vmin = {'CONV':5,'IR':2}
    vmax = {'CONV':15,'IR':8.1}

    for ist in Storms:
        # Obtain dates for a storm
        dates = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes[ist]]
        dates_zip = list( chain.from_iterable( zip(dates,dates)) )
        for imp in MP:
            for ida in DA:
                rmse_zip = list( chain.from_iterable( zip(d_rmse[ist][imp][ida]['xb'],d_rmse[ist][imp][ida]['xa']) ) )
                len_seg = len(dates_zip)
                for i in range(1,len_seg):
                    # specify which segment uses which color/linestyle
                    if i % 2 == 0:
                        line='-'
                    else:
                        line='--'
                    # plot segment
                    if i == 1:
                        ax[ist].plot(dates_zip[i-1:i+1],rmse_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='4',)
                    elif i == 2:
                        if ist == Storms[0]:
                            ax[ist].plot(dates_zip[i-1:i+1],rmse_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='2',label=labels[imp]['xb']) #WSM6_Forecast
                        else:
                            ax[ist].plot(dates_zip[i-1:i+1],rmse_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='2')
                    elif i == 3:
                        if ist == Storms[0]:
                            ax[ist].plot(dates_zip[i-1:i+1],rmse_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='2',label=labels[imp]['xa'])
                        else:
                            ax[ist].plot(dates_zip[i-1:i+1],rmse_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='2')
                    else:
                        ax[ist].plot(dates_zip[i-1:i+1],rmse_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='2')
                # add patch: to get a sense where metrics of IR experiments are
                if 'CONV' in ida and imp == MP[0]:
                    ax[ist].axhspan(vmin['IR'], vmax['IR'], color='grey', alpha=0.2,)

        # axis
        ax[ist].set_xlim([dates[0],dates[-1]])
        ax[ist].set_ylim( vmin[ida],vmax[ida] )
        ax[ist].tick_params(axis='x', labelrotation=12,labelsize=8)
        ax[ist].tick_params(axis='y',labelsize=10)
        ax[ist].axhline(y=0.0,color='k',linestyle='-',linewidth='2')
        ax[ist].grid(True,linestyle='--',alpha=0.5)
        if ist == Storms[0]:
            ax[ist].legend(frameon=True,loc='upper center',fontsize='10',ncol=4) #24
        ax[ist].set_ylabel('IR RMSE (K)',fontsize=10)

    # title
    ax[Storms[0]].set_title( 'IR RMSE of WRF-EnKF Cyclings for '+DA[0]+' Experiments',fontweight="bold",fontsize='12' )
    #ax[Storms[0]].set_title( 'IR Bias: domain-mean{'+r'$\mathbf{\overline{H(X)}}$'+' - Obs}',fontweight="bold",fontsize='12' )

    # Add storm info
    for ist in Storms:
        if Storms.index(ist) == 0:
            fig.text(0.02,0.86,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 1:
            fig.text(0.02,0.64,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 2:
            fig.text(0.02,0.42,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 3:
            fig.text(0.02,0.20,ist, fontsize=12, ha='center', va='center',rotation='vertical')

    #fig.suptitle( suptt+DA[0]+' experiments', fontsize=15, fontweight='bold')

    # title
    #ax.set_title( 'IR Bias: domain-mean{'+r'$\mathbf{\overline{H(X)}}$'+' - Obs}',fontweight="bold",fontsize='15' )
    # saved name
    des_name = small_dir+'Clean_results/SYSTEMS/Vis_analyze/Paper1/IR_RMSE_'+ida+'_exps.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )



# Bias over time
def Plot_bias( d_bias,ida ):

    # ------ Plot Figure -------------------
    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=300)
    grid = fig.add_gridspec(ncols=1,nrows=4,top=0.95,right=0.95,)
    ax = {}
    for ist in Storms:
        ax[ist] = {}
        ir = Storms.index( ist )
        # gridspec inside gridspec
        ax[ist] = fig.add_subplot( grid[ir] )

    # Customize labels
    labels = {}
    for imp in MP:
        if imp == 'THO':
            labels[imp] = {'xa':'THO_DA','xb':'THO_FC'}
        elif imp == 'WSM6':
            labels[imp] = {'xa':'WSM6_DA','xb':'WSM6_FC'}

    # Customize colors
    colors = {'WSM6':'red','THO':'blue'}

    # Customize vmin and vmax
    vmin = {'CONV':-10,'IR':-4.1}
    vmax = {'CONV':10,'IR':4.1}
    # add patch: to get a sense where metrics of IR experiments are
    #if 'CONV' in ida and imp == MP[0]:
    #    ax.axhspan(-4, 4, color='grey', alpha=0.2, label="IR Area")

    for ist in Storms:
        # Obtain dates for a storm
        dates = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes[ist]]
        #sitart_time = datetime.strptime( DAtimes[ist][0],"%Y%m%d%H%M")
        #end_time = datetime.strptime( DAtimes[ist][-1],"%Y%m%d%H%M")
        dates_zip = list( chain.from_iterable( zip(dates,dates)) )
        for imp in MP:
            for ida in DA:
                bias_zip = list( chain.from_iterable( zip(d_bias[ist][imp][ida]['xb'],d_bias[ist][imp][ida]['xa']) ) )
                len_seg = len(dates_zip)
                for i in range(1,len_seg):
                    # specify which segment uses which color/linestyle
                    if i % 2 == 0:
                        line='-'  
                    else:
                        line='--' 
                    # plot segment
                    if i == 1:
                        ax[ist].plot(dates_zip[i-1:i+1],bias_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='4',)
                    elif i == 2:
                        if ist == Storms[0]:
                            ax[ist].plot(dates_zip[i-1:i+1],bias_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='2',label=labels[imp]['xb']) #WSM6_Forecast
                        else:
                            ax[ist].plot(dates_zip[i-1:i+1],bias_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='2')
                    elif i == 3:
                        if ist == Storms[0]:
                            ax[ist].plot(dates_zip[i-1:i+1],bias_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='2',label=labels[imp]['xa'])
                        else:
                            ax[ist].plot(dates_zip[i-1:i+1],bias_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='2')
                    else:
                        ax[ist].plot(dates_zip[i-1:i+1],bias_zip[i-1:i+1],colors[imp],linestyle=line,linewidth='2')
                # add patch: to get a sense where metrics of IR experiments are
                if 'CONV' in ida and imp == MP[0]:
                    ax[ist].axhspan(vmin['IR'], vmax['IR'], color='grey', alpha=0.2,)

        # axis
        ax[ist].set_xlim([dates[0],dates[-1]])
        ax[ist].set_ylim( vmin[ida],vmax[ida] )
        ax[ist].tick_params(axis='x', labelrotation=12,labelsize=8)
        ax[ist].tick_params(axis='y',labelsize=10)
        ax[ist].axhline(y=0.0,color='k',linestyle='-',linewidth='2')
        ax[ist].grid(True,linestyle='--',alpha=0.5)
        if ist == Storms[0]:
            ax[ist].legend(frameon=True,loc='upper center',fontsize='10',ncol=4) #24
        ax[ist].set_ylabel('IR Bias (K)',fontsize=10)

    # title
    ax[Storms[0]].set_title( 'IR Bias of WRF-EnKF Cyclings for '+DA[0]+' Experiments',fontsize='12' )
    #ax[Storms[0]].set_title( 'IR Bias: domain-mean{'+r'$\mathbf{\overline{H(X)}}$'+' - Obs}',fontweight="bold",fontsize='12' )

    # Add storm info
    for ist in Storms:
        if Storms.index(ist) == 0:
            fig.text(0.02,0.86,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 1:
            fig.text(0.02,0.64,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 2:
            fig.text(0.02,0.42,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 3:
            fig.text(0.02,0.20,ist, fontsize=12, ha='center', va='center',rotation='vertical')
    
    #fig.suptitle( suptt+DA[0]+' experiments', fontsize=15, fontweight='bold')
    
    # title
    #ax.set_title( 'IR Bias: domain-mean{'+r'$\mathbf{\overline{H(X)}}$'+' - Obs}',fontweight="bold",fontsize='15' )
    # saved name
    des_name = small_dir+'Clean_results/SYSTEMS/Vis_analyze/Paper1/IR_Bias_'+ida+'_exps.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


def IR_metric( ist,ida,imp ):

    print('Calculating...')
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
    d_metric = {'xb_rmse': np.array(xb_rmse),'xa_rmse':np.array(xa_rmse),'xb_bias': np.array(xb_bias),'xa_bias':np.array(xa_bias) }
    return d_metric



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    MP = ['WSM6','THO']
    DA = ['IR',]

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
    for istorm in Storms:
        Expers[istorm] = {}
        for imp in MP:
            Expers[istorm][imp] = {}
            for ida in DA:
                Expers[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # DAtimes
    DAtimes = generate_times( Storms, start_time_str, end_time_str, 1 )

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

    # Plot the metrics
    if if_bias:
        Plot_bias( d_bias, DA[0] ) 
    if if_rmse:
        Plot_rmse( d_rmse, DA[0] )


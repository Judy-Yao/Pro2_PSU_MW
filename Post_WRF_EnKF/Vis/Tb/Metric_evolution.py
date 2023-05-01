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
import Obspace_compare_IR_txt_bin as IR_obs

# setting font size
plt.rcParams.update({'font.size': 15})

def RMSE(simu, obs):
    return np.sqrt( ((simu - obs) ** 2).mean() )

def Bias(simu, obs):
    return  np.sum((simu - obs),0)/np.size(obs,0)   


# ------------------------------------------------------------------------------------------------------
#           Operation: Read, process, and plot the evolution of IR metrics
# ------------------------------------------------------------------------------------------------------

def Plot_IR_metric( small_dir, Storm, Expers, DAtimes, IR_metric ):

    # ------ Define range of metrics -------------------
    rmse_min = 0#np.amin( [np.amin(IR_metric['xb_rmse']),np.amin(IR_metric['xa_rmse']) ])
    rmse_max = 15#np.amax( [np.amax(IR_metric['xb_rmse']),np.amax(IR_metric['xa_rmse']) ])
    bias_min = -10#np.amin( [np.amin(IR_metric['xb_bias']),np.amin(IR_metric['xa_bias']) ])
    bias_max = 5#np.amax( [np.amax(IR_metric['xb_bias']),np.amax(IR_metric['xa_bias']) ])
    
    # ------ Plot Figure -------------------
    dates = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]   
 
    # Set up figure
    if len(Expers) == 1:
        fig,ax = plt.subplots(1, 2, figsize=(16,8), dpi=300 )
        Color_set = ['red',] #'blue']
        Line_types = [['--','-'],] #['--','-']] 
        Labels = [['THO:Xb','THO:Xa'],] #['WSM6:Xb','WSM6:Xa']]
    elif len(Expers) == 2:
        fig,ax = plt.subplots(2, 2, figsize=(18,8), dpi=300 )
        Color_set = ['red','blue']
        Line_types = [['--','-'],['--','-']] 
        Labels = [['THO:Xb','THO:Xa'],['WSM6:Xb','WSM6:Xa']]
    else:
        raise ValueError("Current algorithm does not handle the color and label setting. Modify it as needed!")

    # Plot figure
    # each row: an experiment; first column: RMSE; second column: Bias
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")

    if len(Expers) == 1:
        iexper = 0
        for j in range(2):
            if j == 0:
                ax[j].plot_date(dates, IR_metric['xb_rmse'][iexper,:], color=Color_set[iexper], linestyle=Line_types[iexper][0])
                ax[j].plot_date(dates, IR_metric['xa_rmse'][iexper,:], color=Color_set[iexper], linestyle=Line_types[iexper][1])
                ax[j].set_ylim( rmse_min,rmse_max )
            elif j == 1:
                ax[j].plot_date(dates, IR_metric['xb_bias'][iexper,:], color=Color_set[iexper], label=Labels[iexper][0], linestyle=Line_types[iexper][0])
                ax[j].plot_date(dates, IR_metric['xa_bias'][iexper,:], color=Color_set[iexper], label=Labels[iexper][1], linestyle=Line_types[iexper][1])
                for i in range(len(dates)-1):
                    ax[j].annotate("", xytext=(dates[i], IR_metric['xb_bias'][iexper,i]), xy=(dates[i], IR_metric['xa_bias'][iexper,i]),arrowprops=dict(arrowstyle="->",color='black',lw=2,ls='--'))
                    ax[j].annotate("", xytext=(dates[i], IR_metric['xa_bias'][iexper,i]), xy=(dates[i+1], IR_metric['xb_bias'][iexper,i+1]),arrowprops=dict(arrowstyle='->',lw=2.5,ls='-'))    
                ax[j].set_ylim( bias_min,bias_max )
                ax[j].axhline(y=0.0,color='k',linestyle='-')
            else:
                raise ValueError("Current algorithm does not handle the color and label setting. Modify it as needed!")
    
            ax[j].set_xlim( start_time, end_time)
            ax[j].tick_params(axis='x', labelrotation=15,labelsize=12)
            
    else:
        for iexper in range(len(Expers)):
            for j in range(2):
                if j == 0:
                    ax[iexper,j].plot_date(dates, IR_metric['xb_rmse'][iexper,:], color=Color_set[iexper], linestyle=Line_types[iexper][0])
                    ax[iexper,j].plot_date(dates, IR_metric['xa_rmse'][iexper,:], color=Color_set[iexper], linestyle=Line_types[iexper][1])
                    ax[iexper,j].set_ylim( rmse_min,rmse_max )

                elif j == 1:
                    ax[iexper,j].plot_date(dates, IR_metric['xb_bias'][iexper,:], color=Color_set[iexper], label=Labels[iexper][0], linestyle=Line_types[iexper][0])
                    ax[iexper,j].plot_date(dates, IR_metric['xa_bias'][iexper,:], color=Color_set[iexper], label=Labels[iexper][1], linestyle=Line_types[iexper][1])
                    for i in range(len(dates)-1):
                        #ax[iexper,j].annotate("", xytext=(dates[i], IR_metric['xb_bias'][iexper,i]), xy=(dates[i], IR_metric['xa_bias'][iexper,i]),arrowprops=dict(arrowstyle="->",color='black',lw=2,ls='--'))
                        ax[iexper,j].annotate("", xytext=(dates[i], IR_metric['xa_bias'][iexper,i]), xy=(dates[i+1], IR_metric['xb_bias'][iexper,i+1]),arrowprops=dict(arrowstyle='->',lw=2.5,ls='-'))
                    ax[iexper,j].set_ylim( bias_min,bias_max )
                    ax[iexper,j].axhline(y=0.0,color='k',linestyle='-')
                else:
                    raise ValueError("Current algorithm does not handle the color and label setting. Modify it as needed!")
                ax[iexper,j].set_xlim( start_time, end_time)
                ax[iexper,j].tick_params(axis='x', labelrotation=15,labelsize=12)
    
    
    # legend and title
    if len(Expers) == 1:
        ax[1].legend(bbox_to_anchor=(1.25, 1.0),frameon=True,loc='upper right',fontsize='12')
        ax[0].set_title( 'RMSE',fontweight="bold",fontsize='15' )
        ax[1].set_title( 'Bias',fontweight="bold",fontsize='15' )
    else:
        ax[0,1].legend(bbox_to_anchor=(1.25, 1.0),frameon=True,loc='upper right',fontsize='12')
        ax[1,1].legend(bbox_to_anchor=(1.25, 1.0),frameon=True,loc='upper right',fontsize='12')
        ax[0,0].set_title( 'RMSE',fontweight="bold",fontsize='15' )
        ax[0,1].set_title( 'Bias',fontweight="bold",fontsize='15' )

    fig.suptitle( Storm+': Metrics of Tb (GOESR Ch8)', fontsize=15, fontweight='bold')

    # Save the figure
    save_des = small_dir+Storm+'/'+Expers[0]+'/Vis_analyze/Tb/IR_metric_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


def Gather_IR_metric( Storm, sensor, Expers, DAtimes, big_dir ):

    xb_rmse = [[] for i in range(len(Expers))]
    xa_rmse = [[] for i in range(len(Expers))]
    xb_bias = [[] for i in range(len(Expers))]
    xa_bias = [[] for i in range(len(Expers))]

    for DAtime in DAtimes:
        for idx in range(len(Expers)):
            Tb_file = big_dir+Storm+'/'+Expers[idx]+'/Obs_Hx/IR/'+DAtime+"/mean_obs_res_d03"+DAtime+'_'+ sensor+'.txt'
            if os.path.isfile( Tb_file ):
                d_all = IR_obs.read_allTb(Tb_file, sensor )
            else:
                raise ValueError(Tb_file+' does not exist!')
            
            xb_rmse[idx].append( RMSE(d_all['meanYb_obs'], d_all['Yo_obs'] ))
            xa_rmse[idx].append( RMSE(d_all['meanYa_obs'], d_all['Yo_obs'] ))
            xb_bias[idx].append( Bias(d_all['meanYb_obs'], d_all['Yo_obs'] ))
            xa_bias[idx].append( Bias(d_all['meanYa_obs'], d_all['Yo_obs'] )) 
    
    dict_IR_metric = {'xb_rmse': np.array(xb_rmse),'xa_rmse':np.array(xa_rmse),'xb_bias': np.array(xb_bias),'xa_bias':np.array(xa_bias) }
    return dict_IR_metric



if __name__ == '__main__':

    Storm = 'MARIA'
    #Expers = ['IR-J_DA+J_WRF+J_init-SP-intel17-THO-24hr-hroi900',]
    Expers = ['J_DA+J_WRF+J_init-SP-intel17-THO-24hr-hroi900','J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900']
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    Plot_IR = True
    Plot_MW = False

    # Time range set up
    start_time_str = '201709160000'
    end_time_str = '201709170000'
    Consecutive_times = True

    if not Consecutive_times:
        DAtimes = ['201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Plot the RMSE of IR Tb over time
    if Plot_IR:
        sensor = 'abi_gr'
        ch_list = ['8',]
        IR_metric = Gather_IR_metric( Storm, sensor, Expers, DAtimes, big_dir )
        Plot_IR_metric( small_dir, Storm, Expers, DAtimes, IR_metric )




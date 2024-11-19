
import os
import glob
import numpy as np
import netCDF4 as nc
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as patches
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from global_land_mask import globe
import math
from datetime import datetime, timedelta
import time
from scipy.interpolate import interp2d

import Util_data as UD
import Util_Vis
import Obspace_compare_IR_txt_bin as IR_obs
import Diagnostics as Diag



def RMSE(simu, obs):
    return np.sqrt( ((simu - obs) ** 2).mean() )

def Bias(simu, obs):
    return  np.sum((simu - obs),0)/np.size(obs,0)

def mean_Yo_Hx(simu, obs):
    return  np.sum((obs - simu),0)/np.size(obs,0)


def IR_metric_one_cycle( ist,imp ):

    bias = {}
    rmse = {}

    # First get metrics for the CONV experiment
    conv_name = UD.generate_one_name( ist,'CONV',imp )
    conv_file = big_dir+ist+'/'+conv_name+'/Obs_Hx/IR/'+DAtime[ist]+"/mean_obs_res_d03_"+DAtime[ist]+'_'+ sensor+'.txt'
    d_conv = IR_obs.read_Tb_obsRes(conv_file, sensor )
    bias['xb'] = Bias(d_conv['meanYb_obs'], d_conv['Yo_obs'] )
    bias['conv'] = Bias(d_conv['meanYa_obs'], d_conv['Yo_obs'] )
    rmse['xb'] = RMSE(d_conv['meanYb_obs'], d_conv['Yo_obs'] )
    rmse['conv'] = RMSE(d_conv['meanYa_obs'], d_conv['Yo_obs'] )
    # Get metrics for the IR experiment 
    ir_name = UD.generate_one_name( ist,'IR',imp )
    ir_file = big_dir+ist+'/'+ir_name+'/Obs_Hx/IR/'+DAtime[ist]+"/mean_obs_res_d03_"+DAtime[ist]+'_'+ sensor+'.txt'
    d_ir = IR_obs.read_Tb_obsRes(ir_file, sensor )
    bias['ir'] = Bias(d_ir['meanYa_obs'], d_ir['Yo_obs'] )
    rmse['ir'] = RMSE(d_ir['meanYa_obs'], d_ir['Yo_obs'] ) 

    return bias,rmse


def Plot_metric():

    # ------ Plot Figure -------------------
    fig,ax = plt.subplots(1, 2, dpi=300, figsize=(10,7)) #8) ) #figsize=(8,8)

    colors = {'WSM6':'red','THO':'blue'}
    lines = {'HARVEY':'-','IRMA':'--','JOSE':'-.','MARIA':':'}

    for ist in Storms:
        for imp in MP:
            x = list(d_bias[ist][imp].keys())  # Keys will be used for x-axis
            y = list(d_bias[ist][imp].values())
            ax[0].plot( x,y,linestyle=lines[ist],color=colors[imp],marker='s',linewidth='2' )
            x = list(d_rmse[ist][imp].keys())  # Keys will be used for x-axis
            y = list(d_rmse[ist][imp].values())
            ax[1].plot( x,y,linestyle=lines[ist],color=colors[imp],marker='s',linewidth='2' )

    # y = 0 line
    ax[0].axhline(y=0, color='gray', linestyle='-',linewidth=1,alpha=0.5)

    # labels
    ax[0].set_ylim([-6,6])
    ax[1].set_ylim([4,12])
    for i in range(2):
        ax[i].tick_params(axis='both', which='major', labelsize=13) 
        ax[i].set_xticklabels(['Xb','Xa:CONV','Xa:CONV+IR'])
    ax[0].set_ylabel('Bias (K)',fontsize=15) 
    ax[1].set_ylabel('RMSE (K)',fontsize=15)

    # legend: bias
    proxy_artists = [plt.Line2D([0], [0], color='red', linestyle=lw) for lw in list(lines.values()) ]
    ax[0].legend(proxy_artists,list(lines.keys()),fontsize='12',loc='upper center',ncol=2)
    # legend: rmse
    proxy_artists = [plt.Line2D([0], [0], color=color) for color in list(colors.values()) ]
    ax[1].legend(proxy_artists,list(colors.keys()),fontsize='12',loc='upper center',ncol=2)

    # title
    ax[0].set_title( 'IR Bias: mean{'+r'$\mathbf{\overline{H(X)}}$'+' - Obs}',fontsize=15 )
    ax[1].set_title( 'IR RMSE: mean{('+r'$\mathbf{\overline{H(X)}}$'+' - Obs)**2}',fontsize=15 )
    fig.suptitle('1st WRF-EnKF cycle',fontsize=18)

    des_name = small_dir+'/SYSTEMS/Vis_analyze/Tb/IR_bias_rmse_1stCycle.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


if __name__ == '__main__':


    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    MP = ['THO','WSM6']

    sensor = 'abi_gr'
    ch_list = ['8',]

    DAtime = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}

    # limitations
    limit = False
    # ------------------------------------------------------

    # Calculate metrics
    d_bias = {}
    d_rmse = {}
    for ist in Storms:
        d_bias[ist] = {}
        d_rmse[ist] = {}
        for imp in MP:
            bias,rmse = IR_metric_one_cycle( ist,imp )
            d_bias[ist][imp] = bias
            d_rmse[ist][imp] = rmse

    # Plot IR metric
    Plot_metric()





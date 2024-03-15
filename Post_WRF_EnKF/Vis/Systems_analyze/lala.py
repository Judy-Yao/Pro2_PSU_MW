
import os,fnmatch # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
from netCDF4 import Dataset
from wrf import getvar
import math
import scipy as sp
import scipy.ndimage
import matplotlib
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import time

import Util_data as UD
from Track import error_eachInit

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['lines.linewidth'] = 2.5#1.5
matplotlib.rcParams['lines.markersize'] = 2.5
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['font.size'] = 15#6

# Specify time range for model data
def t_range_model( Storm ):

    if Storm == 'HARVEY':
        DF_model_end  = '201708270000'
    elif Storm == 'IRMA':
        DF_model_end  = '201709080000'
    elif Storm == 'MARIA':
        DF_model_end  = '201709210000'
    elif Storm == 'JOSE':
        DF_model_end  = '201709100000'
    else:
        DF_model_end  = None

    return DF_model_end

# Specify time range for best-track data
def t_range_btk( Storm ):

    if Storm == 'HARVEY':
        Btk_start = '201708221200' # '201709161800' #'201709030600'
        Btk_end = '201708270000' # '201709210000' #'201709090000'
    elif Storm == 'IRMA':
        Btk_start = '201709030000'
        Btk_end = '201709080000' #201709080000
    elif Storm == 'MARIA':
        Btk_start = '2017091600'#'201709160000'
        Btk_end = '201709210000'
    elif Storm == 'JOSE':
        Btk_start = '201709050000'
        Btk_end = '201709100000'
    else:
        pass
    return Btk_start, Btk_end

# Customize colormap 
def csm_color( ida, imp ):

    num_colors = 14
    if ida == 'conv' and imp == 'WSM6':
        cm_one = cm.spring.reversed()
        discretize_cm = ListedColormap(cm_one(np.linspace(0,1,num_colors)))
        return discretize_cm.colors[0:]
    if ida == 'IR' and imp == 'WSM6':
        cm_one = cm.Reds
        discretize_cm = ListedColormap(cm_one(np.linspace(0,1,num_colors)))
        return discretize_cm.colors[8:]
    if ida == 'IR' and imp == 'TuneWSM6':
        cm_one = cm.spring.reversed()
        discretize_cm = ListedColormap(cm_one(np.linspace(0,1,num_colors)))
        return discretize_cm.colors[8:]
    elif ida == 'IR+MW' and imp == 'WSM6':
        cm_one = cm.spring 
        discretize_cm = ListedColormap(cm_one(np.linspace(0,1,num_colors)))
        return discretize_cm.colors
    elif ida == 'conv' and imp == 'THO':
        cm_one = cm.summer
        discretize_cm = ListedColormap(cm_one(np.linspace(0,1,num_colors)))
        return discretize_cm.colors[8:]
    elif ida == 'IR' and imp == 'THO':
        cm_one = cm.winter
        discretize_cm = ListedColormap(cm_one(np.linspace(0,1,num_colors)))
        return discretize_cm.colors[0:]
    elif ida == 'IR+MW' and imp == 'THO':
        cm_one = cm.Purples.reversed()
        discretize_cm = ListedColormap(cm_one(np.linspace(0,1,num_colors)))
        return discretize_cm.colors[0:]

def calculate_abs_error():

    Exper_initTimes = {}
    for istorm in Storms:
        Exper_initTimes[istorm] = {}
        for imp in MP:
            Exper_initTimes[istorm][imp] = {}
            for ida in DA:
                iExper = Exper_names[istorm][imp][ida]
                if iExper is not None:
                    Exper_initTimes[istorm][imp][ida] = sorted(fnmatch.filter(os.listdir( wrf_dir+'/'+istorm+'/'+iExper+'/wrf_df/' ),'20*00'))
                else:
                    Exper_initTimes[istorm][imp][ida] = None

    # Calculate the error with lead times
    Exper_error = {}
    for istorm in Storms:
        Exper_error[istorm] = {}
        # Read best-track data
        Btk_start, Btk_end = t_range_btk( istorm )
        best_track = UD.btk_in_duration(istorm, Btk_start, Btk_end, hour_step=6)
        # Time range for model data
        DF_model_end = t_range_model( istorm )
        for imp in MP:
            Exper_error[istorm][imp] = {}
            for ida in DA:
                iExper = Exper_names[istorm][imp][ida]
                print(iExper)
                if iExper is not None:
                    print('Calculating HPI error for '+istorm+' :'+Exper_names[istorm][imp][ida])
                    Exper_error[istorm][imp][ida] = error_eachInit(istorm,wrf_dir,iExper,best_track,Exper_initTimes[istorm][imp][ida],DF_model_end,fc_run_hrs)
                else:
                    Exper_error[istorm][imp][ida] = None  

    return Exper_error

def plog_sys_err_evolution( ):

    # Customize colormap

    # ------------- Plot ---------------------------
    fig, ax=plt.subplots(1,3,sharex='all', figsize=(15,6.5), dpi=400)

    #mean_color = ['#D70040','#0F52BA']
    #Line_types = ['-','-']
    # X axis: leading times
    lead_t = list(range(0, fc_run_hrs+1, 6))

    
    #for istorm in Storms:
    for imp in MP:
        for ida in DA:
            colorset = csm_color( ida, imp )
            err_mean_storm = np.zeros((3,len(lead_t)))
            err_mean_woMaria = np.zeros((3,len(lead_t)))
            for istorm in Storms:
                ic = Storms.index( istorm )
                # Plot each storm
                err = Exper_error[istorm][imp][ida]['mean_abs_error']
                ax[0].plot(lead_t,err[0,:],color=colorset[ic],linestyle='-',linewidth=3,alpha=0.5)
                ax[1].plot(lead_t,err[1,:],color=colorset[ic],linestyle='-',linewidth=3,alpha=0.5)
                ax[2].plot(lead_t,err[2,:],color=colorset[ic],linestyle='-',linewidth=3,alpha=0.5) 
                # Sum up
                err_mean_storm = err_mean_storm+err
                if istorm == 'MARIA':
                    continue
                else:
                    err_mean_woMaria = err_mean_woMaria+err
            err_mean_storm = err_mean_storm/len(Storms)
            err_mean_woMaria = err_mean_woMaria/len(Storms)-1
            ax[0].plot(lead_t,err_mean_storm[0,:],color=colorset[ic+1],linestyle='-',linewidth=7)
            ax[1].plot(lead_t,err_mean_storm[1,:],color=colorset[ic+1],linestyle='-',linewidth=7)
            ax[2].plot(lead_t,err_mean_storm[2,:],color=colorset[ic+1],linestyle='-',linewidth=7,label=imp)
            #ax[0].plot(lead_t,err_mean_woMaria[0,:],color=colorset[ic+1],linestyle='--',linewidth=7)
            #ax[1].plot(lead_t,err_mean_woMaria[1,:],color=colorset[ic+1],linestyle='--',linewidth=7)
            #ax[2].plot(lead_t,err_mean_woMaria[2,:],color=colorset[ic+1],linestyle='--',linewidth=7,label=imp+': w/o MARIA')

    # Details
    ax[2].legend(frameon=True,loc='upper left',fontsize='15')
    ax[0].set_ylim([0,200])
    ax[1].set_ylim([0,50])
    ax[2].set_ylim([0,50])
    for i in range(3):
        ax[i].set_xticks( lead_t[::2] )
        ax[i].set_xlim([0,fc_run_hrs])
        ax[i].tick_params(axis='x',labelsize=20)
        ax[i].grid(True,linewidth=1, color='gray', alpha=0.5, linestyle='-')

    # Set titles and axis label
    fig.text(0.52, 0.02, 'Forecast Lead Time (hour)', ha='center',fontsize=20)
    ax[0].set_ylabel('Track Error (km)',fontsize=15)
    ax[1].set_ylabel( 'MSLP Error (hPa)',fontsize=15)
    ax[2].set_ylabel( 'Vmax Error ($\mathregular{ms^{-1}}$)',fontsize=15)
    ax[0].set_title( 'Track (km)',fontsize = 20 )
    ax[1].set_title( 'MSLP (hPa)',fontsize = 20 )
    ax[2].set_title( 'Vmax ($\mathregular{ms^{-1}}$)',fontsize = 20 )
    fig.suptitle('Storm-averaged Mean Absolute Error',fontsize = 20)

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Model/IR_WSM6_IR_THO.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name+'!' )

 
if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY','JOSE','IRMA','MARIA']
    DA = ['IR',]
    MP = ['WSM6','THO'] #

    Plot_abs_error_wt = False
    Plot_abs_error_2d = True
    fc_run_hrs = 60 
    #------------------------------------
    wrf_dir = big_dir

    # Number of kinds of experiments
    num_kinds = len(DA)*len(MP)

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = {}
        for imp in MP:
            Exper_names[istorm][imp] = {}
            for ida in DA:
                Exper_names[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # Calculate error
    Exper_error = calculate_abs_error()

    # Plot absolute error with respect to forecast time
    if Plot_abs_error_wt:
        start_time=time.process_time()
        plog_sys_err_evolution()
        end_time = time.process_time()
        print('time needed: ', end_time-start_time, ' seconds')

    # Plot absoluate error of 2 experiments
    if Plot_abs_error_2d:
        start_time=time.process_time()
		if num_kinds != 2:
			rasie ValueError('This analysis only works for two experiments!')
		else:
        	plot_err_2d()
        end_time = time.process_time()
        print('time needed: ', end_time-start_time, ' seconds')        

    


import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from matplotlib import colors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import pickle
from numpy.ma import masked_array
import matplotlib.dates as mdates
from matplotlib.dates import date2num

import Util_data as UD

# explicit function to normalize array
def normalize_2d(matrix):
    norm = np.linalg.norm(matrix)
    matrix = matrix/norm  # normalized matrix
    return matrix

# scale the data PER DAY the abs(results) ranging from 0 to 1
# with the orignical signs kept
def scale_w_signs(array):
    # 2d array. 1st d: cycle; 2nd d: level
    min_var = np.min(np.abs(array), axis=1, keepdims=True)
    max_var = np.max(np.abs(array), axis=1, keepdims=True)
    #ptp_var = array.ptp(1) # fancy way to calculate range (max - min)
    #ptp_var = ptp_var[:,np.newaxis] # reshape the var with shape from (X,) to (X,1) 
    scaled_abs_var = (np.abs(array) - min_var)/(max_var-min_var)
    scaled_var = np.sign(array) * scaled_abs_var
    return scaled_var

def plot_twoVarincre_timeseries( ave_var1, ave_var2, geoHkm=None):

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300,constrained_layout=True )
    gs = fig.add_gridspec(3, 3, height_ratios=[0.05, 0.85, 0.10], width_ratios=[1, 0.2, 1])
    ax = fig.add_subplot(gs[1, :])

    #ax = plt.subplot(1,1,1)
    # Set up coordinates
    if not interp_P:
        xv = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]
        y_bottom = np.arange(0,5,0.5)
        y_top = np.arange(5,31,1)
        y_range = np.concatenate( (y_bottom,y_top),axis=0 )

        y_axis_rg = range(len(y_range))
        f_interp = interpolate.interp1d( y_range, y_axis_rg)
        yv = f_interp( geoHkm )
    else:
        xv = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]#range( np.shape(ave_norm_overT)[0] )
        yv = range( np.shape(ave_var1)[1])
    
    # set start and end time
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")

    # Plot
    ave_twoVars = (ave_var1 + ave_var2)/2 # mean of two increments

    # Both increments are negative
    condition1 = (ave_var1 < 0) & (ave_var2 < 0)
    twoVars_nega = np.copy(ave_twoVars) # Deep copy!!!!!
    twoVars_nega[~condition1] = np.nan
    twoV_nega = ax.imshow(twoVars_nega,cmap='Blues_r',origin='lower',aspect='auto',
                            vmin=-1,vmax=0,extent=None) #start_time, end_time, yv[0], yv[-1]])
    cax1 = fig.add_subplot(gs[2, 0])
    cba1 = fig.colorbar(twoV_nega,cax=cax1,orientation='horizontal',pad=0.1)
    cba1.set_ticks(np.linspace(-1,0,6))
    cba1.ax.tick_params(labelsize=13)
    cba1.ax.set_xlabel('Negative Agreement',fontsize=15)

    # Both increments are positive
    condition2 = (ave_var1 > 0) & (ave_var2 > 0)
    twoVars_posi = np.copy(ave_twoVars)
    twoVars_posi[~condition2] = np.nan
    twoV_posi = ax.imshow(twoVars_posi,cmap='Reds',origin='lower',aspect='auto',
                            vmin=0,vmax=1,extent=None)#start_time, end_time, yv[0], yv[-1]])
    cax2 = fig.add_subplot(gs[2, -1])
    cba2 = fig.colorbar(twoV_posi,cax=cax2,orientation='horizontal',pad=0.1)
    cba2.set_ticks(np.linspace(0,1,6))
    cba2.ax.tick_params(labelsize=12)
    cba2.ax.set_xlabel('Positive Agreement',fontsize=15)

    # Increments are in opposite direction
    condition3 = ~(condition1 | condition2)
    twoVars_disAgree = np.copy(ave_twoVars)
    twoVars_disAgree[condition3] = 1.2
    twoVars_disAgree[~condition3] = np.nan
    twoV_disAgree = ax.imshow(twoVars_disAgree,cmap=plt.cm.binary,origin='lower',aspect='auto',
                            vmin=1.1,vmax=2,extent=None)#start_time, end_time, yv[0], yv[-1]])

    # Set X/Y labels
    # set Y label
    if not interp_P:
        ylabel_like = [0.0,1.0,2.0,3.0,4.0,5.0,10.0,15.0,20.0]
        yticks = []
        list_y_range = list(y_range)
        for it in ylabel_like:
            yticks.append( list_y_range.index(it) )
        ax.set_yticks( yticks )
        ax.set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
        ax.set_ylabel('Height (KM)',fontsize=15)
        ax.set_ylim(ymin=0,ymax=25) # cut off data above 25km
    else:
        ax.set_yticks( yv[::10] )
        ax.set_yticklabels( [str(it) for it in P_of_interest[::10]],fontsize=15 )
        ax.set_ylabel('Pressure (hPa)',fontsize=15)
    # set X label
    ax.xaxis.set_major_locator(plt.FixedLocator(np.arange(0,len(DAtimes),6)))
    x_label = [time_dt.strftime("%m-%d %H") for time_dt in time_interest_dt[::6]]
    ax.xaxis.set_ticklabels( x_label )
    ax.tick_params(axis='x', labelrotation=45, labelsize=10)

   # Set title
    if specify_area:
        title_name = 'Scaled-EnKF-Increment Agreement: Qvapor and rt_vo in the circle with R='+str(radius_threshold)+' KM'
    else:
        title_name = 'Scaled-EnKF-Increment Agreement: d-mean Qvapor and rt_vo in the circle with R='+str(radius_threshold)+' KM' 

    ax.set_title( title_name,fontweight="bold",fontsize='12' )
    fig.suptitle(Storm+': '+Exper_name, fontsize=10, fontweight='bold')

    # Save the figure
    if not interp_P:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/ML_incrAgree_'+var_name[0]+'_'+var_name[1]+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    else:
        #save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/Interp_increAgree_'+var_name[0]+'_'+var_name[1]+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/Interp_increAgree_Dmean'+var_name[0]+'_CircleMean_'+var_name[1]+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


    return None



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'HARVEY'
    MP = 'TuneWSM6'
    DA = 'IR'
    v_interest = [['QVAPOR','rt_vo'],]

    start_time_str = '201708221200'
    end_time_str = '201708241200'
    Consecutive_times = True

    interp_P = True
    P_of_interest = list(range( 900,10,-20 ))
    interp_H = False
    specify_area = False
    radius_threshold = 300 #km

    If_plot_series = True
    # -------------------------------------------------------    
    Exper_name = UD.generate_one_name( Storm,DA,MP )

    # Identify DA times in the period of interest
    if not Consecutive_times:
        DAtimes = ['201709140000',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Plot the time evolution of domain-averaged increments
    if If_plot_series:
        for var_name in v_interest:
            if interp_P:
                # Read var 1: qvapor
                save1_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/DomainMean/Interp_increment_'+var_name[0]+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
                with open(save1_des,'rb') as file:
                    meta_and_data1 = pickle.load( file )
                ave_var1 = meta_and_data1['ave_var_overT']
                # Read var 2
                save2_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/CircleMean/Interp_increment_'+var_name[1]+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
                with open(save2_des,'rb') as file:
                    meta_and_data2 = pickle.load( file )
                ave_var2 = meta_and_data2['ave_var_overT']
                # scale the data per cycle with the result ranging from -1 to 1
                ave_var1_scale = scale_w_signs( ave_var1 )
                ave_var2_scale = scale_w_signs( ave_var2 )   
                # Plot
                plot_twoVarincre_timeseries( np.transpose(ave_var1_scale), np.transpose(ave_var2_scale) )
            else:
                # Read var 1
                save1_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/ML_increment_'+var_name[0]+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
                with open(save1_des,'rb') as file:
                    meta_and_data1 = pickle.load( file )
                ave_var1 = meta_and_data1['ave_var_overT']
                # Read var 2
                save2_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/ML_increment_'+var_name[1]+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
                with open(save2_des,'rb') as file:
                    meta_and_data2 = pickle.load( file )
                ave_var2 = meta_and_data['ave_var_overT']
                geoHkm_half_eta = meta_and_data2['geoHkm_half_eta']
                # normalize
                ave_var1_scale = scale_w_signs( ave_var1 )
                ave_var2_scale = scale_w_signs( ave_var2 )
                # Plot
                plot_twoVarincre_timeseries( np.transpose(ave_var1_scale), np.transpose(ave_var2_scale), geoHkm_half_eta )         











    

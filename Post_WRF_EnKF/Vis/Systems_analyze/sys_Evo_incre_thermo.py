#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
from scipy import interpolate
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from matplotlib import colors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import pickle

import Read_Obspace_IR as ROIR
from Util_Vis import HydroIncre
import Diagnostics as Diag
import Util_data as UD

# Generate time series
def generate_times( Storms, start_time_str, end_time_str ):

    dict_times = {}
    for istorm in Storms:
        time_diff = datetime.strptime(end_time_str[istorm],"%Y%m%d%H%M") - datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        dict_times[istorm] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times


# ------------------------------------------------------------------------------------------------------
#           Operation: Plot the domain-averaged value of each variable in a time series
# ------------------------------------------------------------------------------------------------------

# Determine the colormap design
def transform_Q( ave_var_overT ):

    # trasform values in ave_var_overT to a linearly increased function
    # -------------------------------------------------------- 
    # (before transformation) -1e-3,-1e-4,-1e-5,0,1e-5,1e-4,1e-3
    # Step 1: 3,4,5,0,-5,-4,-3
    # Step 2: -17,-16,-15,0,15,16,17
    ave_trans_Q = np.zeros( np.shape(ave_var_overT) )
    # Find the extreme values for two ends of the increment (in log)
    nega_incre = []
    posi_incre = []
    for ir in range( np.shape(ave_var_overT)[0] ): #Note how iteration over a 2D array works like this!
        for ic in range( np.shape(ave_var_overT)[1] ):
            #print('original value: ', ave_var_overT[ir,ic])
            if ave_var_overT[ir,ic] < 0:
                nega_tmp = 0-math.log10(abs( ave_var_overT[ir,ic] ))
                nega_incre.append( nega_tmp )
                ave_trans_Q[ir,ic] = nega_tmp - 20
               #print('transformed value: ', ave_trans_Q[ir,ic])
            elif ave_var_overT[ir,ic] > 0:
                posi_tmp = math.log10( ave_var_overT[ir,ic] )
                posi_incre.append( math.log10( ave_var_overT[ir,ic] ))
                ave_trans_Q[ir,ic] = posi_tmp + 20
                #print('transformed value: ', ave_trans_Q[ir,ic])
            else:
                ave_trans_Q[ir,ic] = 0

    nega_min = 3 #int(np.floor( np.min( nega_incre ) )) #eg., -1e-6.41183 > -1e-6 3
    nega_max = 8 #int(np.ceil( np.max( nega_incre ) )) #eg., -1e-11.98404 > -1e-12
    posi_min = -8 #int(np.floor( np.min( posi_incre ) )) #eg., 1e-11.78655 > 1e-12
    posi_max = -3 #int(np.ceil( np.max( posi_incre ) )) #eg., 1e-7.4489 > 1e-7 -3
    print('Min/Max of exponent in the negative range: ', nega_min,nega_max )
    print('Min/Max of exponent in the positive range: ', posi_min,posi_max )

    # Calculate the number of splitted colorbar
    # from center to the positive end
    n_posi = posi_max-posi_min+1#int(posi_max-posi_min+1) # one bar is white between 0 and the smallest positive value
    n_nega = nega_max-nega_min+1#int(nega_max-nega_min+1) # one bar is white between 0 and the largest negative value
    print( 'Number of positive bars:',n_posi )
    print( 'Number of negative bars:',n_nega )
    return [n_nega,nega_min,nega_max,n_posi,posi_min,posi_max,ave_trans_Q]


# Plot increment itself
def plot_var_incre_timeseries( ave_var_overT,ave_T_profile,N_times,geoHkm=None ):

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)
    # Set up coordinates
    if not interp_P:
        xv = range( np.shape(ave_var_overT)[0] ) 
        #y_bottom = np.arange(0,10,1)
        #y_middle = np.arange(10,15,0.5)
        #y_top = np.arange(15,31,1)
        #y_range = np.concatenate( (y_bottom,y_middle,y_top),axis=0 ) # control the scale
        y_bottom = np.arange(0,5,0.5)
        y_top = np.arange(5,31,1)
        y_range = np.concatenate( (y_bottom,y_top),axis=0 )

        y_axis_rg = range(len(y_range))
        f_interp = interpolate.interp1d( y_range, y_axis_rg)
        yv = f_interp( geoHkm )
        xcoor, ycoor = np.meshgrid( xv, yv )
    else:
        xv = range( np.shape(ave_var_overT)[0] )
        yv = range( np.shape(ave_var_overT)[1])
        xcoor, ycoor = np.meshgrid( xv, yv )

    if 'Q' in var_name:
        transformed = transform_Q( ave_var_overT )
        # number of colorbars
        n_nega = transformed[0]
        n_posi = transformed[3]
        Hydromap = HydroIncre( n_nega,n_posi )
        # extreme values of colorbar
        nega_min = transformed[1]
        nega_min_linear = nega_min - 20
        nega_max = transformed[2]
        nega_max_linear = nega_max - 20
        posi_min = transformed[4]
        posi_min_linear = posi_min + 20
        posi_max = transformed[5]
        posi_max_linear = posi_max + 20

        # transform again to make ave_trans_Q "linearly increased"
        # From -17,-16,-15,0,15,16,17 to -3,-2,-1,0,1,2,3
        ave_trans_Q = transformed[6]
        #print( ave_trans_Q )
        double_trans_Q =  np.zeros( np.shape(ave_trans_Q) )
        for ir in range( np.shape(ave_trans_Q)[0] ):
            for ic in range( np.shape(ave_trans_Q)[1] ):
                if ave_trans_Q[ir,ic] > 0:
                    double_trans_Q[ir,ic] = (ave_trans_Q[ir,ic]-(posi_min_linear-1))
                elif ave_trans_Q[ir,ic] < 0:
                    double_trans_Q[ir,ic] = (ave_trans_Q[ir,ic]-(nega_max_linear+1))
                else:
                    double_trans_Q[ir,ic] = 0
        print(np.amin(double_trans_Q))
        print(np.amax(double_trans_Q))
        #print( double_trans_Q )
        #bounds = list(range(nega_min_linear,nega_max_linear+1))+[0]+list(range(posi_min_linear,posi_max_linear+1))
        #bounds = list(range( int(np.floor( np.amin( double_trans_Q ) )), int(np.ceil( np.amax( double_trans_Q ) )+1)))
        bounds = [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]#list(range( int(np.floor( np.amin( double_trans_Q ) )), int(np.ceil( np.amax( double_trans_Q ) )+1)))
        incre_contourf = ax.contourf( xcoor, ycoor, np.transpose( double_trans_Q ), cmap=Hydromap, levels=bounds,extend='both')
        # Add color bar below the plot
        color_bar = fig.colorbar(incre_contourf,orientation = 'horizontal',pad=0.15,ticks=bounds)
        nega_side = list(range(nega_min,nega_max+1))
        nega_label = ['-10^'+str(it) for it in nega_side]
        posi_side = list(range(posi_min, posi_max+1))
        posi_label = ['10^'+str(it) for it in posi_side]
        color_bar.ax.set_xticklabels( nega_label+['0']+posi_label, rotation=45)
        color_bar.ax.tick_params(labelsize=12)
        color_bar.ax.set_xlabel('Increment: Xa-Xb',fontsize=15)
    else:
        if_all_posi = all( it[0] > 0 for it in ave_var_overT.tolist() )
        if_all_nega = all( it[0] < 0 for it in ave_var_overT.tolist() )
        if if_all_posi or if_all_nega is True:
            cmap = 'jet'
        else:
            cmap = 'bwr'
        bounds = np.linspace( int(np.floor(ave_var_overT.min())),int(np.ceil(ave_var_overT.max())),10 )
        bounds_format = [ "{0:.3f}".format( item ) for item in bounds]
        incre_contourf = ax.contourf( xcoor, ycoor, np.transpose( ave_var_overT ), cmap=cmap, vmin=ave_var_overT.min(), vmax=ave_var_overT.max(),levels=bounds_format,extend='both')
        # Add color bar below the plot
        color_bar = fig.colorbar(incre_contourf,orientation = 'horizontal',pad=0.15,ticks=bounds)
        bounds_str =  [ str(item) for item in bounds_format ]
        color_bar.ax.set_xticklabels( bounds_str, rotation=45)
        color_bar.ax.tick_params(labelsize=12)
        color_bar.ax.set_xlabel('Domain-mean Increment',fontsize=15)

    # Plot T profile
    if not interp_P:
        pass
    else:
        pass
        #T_contour = ax.contour( xcoor[0:-5,:], ycoor[0:-5,:], np.transpose(ave_T_profile[:,0:-5]-273.15),colors='k')
        #ax.clabel(T_contour, inline=True)

    # Set X/Y labels
    # set X label
    ax.set_xticks( xv[::4] )
    ax.set_xlabel('EnKF Cycle',fontsize=15)
    ax.tick_params(axis='x', labelsize=15)
    # set Y label
    if not interp_P:
        ylabel_like = [0.0,1.0,2.0,3.0,4.0,5.0,10.0,15.0,20.0]
        yticks = []
        list_y_range = list(y_range)
        for it in ylabel_like:
            yticks.append( list_y_range.index(it) )
        ax.set_yticks( yticks )
        ax.set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
        ax.set_ylabel('Height (km)',fontsize=15)
        ax.set_ylim(ymin=0,ymax=25) # cut off data above 25km
    else:
        ax.set_yticks( yv[::10] )
        ax.set_yticklabels( [str(it) for it in P_of_interest[::10]],fontsize=15 )
        ax.set_ylabel('Pressure (hPa)',fontsize=15)

    # Set title
    if 'Q' in var_name:
        title_name = 'EnKF Increment: '+var_name+' (kg/kg)'
    elif var_name == 'T':
        title_name = 'EnKF Increment: '+var_name+' (k)'
    else:
        title_name = 'EnKF Increment: '+var_name+' (m/s)'
    ax.set_title( title_name,fontweight="bold",fontsize='12' )
    fig.suptitle('Storms: '+DA[0]+'_'+MP[0], fontsize=15, fontweight='bold')

    # Save the figure
    if not interp_P:
        save_des = small_dir+'/SYSTEMS/Vis_analyze/EnKF/'+DA[0]+'_'+MP[0]+'_ML_increment_'+var_name+'_'+str(N_times)+'cycles.png'
    else:
        save_des = small_dir+'/SYSTEMS/Vis_analyze/EnKF/'+DA[0]+'_'+MP[0]+'_Interp_increment_'+var_name+'_'+str(N_times)+'cycles.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()
    return None

# ------------------------------------------------------------------------------------------------------
#           Operation: Read the domain-averaged increment of individual storm and average over storms
# ------------------------------------------------------------------------------------------------------
def eachVar_timeSeries( Exper_names,dict_times,var_name ): 

    if 'Q' in var_name:
        nLevel = 42

    # Check the length of DAtimes is the same for all storms
    for i in range(len(Storms)-1):
        if len(dict_times[Storms[i]]) != len(dict_times[Storms[i+1]]):
            raise ValueError('The length of times is not equal between experiments')

    N_times = len(dict_times[Storms[0]])

    # Array for the mean over all storms
    if interp_P: 
        ave_var = np.zeros( [N_times,len(P_of_interest)] )
        ave_T = np.zeros( [N_times,len(P_of_interest)] )
    else:
        ave_var = np.zeros( [N_times,nLevel] )
        ave_T = np.zeros( [N_times,nLevel] )
        geoHkm = np.zeros( [N_times,nLevel] )        

    # Loop thru storms
    for imp in MP:
        for ida in DA:
            for ist in Storms:
                print('Loading data for '+ ist)
                if interp_P:
                    if not specify_area:
                        save_des = small_dir+ist+'/'+Exper_names[ist][imp][ida]+'/Data_analyze/EnKF/DomainMean/Interp_increment_'+var_name+'_'+start_time_str[ist]+'_'+end_time_str[ist]+'.pickle'
                    else:
                        save_des = small_dir+ist+'/'+Exper_names[ist][imp][ida]+'/Data_analyze/EnKF/CircleMean/Interp_increment_'+var_name+'_'+start_time_str[ist]+'_'+end_time_str[ist]+'.pickle'
                    print(save_des)
                    # Read data from a pickle file
                    with open(save_des,'rb') as file:
                        load_data = pickle.load( file )
                    ave_var = ave_var + load_data['ave_var_overT'] 
                    ave_T = ave_T + load_data['ave_T_profile']
                    # Make average over storms
                    ave_var = ave_var/len(Storms)
                    ave_T = ave_T/len(Storms)
                else:
                    save_des = small_dir+ist+'/'+Exper_names[ist][imp][ida]+'/Data_analyze/EnKF/ML_increment_'+var_name+'_'+start_time_str[ist]+'_'+end_time_str[ist]+'.pickle'
                    # Read data from a pickle file
                    with open(save_des,'r') as file:
                        load_data = pickle.load( file )
                    ave_var = ave_var + load_data['ave_var_overT'] 
                    ave_T = ave_T + load_data['ave_T_profile'] 
                    geoHkm = geoHkm + load_data['geoHkm_half_eta'] 
                    # Make average over storms
                    ave_var = ave_var/len(Storms)
                    ave_T = ave_T/len(Storms)
                    geoHkm = geoHkm/len(Storms)

    # Plot the increament in a time series
    if If_plot_series:
        if interp_P:
            plot_var_incre_timeseries( ave_var,ave_T,N_times )
        else:
            plot_var_incre_timeseries( ave_var,ave_T,N_times,geoHkm )

    return None





if __name__ == '__main__':

    big_dir = '/scratch_S2/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storms = ['IRMA','JOSE','MARIA','HARVEY']
    MP = ['TuneWSM6',]
    DA = ['IR',]
    v_interest = ['QVAPOR',]#'QICE','QCLOUD','QRAIN','QSNOW','QGRAUP']
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708241200','IRMA':'201709050000','JOSE':'201709070000','MARIA':'201709180000'}
    Consecutive_times = True

    interp_P = True
    P_of_interest = list(range( 900,10,-20 ))
    interp_H = False
    specify_area = False

    If_plot_series = True
    # -------------------------------------------------------    

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

    # Plot time series of domain-average increment for each variable
    if If_plot_series:
        start_time=time.process_time()
        for var_name in v_interest:
            print('Calculate '+var_name+'...')
            eachVar_timeSeries( Exper_names,dict_times,var_name )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')












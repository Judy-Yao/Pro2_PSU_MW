#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os,sys,stat # functions for interacting with the operating system
import subprocess
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
from Util_Vis import HydroIncre

# Find the increment (posterior minus prior) in nc format and write it into another nc file
def Find_ncdiff( wrf_dir ):


    # Write a bash script 
    file_dir = wrf_dir+'find_increment.sh'
    with open ( file_dir, 'w') as rsh:
        rsh.write('''\
#! /bin/bash

# ncview is available to intel17 set
module restore intel17 

work_dir=$1
cd $work_dir
echo 'Enter ' $work_dir

# write the increment (posterior minus prior) in nc format
ncdiff wrf_enkf_output_d03_mean wrf_enkf_input_d03_mean wrf_d03_mean_increment 
''')

    # Change the permission 
    os.chmod( file_dir, stat.S_IRWXU)

    # Execute the bash script
    status = subprocess.call([wrf_dir+'find_increment.sh', wrf_dir])
    if status == 0:
        print("Succeed to find the increment!")
    else:
        raise ValueError('Failed to find the increment!')





# ------------------------------------------------------------------------------------------------------
#           Operation: Plot the domain-averaged value of each variable in a time series
# ------------------------------------------------------------------------------------------------------
# Plot normalized increment
def plot_var_incre_norm_timeseries( small_dir, Storm, Exper_name, DAtimes, var, baseP, ave_norm_overT ):

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)

    # Create a coordinate matrix from coordinate vectors
    xv = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]#range( np.shape(ave_norm_overT)[0] )
    yv = range( np.shape(ave_norm_overT)[1])
    xcoor, ycoor = np.meshgrid( xv, yv )
    bounds = np.linspace(-0.5, 0.5, 11)
    incre_contourf = ax.contourf( xcoor, ycoor, np.transpose( ave_norm_overT), cmap='RdBu_r', vmin=np.min(ave_norm_overT), vmax=np.max(ave_norm_overT),levels=bounds,extend='both')
    # Add color bar below the plot
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("bottom", size="5%", pad=0.8) #divider.new_vertical(size='5%', pad=0.8, pack_start = True)
    color_bar = fig.colorbar(incre_contourf,orientation = 'horizontal',pad=0.15)
    color_bar.ax.tick_params(labelsize=15)
    color_bar.ax.set_xlabel('Normalized Increment',fontsize=15)

    # Set X/Y labels
    # set X label
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M") 
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M") 
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45)
    # set Y label
    ax.set_yticks( yv[::5] )
    ax.set_yticklabels( [str(it) for it in baseP[::5]],fontsize=15 )   
    ax.set_ylabel('Pressure (hPa)',fontsize=15)
    ax1 = ax.twinx()
    ax1.set_yticks( yv[::5] );
    ax1.set_yticklabels( [str(it) for it in yv[::5]],fontsize=15 ) 
    ax1.set_ylabel('Model Level',fontsize=15)
    # Set title
    ax.set_title( 'Domain-averaged '+var,fontweight="bold",fontsize='15' )

    # Save the figure
    save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/norm_increment_'+var+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()

# Determine the colormap design
def transform_Q( ave_var_overT ):
   
    # trasform values in ave_var_overT to a linearly increased function
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
   
    nega_min = int(np.floor( np.min( nega_incre ) )) #eg., -1e-6.41183 > -1e-6
    nega_max = int(np.ceil( np.max( nega_incre ) )) #eg., -1e-11.98404 > -1e-12
    posi_min = int(np.floor( np.min( posi_incre ) )) #eg., 1e-11.78655 > 1e-12
    posi_max = int(np.ceil( np.max( posi_incre ) )) #eg., 1e-7.4489 > 1e-7
    print('Min/Max of exponent in the negative range: ', nega_min,nega_max )
    print('Min/Max of exponent in the positive range: ', posi_min,posi_max )

    # Calculate the number of splitted colorbar
    # from center to the positive end
    n_posi = int(posi_max-posi_min+1) # one bar is white between 0 and the smallest positive value
    n_nega = int(nega_max-nega_min+1) # one bar is white between 0 and the largest negative value
    print( 'Number of positive bars:',n_posi )
    print( 'Number of negative bars:',n_nega )
    return [n_nega,nega_min,nega_max,n_posi,posi_min,posi_max,ave_trans_Q]

# Plot increment itself
def plot_var_incre_timeseries( small_dir, Storm, Exper_name, DAtimes, var, P_interest, ave_var_overT ):

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)

    # Create a coordinate matrix from coordinate vectors
    xv = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]#range( np.shape(ave_norm_overT)[0] )
    yv = range( np.shape(ave_var_overT)[1])
    xcoor, ycoor = np.meshgrid( xv, yv )
    
    if 'Q' in var: 
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

        ave_trans_Q = transformed[6]
        bounds = list(range(nega_min_linear,nega_max_linear+1))+[0]+list(range(posi_min_linear,posi_max_linear+1))
        incre_contourf = ax.contourf( xcoor, ycoor, np.transpose( ave_trans_Q ), cmap=Hydromap, vmin=nega_min_linear, vmax=posi_max_linear,levels=bounds,extend='both')
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
        color_bar.ax.set_xlabel('Increment: Xa-Xb',fontsize=15)

    # Set X/Y labels
    # set X label
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45)
    # set Y label
    ax.set_yticks( yv[::10] )
    ax.set_yticklabels( [str(it) for it in P_interest[::10]],fontsize=15 )
    ax.set_ylabel('Pressure (hPa)',fontsize=15)
    #ax1 = ax.twinx()
    #ax1.set_yticks( yv[::5] );
    #ax1.set_yticklabels( [str(it) for it in yv[::5]],fontsize=15 )
    #ax1.set_ylabel('Model Level',fontsize=15)
    
    # Set title
    if 'Q' in var:
        title_name = Storm+'('+Exper_name+')'+': domain-averaged '+var+' (KG/KG)'
    elif var == 'T':
        title_name = Storm+'('+Exper_name+')'+': domain-averaged '+var+' (K)'
    else:
        title_name = Storm+'('+Exper_name+')'+': domain-averaged '+var+' (M/S)'
    ax.set_title( title_name,fontweight="bold",fontsize='15' )

    # Save the figure
    save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/increment_'+var+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    #plt.show()
    plt.close()



def eachVar_plot( big_dir, small_dir, Storm, Exper_name, DAtimes, v_interest ):

    for var_name in v_interest:
        print('Checking '+var_name+'...')
        
        # !!! Manullay set !!!
        num_levels = 42

        # Model level to pressure ( for figure display )
        ## Notice that PB doesn't vary much before and after assimilaiton
        P_hpa_overT = np.zeros( [len(DAtimes),num_levels] )
        for t_idx in range( len(DAtimes) ):
            print( 'At ' + DAtimes[t_idx] )
            ncdir = nc.Dataset(  big_dir+Storm+'/'+Exper_name+'/fc/'+DAtimes[t_idx]+'/wrf_enkf_output_d03_mean', 'r' )
            PB = ncdir.variables['PB'][0,:,:,:]
            P = ncdir.variables['P'][0,:,:,:]
            P_hpa = (PB + P)/100
            P_hpa_overT[t_idx,:] = np.mean( P_hpa.reshape( P_hpa.shape[0],-1),axis=1) # hPa

        # Construct a new array (using interpolation)
        P_of_interest = list(range( 995,49,-20 ))
        ave_var_overT = np.zeros( [len(DAtimes),len(P_of_interest)] )
        
        for DAtime in DAtimes:
            # Read variables of interest
            wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_d03_mean_increment'
            ncdir = nc.Dataset(wrf_file, 'r')

            if var_name == 'UV':
                u = ncdir.variables['U'][0,:,:,:]
                u_mean = ( u[:,:,0:-1] + u[:,:,1:])/2
                v = ncdir.variables['V'][0,:,:,:]
                v_mean = ( v[:,0:-1,:] + v[:,1:,:])/2
                var = np.sqrt( np.power(u_mean,2) + np.power(v_mean,2) )
            elif var_name == 'W':
                w = ncdir.variables[var_name][0,:,:,:] # level (staggered), lat, lon
                var = ( w[0:-1,:,:] + w[1:,:,:])/2 # interpolate to mass of the grid
            elif var_name == 'T':
                P1000MB=100000
                R_D=287
                CP=7*R_D/2
                P = ncdir.variables['P'][0,:,:,:]
                PB = ncdir.variables['PB'][0,:,:,:]
                Pres = PB + P
                t = ncdir.variables[var_name][0,:,:,:]
                var = (t+300.0)*( (Pres/P1000MB)**(R_D/CP) )
            else:
                var = ncdir.variables[var_name][0,:,:,:] # level,lat,lon 
         
            # Average the value over the whole domain
            var_mean = np.mean( var.reshape( var.shape[0],-1),axis=1)  # "-1" means the last two dimensions are multiplied
           
            # Interpolate to P level of interest
            f_interp = interpolate.interp1d( P_hpa_overT[t_idx,:].reshape( np.shape(var_mean) ), var_mean )
            var_interp = f_interp( P_of_interest )

            # Process for mixing ratio
            idx_zero = np.where( abs(var_interp) <= 1e-8 )[0]
            var_interp[idx_zero] = 0

            # Get the time index
            t_idx = np.where([DAtime == it for it in DAtimes])[0]
            ave_var_overT[t_idx,:] = var_interp
       
        # Plot the increament in a time series
        plot_var_incre_timeseries( small_dir, Storm, Exper_name, DAtimes, var_name, P_of_interest, ave_var_overT )



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # Configuration
    Storm = 'MARIA'
    Exper_name = 'J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900'
    v_interest = [ 'UV','W','T','QVAPOR','QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP']
    start_time_str = '201709160000'
    end_time_str = '201709160900'
    Consecutive_times = True
    if_ncdiff = False
    if_plot = True


    # Identify DA times in the period of interest
    if not Consecutive_times:
        DAtimes = ['201709140000',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Loop through each time to calculate the increment 
    if if_ncdiff:
        start_time=time.process_time()
        for DAtime in DAtimes:
            print("Getting the increment (posterior - prior) at "+DAtime)
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            Find_ncdiff( wrf_dir )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')        

    if if_plot:
        start_time=time.process_time()
        eachVar_plot( big_dir, small_dir, Storm, Exper_name, DAtimes, v_interest )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds') 



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

import Read_Obspace_IR as ROIR
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
   
    nega_min = int(np.floor( np.min( nega_incre ) )) #eg., -1e-6.41183 > -1e-6
    nega_max = int(np.ceil( np.max( nega_incre ) )) #eg., -1e-11.98404 > -1e-12
    posi_min = int(np.floor( np.min( posi_incre ) )) #eg., 1e-11.78655 > 1e-12
    posi_max = int(np.ceil( np.max( posi_incre ) )) #eg., 1e-7.4489 > 1e-7
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
def plot_var_incre_timeseries( ave_var_overT,ave_T_profile,geoHkm=None ):

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)
    # Set up coordinates
    if not interp_P:
        xv = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]
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
        xv = [datetime.strptime( it,"%Y%m%d%H%M") for it in DAtimes]#range( np.shape(ave_norm_overT)[0] )
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
        bounds = list(range( int(np.floor( np.amin( double_trans_Q ) )), int(np.ceil( np.amax( double_trans_Q ) )+1)))
        #print(bounds)
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
        T_contour = ax.contour( xcoor[0:-5,:], ycoor[0:-5,:], np.transpose(ave_T_profile[:,0:-5]-273.15),colors='k')
        ax.clabel(T_contour, inline=True)

    # Set X/Y labels
    # set X label
    start_time = datetime.strptime( DAtimes[0],"%Y%m%d%H%M")
    end_time = datetime.strptime( DAtimes[-1],"%Y%m%d%H%M")
    ax.set_xlim( start_time, end_time)
    ax.tick_params(axis='x', labelrotation=45)
    # set Y label
    if not interp_P:
        ylabel_like = [0.0,1.0,2.0,3.0,4.0,5.0,10.0,15.0,20.0]
        #ylabel_like = [0.0,5.0,10.0,11.0,12.0,13.0,14.0,15.0,20.0,25.0,30.0]
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
    #ax1 = ax.twinx()
    #ax1.set_yticks( yv[::5] );
    #ax1.set_yticklabels( [str(it) for it in yv[::5]],fontsize=15 )
    #ax1.set_ylabel('Model Level',fontsize=15)
    
    # Set title
    if 'Q' in var_name:
        title_name = 'EnKF Increment: '+var_name+' (KG/KG)'
    elif var_name == 'T':
        title_name = 'EnKF Increment: '+var_name+' (K)'
    else:
        title_name = 'EnKF Increment: '+var_name+' (M/S)'
    ax.set_title( title_name,fontweight="bold",fontsize='12' )
    fig.suptitle(Storm+': '+Exper_name, fontsize=10, fontweight='bold')

    # Save the figure
    if not interp_P:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/ML_increment_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    else:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/Interp_increment_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


def eachVar_plot( ):

    if 'Q' in var_name:
        nLevel = 42 

    if interp_P:
        # ---------- Interpolate to specified pressure levels ----------
        ## Notice that PB doesn't vary much before and after assimilaiton
        P_hpa_overT = np.zeros( [len(DAtimes),nLevel] )
        for t_idx in range( len(DAtimes) ):
            print( 'At ' + DAtimes[t_idx] )
            ncdir = nc.Dataset(  big_dir+Storm+'/'+Exper_name+'/fc/'+DAtimes[t_idx]+'/wrf_enkf_input_d03_mean', 'r' )
            PB = ncdir.variables['PB'][0,:,:,:]
            P = ncdir.variables['P'][0,:,:,:]
            P_hpa = (PB + P)/100
            P_hpa_overT[t_idx,:] = np.mean( P_hpa.reshape( P_hpa.shape[0],-1),axis=1) # hPa
        # Construct a new array (using interpolation)
        ave_var_overT = np.zeros( [len(DAtimes),len(P_of_interest)] )
        ave_T_profile = np.zeros( [len(DAtimes),len(P_of_interest)] )
    else:
        for t_idx in range( len(DAtimes) ):
            print( 'At ' + DAtimes[t_idx] )
            ncdir = nc.Dataset(  big_dir+Storm+'/'+Exper_name+'/fc/'+DAtimes[t_idx]+'/wrf_enkf_input_d03_mean', 'r' )
            PHB = ncdir.variables['PHB'][0,:,:,:]
            PH = ncdir.variables['PH'][0,:,:,:]
            geoHkm = (PHB+PH)/9.8/1000 # in km
            geoHkm = geoHkm.reshape( geoHkm.shape[0],-1)
            geoHkm_Dmean = np.mean( geoHkm, axis=1 )
            geoHkm_half_eta = (geoHkm_Dmean[:-1]+geoHkm_Dmean[1:])/2
            geoHkm_half_eta = np.ma.getdata(geoHkm_half_eta) 
        # Construct a new array at model level
        ave_var_overT = np.zeros( [len(DAtimes),nLevel] )
        ave_T_profile = np.zeros( [len(DAtimes),nLevel] ) 

    # Read increment of variable of interest
    for DAtime in DAtimes:
        if 'Q' in var_name:
            # Read mixting ratios of interest
            wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_d03_mean_increment'
            ncdir = nc.Dataset(wrf_file, 'r')
            var = ncdir.variables[var_name][0,:,:,:] # level,lat,lon 
        elif var_name == 'T':     
            P1000MB=100000
            R_D=287
            CP=7*R_D/2
            xa_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
            xa_ncdir = nc.Dataset(xa_file, 'r')
            xa_Pres = xa_ncdir.variables['P'][0,:,:,:] + xa_ncdir.variables['PB'][0,:,:,:]
            xa_t = xa_ncdir.variables['T'][0,:,:,:]
            xa_T = (xa_t+300.0)*( (xa_Pres/P1000MB)**(R_D/CP) )
            xb_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_input_d03_mean'
            xb_ncdir = nc.Dataset(xb_file, 'r')
            xb_Pres = xb_ncdir.variables['P'][0,:,:,:] + xb_ncdir.variables['PB'][0,:,:,:]
            xb_t = xb_ncdir.variables['T'][0,:,:,:]
            xb_T = (xb_t+300.0)*( (xb_Pres/P1000MB)**(R_D/CP) ) 
            var = xa_T-xb_T
        else:
            raise ValueError('Invalid variable!')

        # Read T profile
        P1000MB=100000
        R_D=287
        CP=7*R_D/2
        xa_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
        xa_ncdir = nc.Dataset(xa_file, 'r')
        xa_Pres = xa_ncdir.variables['P'][0,:,:,:] + xa_ncdir.variables['PB'][0,:,:,:]
        xa_t = xa_ncdir.variables['T'][0,:,:,:]
        T = (xa_t+300.0)*( (xa_Pres/P1000MB)**(R_D/CP) ) 
        T_mean = np.mean( T.reshape( T.shape[0],-1),axis=1)      
            
        # Average the value over the whole domain
        var_mean = np.mean( var.reshape( var.shape[0],-1),axis=1)  # "-1" means the last two dimensions are multiplied
        
        if not interp_P:
            # Get the time index
            t_idx = np.where([DAtime == it for it in DAtimes])[0]
            idx_zero = np.where( abs(var_mean) <= 1e-8 )[0]
            var_mean[idx_zero] = 0
            ave_var_overT[t_idx,:] = var_mean
            ave_T_profile[t_idx,:] = T_mean
        else:   
            # Interpolate to P level of interest
            fT_interp = interpolate.interp1d( P_hpa_overT[t_idx,:].reshape( np.shape(T_mean) ), T_mean )
            T_interp = fT_interp( P_of_interest )
            f_interp = interpolate.interp1d( P_hpa_overT[t_idx,:].reshape( np.shape(var_mean) ), var_mean )
            var_interp = f_interp( P_of_interest )
            # Process for mixing ratio
            idx_zero = np.where( abs(var_interp) <= 1e-8 )[0]
            var_interp[idx_zero] = 0
            # Get the time index
            t_idx = np.where([DAtime == it for it in DAtimes])[0]
            ave_var_overT[t_idx,:] = var_interp
            ave_T_profile[t_idx,:] = T_interp

    # Plot the increament in a time series
    if not interp_P:
        plot_var_incre_timeseries( ave_var_overT,ave_T_profile,geoHkm_half_eta )
    else:
        plot_var_incre_timeseries( ave_var_overT,ave_T_profile )

# ------------------------------------------------------------------------------------------------------
#           Object: increments per snapshot 
# ------------------------------------------------------------------------------------------------------

def plot_snapshot(big_dir, small_dir, Storm, Exper_name, var_name, DAtime, P_of_interest, d_model):

    # Read WRF domain
    wrf_file = wrf_dir+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )
    
    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = ROIR.read_TCvitals(small_dir+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,2.5), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    for isub in range(3):
        ax.flat[isub].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.flat[isub].coastlines(resolution='10m', color='black',linewidth=0.5)
        cs = ax.flat[isub].scatter(d_model['lon'],d_model['lat'],1.5,c=d_model['incre'][isub,:],\
                #edgecolors='none', cmap='RdBu_r', vmin=min_corr, vmax=max_corr, transform=ccrs.PlateCarree())
                edgecolors='none', cmap='RdBu_r',transform=ccrs.PlateCarree())
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
            ax.flat[isub].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

    # Colorbar
    caxes = fig.add_axes([0.2, 0.1, 0.6, 0.02])
    cbar = fig.colorbar(cs, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)

    #subplot title
    font = {'size':8,}
    for isub in range(3):
        ax.flat[isub].set_title( str(P_of_interest[isub])+' hPa', font, fontweight='bold')

    #title for all
    fig.suptitle(Storm+': '+Exper_name+'(Qvapor Xa-Xb)', fontsize=8, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(3):
        gl = ax[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top = False
        gl.xlabels_bottom = True
        if j==0:
            gl.ylabels_left = True
            gl.ylabels_right = False
        else:
            gl.ylabels_left = False
            gl.ylabels_right = False

        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 4}
        gl.ylabel_style = {'size': 6}

    # Save the figure
    save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/EnKF/increment_'+var_name+'_'+DAtime+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


def incre_snapshot( DAtime, wrf_dir, small_dir, Storm, Exper_name, var_name ):
    
    # Dimension of the domain
    xmax = 297
    ymax = 297
    nLevel = 42
    len_idx_x = xmax*ymax

    # Read the model constant variables
    mean_xa = wrf_dir + '/wrf_enkf_output_d03_mean'
    ncdir = nc.Dataset( mean_xa, 'r')
    # pressure levels
    PB = ncdir.variables['PB'][0,:,:,:]
    P = ncdir.variables['P'][0,:,:,:]
    P_hpa = (PB + P)/100 # 0 dimension: bottom to top
    P_hpa = P_hpa.reshape(nLevel,len_idx_x)
    # lon and lat
    lon_all = ncdir.variables['XLONG'][0,:,:].flatten()
    lat_all = ncdir.variables['XLAT'][0,:,:].flatten()

    # Read the increment
    if 'Q' in var_name:
        # Read mixting ratios of interest
        wrf_file =  wrf_dir + '/wrf_d03_mean_increment'
        ncdir = nc.Dataset(wrf_file, 'r')
        var = ncdir.variables[var_name][0,:,:,:] # level,lat,lon
        var = var.reshape(nLevel,len_idx_x)
    else:
        raise ValueError('Invalid variable!')

    # Specify pressure levels of interest
    P_of_interest = [100,500,850]
    
    # Interpolate the var to the pressure levels of interest
    start_time=time.process_time()
    var_P = np.zeros( (len(P_of_interest),len_idx_x),)
    for im in range(len_idx_x):
        f_interp = interpolate.interp1d( P_hpa[:,im], var[:,im])
        var_P[:,im] = f_interp( P_of_interest )
    end_time = time.process_time()
    print ('time needed for the interpolation: ', end_time-start_time, ' seconds')
    d_model = {'lon':lon_all,'lat':lat_all,'incre': var_P}

    # Plot
    plot_snapshot(big_dir, small_dir, Storm, Exper_name, var_name, DAtime, P_of_interest, d_model)
    return None




if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'JOSE'
    Exper_name = 'J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'
    v_interest = [ 'QVAPOR',]

    start_time_str = '201709050000'
    end_time_str = '201709060000'
    Consecutive_times = True
    
    interp_P = False
    P_of_interest = list(range( 995,49,-20 ))

    If_ncdiff = False
    If_plot_snapshot = False
    If_plot = True
    # -------------------------------------------------------    

    # Identify DA times in the period of interest
    if not Consecutive_times:
        DAtimes = ['201709140000',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Loop through each time to calculate the increment 
    if If_ncdiff:
        start_time=time.process_time()
        for DAtime in DAtimes:
            print("Getting the increment (posterior - prior) at "+DAtime)
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            Find_ncdiff( wrf_dir )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')        

    # Plot the increment per snapshot
    if If_plot_snapshot:
        print('------------ Plot the increment --------------')
        for DAtime in DAtimes:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            print('At '+DAtime)
            for var_name in v_interest:
                print('Plot '+var_name+'...')
                incre_snapshot( DAtime, wrf_dir, small_dir, Storm, Exper_name, var_name )

    # Plot the time evolution of domain-averaged increments
    if If_plot:
        start_time=time.process_time()
        for var_name in v_interest:
            print('Plot '+var_name+'...')
            eachVar_plot( )
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds') 
        
        



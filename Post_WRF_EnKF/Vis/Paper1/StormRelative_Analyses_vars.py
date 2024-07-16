import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
from wrf import getvar, ll_to_xy, to_np, interplevel
import matplotlib
from scipy import interpolate
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import subprocess
import pickle

import Util_data as UD
import EnKF_minSLP_track as SC #StormCenter

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['lines.linewidth'] = 2.5#1.5
matplotlib.rcParams['lines.markersize'] = 2.5
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['font.size'] = 8


# ------------------------------------------------------------------------------------------------------
#           Operation: Calculate, Read, and Save
# ------------------------------------------------------------------------------------------------------

# Calculate mean absolute error
def mean_absolute_error(obs,model):
    return np.mean(np.abs(np.subtract(obs, model)))

# Generate time series
def generate_times( Storms, start_time_str, end_time_str, interval ):

    dict_times = {}
    for istorm in Storms:
        time_diff = datetime.strptime(end_time_str[istorm],"%Y%m%d%H%M") - datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+interval, interval))]
        dict_times[istorm] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times

# Read fields from WRF outputs
def read_var( xa_ncdir,ivar ):

    if ivar == 'U': #U-component of Wind on Mass Points
        var = getvar( xa_ncdir,'ua',units='m s-1')
        var = to_np( var ) #Convert to NumPy array
    elif ivar == 'V': #V-component of Wind on Mass Points
        var = getvar( xa_ncdir,'va',units='m s-1')
        var = to_np( var ) #Convert to NumPy array
    elif ivar == 'WIND':
        var = getvar( xa_ncdir,'wspd_wdir',units='m s-1')
        var = to_np( var ) #Convert to NumPy array
        var = var[0,:,:] # wind speed only
    elif ivar == 'W': #W-component of Wind on Mass Points
        var = getvar( xa_ncdir,'wa',units='m s-1')
        var = to_np( var ) #Convert to NumPy array
    elif ivar == 'T': #Temperature in Kelvin
        var = getvar( xa_ncdir,'temp',units='K')
        var = to_np( var ) #Convert to NumPy array
    elif ivar == 'QVAPOR':
        var = xa_ncdir.variables[ivar][0,:,:,:] # level,lat,lon
    elif ivar == 'REFL_10CM':
        var = xa_ncdir.variables[ivar][0,:,:,:] # level,lat,lon
    else:
        pass
    return var


# calculate the azimuthal mean at different radius 
# inside a circle centered at the storm center
def calc_azimuthal_mean( itp_var, center):
    
    # Creates two arrays: one for the row indices (y) and one for the column indices (x)
    row_idx, col_idx = np.indices(itp_var.shape)

    # calculate the radius at each grid point
    x = col_idx - center.values[0] #center.values[0] will contain the X (west_east) values.
    y = row_idx - center.values[1]
    r_grid = d_reso*np.sqrt(x**2 + y**2) # resolution

    # Determine bins of radius 
    r_bins = np.arange(0, max_radius, r_interval)
    var_azi_mean = np.zeros(len(r_bins) - 1, dtype=np.float64)

    # average for a certain bin
    for i in range(len(r_bins) - 1):
        mask = (r_grid >= r_bins[i]) & (r_grid < r_bins[i+1])
        # Check if the mask is not empty to avoid computing mean of an empty slice
        if np.any(mask):
            var_azi_mean[i] = np.nanmean(itp_var[mask])
        else:
            var_azi_mean[i] = np.nan

    return var_azi_mean


def itp2Pres_AzimuMean( xa_dir,r_bins,center_ij ):

    # Read pressure levels
    xa_ncdir = nc.Dataset( xa_dir, 'r')
    pres = getvar(xa_ncdir, "pressure") # hPa

    # Interpolate
    itp_vars_mean = {}

    # Loop thru vars 
    for ivar in var_names:
        # Initiate containers
        if interp_P: # ---------- Interpolate to specified pressure levels ----------
            # Construct a new array (using interpolation)
            itp_vars_mean[ivar] = np.full( (len(P_interest),len(r_bins)-1), np.nan )
        else:
            # Construct a new array at model level
            itp_vars_mean[ivar] = np.full( (nLevel,len(r_bins)-1), np.nan )
        # Read variable
        var = read_var( xa_ncdir,ivar )
        # Interpolate 
        if interp_P:
            for ip in P_interest:
                # interpolate to specified pressure level
                itp_var_xary = interplevel(var, pres, ip)
                itp_var = to_np(itp_var_xary)
                # calculate the azimuthal mean at that level
                itp_vars_mean[ivar][P_interest.index(ip),:] = calc_azimuthal_mean( itp_var, center_ij )
        else:
            pass

    return itp_vars_mean


# Azimuthal mean calculation of 3D variables 
def azimuMean_3Dvar( Storm,Exper_name,DAtimes):

    start_time=time.process_time()

    # Dimension
    xmax = 297
    ymax = 297
    nLevel = 42

    # bins along radius 
    r_bins = np.arange(0, max_radius, r_interval)

    # Find the storm center
    if model_center:
        d_model =  SC.model_minSLP( big_dir,Storm,Exper_name,DAtimes,True) # xb_slp,xa_lat,xa_lon
    else: # best track
        # Read the best track position every 6 hours
        d_btk6Hrs = UD.read_bestrack(Storm)
        # Read the hourly best track position 
        d_btkHr = UD.interpolate_locations( DAtimes, d_btk6Hrs) # time,lat,lon

    # Loop thru EnKF cycles
    for DAtime in DAtimes:

        print('At DAtime: '+DAtime)
        xa_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
        # Get the time index
        t_idx = DAtimes.index( DAtime )

        # Find the storm center in grid space
        xa_ncdir = nc.Dataset( xa_dir, 'r')
        if model_center:
            center_ij = ll_to_xy(xa_ncdir,d_model['xa_lat'][t_idx],d_model['xa_lon'][t_idx])
        else:
            center_ij = ll_to_xy(xa_ncdir,d_btkHr['lat'][t_idx],d_btkHr['lon'][t_idx]) 

        # Read and calculate the mean
        itp_vars_mean = itp2Pres_AzimuMean( xa_dir,r_bins,center_ij)

        # Save data
        # Metadata
        current_datetime = datetime.now()
        formatted_datetime = current_datetime.strftime('%Y-%m-%d %H:%M:%S')
        # create data
        for var_name in var_names:
            if interp_P:
                metadata = {'created_at':formatted_datetime, 'Interpolated_to': 'Pressure (hPa)','Interpolated_at':P_interest,'max_radius':max_radius,'r_interval':r_interval}
                if model_center:
                    save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/Az_Mean/ModelCenter_Interp_'+var_name+'_'+DAtime+'.pickle'
                else:
                    save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/Az_Mean/BTKCenter_Interp_'+var_name+'_'+DAtime+'.pickle'
                # create a dictionary with metadata and data
                meta_and_data = {'metadata':metadata,'var_AZmean':itp_vars_mean[var_name]}
            else:
                pass

            # Write the dictionary to a pickle file
            with open(save_des,'wb') as file:
                pickle.dump( meta_and_data, file )
            print( 'Saving the data: ', save_des )


    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    return None


# Load azimuthal mean of 3D fields
def Load_AZmean_3Dvar( Storm,exp_name,DAtimes,ivar):

    var_mean = {}
    for DAtime in DAtimes:
        # file name
        if interp_P:
            if model_center:
                save_des = small_dir+Storm+'/'+exp_name+'/Data_analyze/EnKF/Az_Mean/ModelCenter_Interp_'+ivar+'_'+DAtime+'.pickle'
            else:
                save_des = small_dir+Storm+'/'+exp_name+'/Data_analyze/EnKF/Az_Mean/BTKCenter_Interp_'+ivar+'_'+DAtime+'.pickle'
        else:
            pass 
        # load file
        with open(save_des,'rb') as file:
            meta_and_data = pickle.load( file )
        # store in a dict
        var_mean[DAtime] = meta_and_data['var_AZmean'] 
    
    return var_mean

# For each WRF-EnKF experiment, average the az_mean over all cycles
def Average_allCycles( Storm, DAtimes, array ):

    mean_cycles = {}
    for ivar in var_names:
        sum_arr = 0    
        for DAtime in DAtimes:
            sum_arr = sum_arr + array[ivar][DAtime]
        mean_cycles[ivar] = sum_arr/len(DAtimes)
    return mean_cycles

# For each kind of MP*DA experiment, average the az_mean over all storms
def Average_allStorms( az_cycle_mean ):

    mean_storms = {}
    for imp in MP:
        mean_storms[imp] = {}
        for ida in DA:
            mean_storms[imp][ida] = {}
            for iv in var_names:
                sum_arr = 0
                for ist in Storms:
                    sum_arr = sum_arr + az_cycle_mean[ist][imp][ida][iv]
                mean_storms[imp][ida][iv] = sum_arr/len(Storms)
    return mean_storms

# ------------------------------------------------------------------------------------------------------
#           Operation: Plot
# ------------------------------------------------------------------------------------------------------

# test plot
def test_plot( ist,imp,ida,ivar,array ):

    # Set up figure
    fig,ax = plt.subplots( figsize=(6.5,8.5),dpi=200 )

    # Set x and y coord 
    r_bins = np.arange(0, max_radius, r_interval)
    rad = (r_bins[:-1]+r_bins[1:])/2
    if interp_P:
        yv = P_interest
    # make a mesh grid
    xcoor, ycoor = np.meshgrid( rad, yv )

    # Plot
    if ist == 'All':
        bounds = np.arange(4,20+1,2)
        wind_contourf = ax.contourf( xcoor,ycoor,array[imp][ida][ivar],levels=bounds,extend='both',cmap='hot_r')
    else:
        wind_contourf = ax.contourf( xcoor, ycoor, array[ist][imp][ida][ivar], cmap='hot_r')

    # invert the y-axis
    ax.invert_yaxis()
    
    cbar = fig.colorbar(wind_contourf)

    # Save figures
    figure_des=ist+'_'+imp+'_'+ida+'_'+ivar+'.png'#plot_dir+'AZmean'+var_name+'_cycle0.png'
    plt.savefig(figure_des, dpi=400)
    print('Saving the figure: ', figure_des)

# plot the systematic analysis: azimuthal mean over all storms and all cycles
# layout:
# WSM6: conv, IR, IR+MW
# THO: conv, IR, IR+MW
def plot_storm_mean( ivar ):

    # Set up figure
    fig = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    grids = fig.add_gridspec(ncols=3,nrows=2,top=0.93,left=0.12,hspace=0.04,wspace=0.03)
    ax = {}
    for imp in MP:
        ax[imp] = {}
        for ida in DA:
            ax[imp][ida] = fig.add_subplot( grids[MP.index(imp),DA.index(ida)] )

    # Set x and y coord 
    r_bins = np.arange(0, max_radius, r_interval)
    rad = (r_bins[:-1]+r_bins[1:])/2
    if interp_P:
        yv = P_interest
    # make a mesh grid
    xcoor, ycoor = np.meshgrid( rad, yv )

    # Customization
    if ivar == 'WIND': 
        bounds = np.arange(6,20+1,2)

    # Plot 
    for imp in MP:
        for ida in DA:
            ax_ctf = ax[imp][ida].contourf( xcoor,ycoor,az_storm_mean[imp][ida][ivar],levels=bounds,extend='both',cmap='hot_r')
    
    # Create a colorbar above the first row of subplots
    cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8]) #fig.add_axes([0.925, 0.52, 0.03, 0.43])
    cbar = fig.colorbar(ax_ctf, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Wind Speed (m $\mathregular{s^{-1}}$)')
    #cbar.set_ticks([0, 5, 10, 15, 20])
    #cbar.set_ticklabels(['0%', '5%', '10%', '15%', '20%'])

    # axes attributes
    for imp in MP:
        for ida in DA:
            # y axis
            ax[imp][ida].set_ylim([900,100])
            y_ticks = [900,700,500,300,100]
            ax[imp][ida].set_yticks( y_ticks )
            if ida == 'CONV':
                ax[imp][ida].set_yticklabels([str(it) for it in y_ticks])
                ax[imp][ida].set_ylabel('Pressure (hPa)')
            else:
                ax[imp][ida].set_yticklabels([])
            # x axis
            ax[imp][ida].set_xlim([0,200])
            if ida == DA[0]:
                x_ticks = [0,50,100,150,200]
                ax[imp][ida].set_xticks( x_ticks )
                ax[imp][ida].set_xticklabels([str(it) for it in x_ticks])
            else:
                x_ticks = [50,100,150,200]
                ax[imp][ida].set_xticks( x_ticks )
                ax[imp][ida].set_xticklabels([str(it) for it in x_ticks])
            if imp == MP[1]:
                ax[imp][ida].set_xlabel('Radius (km)')
                ax[imp][ida].set_xticklabels([str(it) for it in x_ticks])
            else:
                ax[imp][ida].set_xticklabels([])
   
    # Add DA information
    for ida in DA:
        if DA.index(ida) == 0:
            fig.text(0.25,0.95,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 1:
            fig.text(0.51,0.95,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 2:
            fig.text(0.78,0.95,ida, fontsize=12, ha='center', va='center')

    # Add MP information
    fig.text(0.03,0.74,MP[0], fontsize=11, ha='center', va='center',rotation='vertical')
    fig.text(0.03,0.30,MP[1], fontsize=11, ha='center', va='center',rotation='vertical')

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/sys_AZmean_'+ivar+'.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


# plot the azimuthal mean over all cycles
# layout:
# CONV: 4 storms
# IR: 4 storms
# IR+MW: 4 storms
def plot_cycle_mean( ivar,imp ):

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200) # standard: 6.5,8.5
    grids = fig.add_gridspec(ncols=4,nrows=3,top=0.93,left=0.12,hspace=0.04,wspace=0.03)
    ax = {}
    for ida in DA:
        ax[ida] = {}
        for ist in Storms:
            ax[ida][ist] = fig.add_subplot( grids[DA.index(ida),Storms.index(ist)] )

    # Set x and y coord 
    r_bins = np.arange(0, max_radius, r_interval)
    rad = (r_bins[:-1]+r_bins[1:])/2
    if interp_P:
        yv = P_interest
    # make a mesh grid
    xcoor, ycoor = np.meshgrid( rad, yv )

    # Customization

    # Plot 
    for ida in DA:
        for ist in Storms:
            if ivar == 'WIND':
                if ist == 'IRMA':
                    bounds = np.arange(10,45+1,5)
                    irma_ctf = ax[ida][ist].contourf( xcoor,ycoor,az_cycle_mean[ist][imp][ida][ivar],levels=bounds,extend='both',cmap='magma_r')
                elif ist == 'HARVEY': 
                    bounds = np.arange(0,10+1,1.5)
                    harvey_ctf = ax[ida][ist].contourf( xcoor,ycoor,az_cycle_mean[ist][imp][ida][ivar],levels=bounds,extend='max',cmap='cividis_r')
                else:
                    bounds = np.arange(4,18+1,2)
                    others_ctf = ax[ida][ist].contourf( xcoor,ycoor,az_cycle_mean[ist][imp][ida][ivar],levels=bounds,extend='both',cmap='hot_r')
            # invert the y-axis
            ax[ida][ist].invert_yaxis()

    # Colorbar
    cbar_ax = fig.add_axes([0.92, 0.12, 0.02, 0.8]) #fig.add_axes([0.925, 0.52, 0.03, 0.43])
    cbar_irma = fig.colorbar(irma_ctf, cax=cbar_ax, orientation='vertical')
    cbar_irma.set_label('Wind Speed (m $\mathregular{s^{-1}}$)')

    cbar_ax = fig.add_axes([0.10, 0.05, 0.4, 0.02])
    cbar_harvey = fig.colorbar(harvey_ctf, cax=cbar_ax, orientation='horizontal')
    cbar_harvey.set_label('Wind Speed (m $\mathregular{s^{-1}}$)')

    cbar_ax = fig.add_axes([0.52, 0.05, 0.40, 0.02])
    cbar_others = fig.colorbar(others_ctf, cax=cbar_ax, orientation='horizontal')
    cbar_others.set_label('Wind Speed (m $\mathregular{s^{-1}}$)')

    # axes attributes
    for ida in DA:
        for ist in Storms:
            # y axis
            ax[ida][ist].set_ylim([900,100])
            y_ticks = [900,700,500,300,100]
            ax[ida][ist].set_yticks( y_ticks )
            if ist == 'HARVEY':
                ax[ida][ist].set_yticklabels([str(it) for it in y_ticks])
                ax[ida][ist].set_ylabel('Pressure (hPa)')
            else:
                ax[ida][ist].set_yticklabels([])
            # x axis
            ax[ida][ist].set_xlim([0,200])
            if ist == DA[0]:
                x_ticks = [0,50,100,150,200]
                ax[ida][ist].set_xticks( x_ticks )
                ax[ida][ist].set_xticklabels([str(it) for it in x_ticks])
            else:
                x_ticks = [50,100,150,200]
                ax[ida][ist].set_xticks( x_ticks )
                ax[ida][ist].set_xticklabels([str(it) for it in x_ticks])
            if ida == DA[-1]:
                ax[ida][ist].set_xlabel('Radius (km)')
                ax[ida][ist].set_xticklabels([str(it) for it in x_ticks])
            else:
                ax[ida][ist].set_xticklabels([])

    # Add storm information
    for ist in Storms:
        if Storms.index(ist) == 0:
            fig.text(0.21,0.95,ist, fontsize=12, ha='center', va='center')
        elif Storms.index(ist) == 1:
            fig.text(0.42,0.95,ist, fontsize=12, ha='center', va='center')
        elif Storms.index(ist) == 2:
            fig.text(0.61,0.95,ist, fontsize=12, ha='center', va='center')
        else:
            fig.text(0.8,0.95,ist, fontsize=12, ha='center', va='center')

    # Add DA information
    fig.text(0.03,0.79,DA[0], fontsize=11, ha='center', va='center',rotation='vertical')
    fig.text(0.03,0.52,DA[1], fontsize=11, ha='center', va='center',rotation='vertical')
    fig.text(0.03,0.25,DA[2], fontsize=11, ha='center', va='center',rotation='vertical')


    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/Storms_AZmean_'+imp+'_'+ivar+'.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )





if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY','JOSE','MARIA','IRMA']
    MP = ['WSM6','THO']
    DA = ['CONV','IR','IR+MW']
    # variables of interest
    var_names= ['WIND',]
    # time period
    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708231200','IRMA':'201709040000','JOSE':'201709060000','MARIA':'201709170000'}
    cycles = 25
    lead_t = list(range(0, cycles, 1))

    # vertical interpolation
    interp_P = True
    P_interest = list(range( 900,80,-20 )) # 900 to 10 hPa
    interp_H = False

    # circle
    model_center = True
    d_reso = 3 # km
    max_radius = 250 #km
    r_interval = 10 # km

    # Operation
    #!!!!!!!!!!!!!!!!!!!!!!
    calculate_ave = False # default
    #!!!!!!!!!!!!!!!!!!!!!!
    mean_all_cycles = True
    mean_all_storms = True

    # Plot
    #------------------------------------
    wrf_dir = big_dir

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = {}
        for imp in MP:
            Exper_names[istorm][imp] = {}
            for ida in DA:
                Exper_names[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # Identify DA times in the period of interest
    d_hrs = generate_times( Storms, start_time_str, end_time_str, 1 )

    # Read and calculate azimuthal averages
    if calculate_ave:
        for ist in Storms:
            for imp in MP:
                for ida in DA:
                    # make dir
                    save_dir = small_dir+ist+'/'+Exper_names[ist][imp][ida]+'/Data_analyze/EnKF/Az_Mean/'
                    savedir_exists = os.path.exists( save_dir )
                    if savedir_exists == False:
                        os.mkdir(save_dir)
                    # calculate
                    print('Calculating the azimuthal average: '+ist+' '+imp+' '+ida)
                    azimuMean_3Dvar( ist, Exper_names[ist][imp][ida], d_hrs[ist])

    # Load pickled data
    az_mean = {}
    for ist in Storms:
        az_mean[ist] = {}
        for imp in MP:
            az_mean[ist][imp] = {}
            for ida in DA:
                az_mean[ist][imp][ida] = {}
                for ivar in var_names:
                    az_mean[ist][imp][ida][ivar] = Load_AZmean_3Dvar( ist,Exper_names[ist][imp][ida],d_hrs[ist],ivar) 

    # average az_mean_var (levels*radius bins) for all cycles
    if mean_all_cycles:
        az_cycle_mean = {}
        for ist in Storms:
            az_cycle_mean[ist] = {}
            for imp in MP:
                az_cycle_mean[ist][imp] = {}
                for ida in DA:
                    az_cycle_mean[ist][imp][ida] = Average_allCycles( ist,d_hrs[ist],az_mean[ist][imp][ida] )
                    #for ivar in var_names:
                    #    test_plot( ist,imp,ida,ivar,az_cycle_mean )

        # plot
        for imp in MP:
            for ivar in var_names:
                plot_cycle_mean( ivar,imp )


    # average az_mean_var (levels*radius bins) for all cycles and for all storms
    if mean_all_cycles & mean_all_storms:
        az_storm_mean = Average_allStorms( az_cycle_mean )
        # plot
        for ivar in var_names:
            plot_storm_mean( ivar )


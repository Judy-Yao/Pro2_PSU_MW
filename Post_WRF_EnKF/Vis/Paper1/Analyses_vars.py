import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
from wrf import getvar, ll_to_xy, to_np
import matplotlib
from scipy import interpolate
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
#from matplotlib import colors
#from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import scipy as sp
#import scipy.ndimage
import time
import subprocess
import pickle

import Util_data as UD
#import metpy.calc as mpcalc
#from metpy.units import units




def find_element_index(reshaped_list, element):
    for i, row in enumerate(reshaped_list):
        if element in row:
            return (i, row.index(element))
    return None

# 2 nested dictionary
def find_value_range(nested_dict):
    all_values = []
    
    for outer_key in nested_dict:
        for inner_key in nested_dict[outer_key]:
            all_values.extend(nested_dict[outer_key][inner_key])
    
    if not all_values:
        return None, None
    
    min_value = np.min(all_values)
    max_value = np.max(all_values)
    
    return min_value, max_value

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
def read_var( xa_file,ivar ):

    xa_ncdir = nc.Dataset(xa_file, 'r')
    if ivar == 'U': #U-component of Wind on Mass Points
        var = getvar( xa_ncdir,'ua',units='m s-1')
        var = to_np( var ) #Convert to NumPy array
    elif ivar == 'V': #V-component of Wind on Mass Points
        var = getvar( xa_ncdir,'va',units='m s-1')
        var = to_np( var ) #Convert to NumPy array
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
    return var.reshape( var.shape[0],-1)
    

def interp2Pressure( xa_dir,idx_x ):

    # Read pressure levels
    ncdir = nc.Dataset( xa_dir, 'r')
    PB = ncdir.variables['PB'][0,:,:,:]
    P = ncdir.variables['P'][0,:,:,:]
    P_hpa = (PB + P)/100
    P_hpa = P_hpa.reshape( P_hpa.shape[0],-1)
    P_hpa = P_hpa[:,idx_x] # only select value in the specified area
    # Quality control: the P_hpa of lowest level is less than 900 mb 
    idx_bad = np.where( P_hpa[0,:] < 900 )[0]
    idx_all = range( P_hpa.shape[1] )
    idx_good = np.delete(idx_all, idx_bad)
    good_P_hpa = P_hpa[:,idx_good]
    # Interpolate 
    var_interp = {}
    var_interp_mean = {}
    for ivar in var_names:
        var_interp[ivar] = np.zeros( [len(P_of_interest),len(idx_good)] )
        var_interp_mean[ivar] = np.zeros( [len(P_of_interest),1] )
    # Loop thru vars 
    for ivar in var_names:
        # read variable
        var = read_var( xa_dir,ivar )
        # only index points that meet the requirement
        good_var = var[:,idx_good]
        # loop thru points
        for im in range( len(idx_good) ):
            f_interp = interpolate.interp1d( good_P_hpa[:,im], good_var[:,im] )
            var_interp[ivar][:,im]= f_interp( P_of_interest )
        # perform domain mean
        var_interp_mean[ivar] = np.nanmean( var_interp[ivar],axis=1 )

    return var_interp_mean


# Area-mean calculation of 3D variables 
def each3DVar_timeSeries_cal( Storm,Exper_name,DAtimes):

    start_time=time.process_time()

    # Dimension
    xmax = 297
    ymax = 297
    nLevel = 42

    ave_var_overT = {}
    for ivar in var_names:
        if interp_P: # ---------- Interpolate to specified pressure levels ----------
            # Construct a new array (using interpolation)
            ave_var_overT[ivar] = np.zeros( [len(P_of_interest),len(DAtimes)] )
        else:
            # Construct a new array at model level
            ave_var_overT[ivar] = np.zeros( [nLevel,len(DAtimes)] )

    if specify_area:
        # Read the best track position every 6 hours
        d_btk = UD.read_bestrack(Storm)
        # Read the hourly best track position 
        dict_btk = UD.interpolate_locations( DAtimes, d_btk)

    # Loop thru EnKF cycles
    for DAtime in DAtimes:

        print('At DAtime: '+DAtime)
        xa_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
        # Get the time index
        t_idx = DAtimes.index( DAtime )

        # might only use data in a specified area
        if specify_area:
            # Find the best-track position 
            btk_dt = [it_str for it_str in dict_btk['time'] ]#[datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]
            bool_match = [DAtime == it for it in btk_dt]
            if True in bool_match:
                idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
            else:
                idx_btk  = None
            # convert from ll to xy
            ncdir = nc.Dataset( xa_dir )
            tc_ij = ll_to_xy(ncdir, dict_btk['lat'][idx_btk], dict_btk['lon'][idx_btk])
            # What ll_to_xy returns is not the xy coordinate itself but the grid index starting from 0. 
            # (https://forum.mmm.ucar.edu/threads/what-does-wrf-python-function-ll_to_xy-returns.12248/)
            tc_i = tc_ij.values[0]
            tc_j = tc_ij.values[1]
            idx_x = UD.find_circle_area_model_ij( xa_dir, tc_i, tc_j, radius_threshold, 3)
        else:
            idx_x = np.arange(xmax*ymax)

        # Read and calculate the mean
        if interp_P: # Interpolate to P level of interest
            var_interp_mean = interp2Pressure( xa_dir,idx_x )
            for ivar in var_names:
                ave_var_overT[ivar][:,t_idx] = var_interp_mean[ivar]
        else:
            # Read height
            ncdir = nc.Dataset( xa_dir, 'r')
            PHB = ncdir.variables['PHB'][0,:,:,:]
            PH = ncdir.variables['PH'][0,:,:,:]
            geoHkm = (PHB+PH)/9.8/1000 # in km
            geoHkm = geoHkm.reshape( geoHkm.shape[0],-1)
            geoHkm = geoHkm[:,idx_x]
            geoHkm_Dmean = np.nanmean( gLoad_domainMean_3DVar_timeSerieseoHkm, axis=1 )
            geoHkm_half_eta = (geoHkm_Dmean[:-1]+geoHkm_Dmean[1:])/2
            geoHkm_half_eta = np.ma.getdata(geoHkm_half_eta)
            # Loop thru variables
            for ivar in var_names:
                # read variable
                var = read_var( xa_dir )
                # perform domain mean
                ave_var_overT[ivar][:,t_idx] = np.nanmean( var,axis=1 )
            
    # Save data
    # Metadata
    current_datetime = datetime.now()
    formatted_datetime = current_datetime.strftime('%Y-%m-%d %H:%M:%S')
    # create data
    for var_name in var_names:
        if interp_P:
            if specify_area:
                metadata = {'created_at':formatted_datetime, 'Interpolated_to': 'Pressure (hPa)','Interpolated_at':P_of_interest,'radius_threshold':radius_threshold}
                save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/CircleMean/Interp_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
            else:
                metadata = {'created_at':formatted_datetime, 'Interpolated_to': 'Pressure (hPa)','Interpolated_at':P_of_interest}
                save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/DomainMean/Interp_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
            # create a dictionary with metadata and data
            meta_and_data = {'metadata':metadata,'ave_var_overT':ave_var_overT[ivar]}
        else:
            if specify_area:
                metadata = {'created_at':formatted_datetime,'radius_threshold':radius_threshold}
                save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/CircleMean/ML_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
            else:
                metadata = {'created_at':formatted_datetime,}
                save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/EnKF/DomainMean/ML_'+var_name+'_'+DAtimes[0]+'_'+DAtimes[-1]+'.pickle'
            meta_and_data = {'metadata':metadata,'ave_var_overT':ave_var_overT[ivar],'geoHkm_half_eta':geoHkm_half_eta}
    
        # Write the dictionary to a pickle file
        with open(save_des,'wb') as file:
            pickle.dump( meta_and_data, file )
        print( 'Saving the data: ', save_des )

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    return None

# Load time-series domain-averaged 3Dvar
def Load_domainMean_3DVar_timeSeries( storm,exp_name,ivar ):

    # dir to file
    if specify_area:
        save_des = small_dir+storm+'/'+exp_name+'/Data_analyze/EnKF/CircleMean/Interp_'+ivar+'_'+start_time_str[storm]+'_'+end_time_str[storm]+'.pickle'
    else:
        save_des = small_dir+storm+'/'+exp_name+'/Data_analyze/EnKF/DomainMean/Interp_'+ivar+'_'+start_time_str[storm]+'_'+end_time_str[storm]+'.pickle'
    # load file
    with open(save_des,'rb') as file:
        meta_and_data = pickle.load( file )
    dmean = meta_and_data['ave_var_overT']

    return dmean

# Calculate differences: Exp - Control
def calc_exp_diff( dmean_3d, list_mp, list_da ):

    diff = {}
   
    for ist in Storms:
        diff[ist] = {}
        # loop thru MP list
        for imp in list_mp:
            diff[ist][imp] = {}
            # loop thru DA list
            for ida in list_da:
                diff[ist][imp][ida] = {}
                # loop thru var
                for ivar in var_names:
                    if exp_DAdiff and not exp_MPdiff:
                        diff[ist][imp][ida][ivar] = dmean_3d[ist][imp][ida][ivar] - dmean_3d[ist][imp][exp_control][ivar]
                    if exp_MPdiff and not exp_DAdiff:
                        diff[ist][imp][ida][ivar] = dmean_3d[ist][imp][ida][ivar] - dmean_3d[ist][exp_control][ida][ivar]
    return diff

# Stack arrays in nested dictionaries of experiment difference
# From Storm*MP*DA*Var to Var*MP*DA
def stack_Storms_diff( dicts, list_mp, list_da ):

    stack_storms = {}
    for ivar in var_names:
        stack_storms[ivar] = {}
        for imp in list_mp:
            stack_storms[ivar][imp] = {}
            for ida in list_da:
                if interp_P:
                    stack_sts = np.zeros((len(P_of_interest),1)) 
                else:
                    pass
                # stack
                for ist in Storms:
                    if os.path.exists( wrf_dir+'/'+ist+'/'+Exper_names[ist][imp][ida] ):
                        stack_sts = np.concatenate((stack_sts,dicts[ist][imp][ida][ivar]),axis=1)
                stack_storms[ivar][imp][ida] = stack_sts[:,1:] #remove the first empty column

    return stack_storms



def set_range( ivar):

    if ivar == 'U': #U-component of Wind on Mass Points
        return -3,3 #3
    elif ivar == 'V': #V-component of Wind on Mass Points
        return -3, 3 #3.5
    elif ivar == 'W': #W-component of Wind on Mass Points
        return -0.05, 0.05 #0.05
    elif ivar == 'T': #Temperature in Kelvin
        return -2,2 #-2,3
    elif ivar == 'QVAPOR':
        return  -0.002,0.001
    elif ivar == 'REFL_10CM':
        return -22,22 #-15, 18
    else:
        pass


# ------------------------------------------------------------------------------------------------------
#           Operation: Plot
# ------------------------------------------------------------------------------------------------------

# Layout:
# WSM6: U,V,W
# WSM6: T,QVAPOR,REFL_10CM
# THO: U,V,W
# THO: T,QVAPOR,REFL_10CM
def plot_sameDA_diffMP( h_count,ida ):

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=2,top=0.93,left=0.1,hspace=0.09)
    
    # reshape 1by6 list to 2by3 list
    vars_2by3 = [var_names[i:i+3] for i in range(0, len(var_names), 3)]

    ax = {}
    for imp in MP:
        ax[imp] = {}

    # gridspec inside gridspec
    wsm6_grid = outer_grids[0].subgridspec(2, 3, wspace=0.05, hspace=0.15)
    tho_grid = outer_grids[1].subgridspec(2, 3, wspace=0.05, hspace=0.15)

    for iv in var_names:
        loc = find_element_index(vars_2by3,iv)
        ax['WSM6'][iv] = fig.add_subplot( wsm6_grid[loc[0],loc[1]] )
        ax['THO'][iv] = fig.add_subplot( tho_grid[loc[0],loc[1]] )


    ## y axis: vertical coordinate
    if interp_P:
        y_axis_rg = P_of_interest
    else:
        pass

    # Plot 
    for imp in MP:
        for iv in var_names:
            min_var, max_var = set_range( iv )
            x_axis_rg = np.linspace(min_var, max_var,num_bins)
            X,Y = np.meshgrid(x_axis_rg,y_axis_rg)
            im = ax[imp][iv].pcolor(X,Y,h_count[iv][imp][ida], cmap='gist_heat_r',vmin=0,vmax=20)
            # invert the y-axis
            ax[imp][iv].invert_yaxis()

    # Create a colorbar above the first row of subplots
    cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8]) #fig.add_axes([0.925, 0.52, 0.03, 0.43])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cbar.set_ticks([0, 5, 10, 15, 20])
    cbar.set_ticklabels(['0%', '5%', '10%', '15%', '20%'])

    # axes attributes
    for imp in MP:
        for iv in var_names:
            ax[imp][iv].set_ylim([900,100])
            y_ticks = [900,700,500,300,100]
            ax[imp][iv].set_yticks( y_ticks )
            if iv == vars_2by3[0][0] or iv == vars_2by3[1][0]:
                ax[imp][iv].set_yticklabels([str(it) for it in y_ticks])
            else:
                ax[imp][iv].set_yticklabels([])


    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/sys_pdf_Xa_diff_'+ida+'_'+exp_control+'.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


# Layout:
# IR: U,V,W
# IR: T,QVAPOR,REFL_10CM
# IR+MW: U,V,W
# IR+MW: T,QVAPOR,REFL_10CM
def plot_sameMP_diffDA( h_count,imp,list_da ):

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=2,top=0.93,left=0.1,hspace=0.09)

    # reshape 1by6 list to 2by3 list
    vars_2by3 = [var_names[i:i+3] for i in range(0, len(var_names), 3)]

    ax = {}
    for ida in list_da:
        ax[ida] = {}

    # gridspec inside gridspec
    wsm6_grid = outer_grids[0].subgridspec(2, 3, wspace=0.05, hspace=0.15)
    tho_grid = outer_grids[1].subgridspec(2, 3, wspace=0.05, hspace=0.15)

    for iv in var_names:
        loc = find_element_index(vars_2by3,iv)
        ax['IR'][iv] = fig.add_subplot( wsm6_grid[loc[0],loc[1]] )
        ax['IR+MW'][iv] = fig.add_subplot( tho_grid[loc[0],loc[1]] )


    ## y axis: vertical coordinate
    if interp_P:
        y_axis_rg = P_of_interest
    else:
        pass

    # Plot 
    for ida in list_da:
        for iv in var_names:
            min_var, max_var = set_range( iv )
            x_axis_rg = np.linspace(min_var, max_var,num_bins)
            X,Y = np.meshgrid(x_axis_rg,y_axis_rg)
            im = ax[ida][iv].pcolor(X,Y,h_count[iv][imp][ida], cmap='gist_heat_r',vmin=0,vmax=20)
            # invert the y-axis
            ax[ida][iv].invert_yaxis()

    # Create a colorbar above the first row of subplots
    cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8]) #fig.add_axes([0.925, 0.52, 0.03, 0.43])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cbar.set_ticks([0, 5, 10, 15, 20])
    cbar.set_ticklabels(['0%', '5%', '10%', '15%', '20%'])

    # axes attributes
    for ida in list_da:
        for iv in var_names:
            ax[ida][iv].set_ylim([900,100])
            y_ticks = [900,700,500,300,100]
            ax[ida][iv].set_yticks( y_ticks )
            if iv == vars_2by3[0][0] or iv == vars_2by3[1][0]:
                ax[ida][iv].set_yticklabels([str(it) for it in y_ticks])
            else:
                ax[ida][iv].set_yticklabels([])


    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/sys_pdf_Xa_diff_'+exp_control+'_'+imp+'.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

# Layout:
# MP1 - MP0: U,V,W
# MP1 - MP0: T,QVAPOR,REFL_10CM
def plot_sameDA_betweenDiffMP( h_count,ida ):

    # Set up figure
    fig = plt.figure( figsize=(6.5,4.5),dpi=200) # standard: 6.5,8.5
    grids = fig.add_gridspec(ncols=3,nrows=2,top=0.93,left=0.15,hspace=0.15, wspace=0.05)

    # reshape 1by6 list to 2by3 list
    vars_2by3 = [var_names[i:i+3] for i in range(0, len(var_names), 3)]

    ax = {}
    for imp in list_mp:
        ax[imp+'_'+exp_control] = {}

    for iv in var_names:
        loc = find_element_index(vars_2by3,iv)
        ax[imp+'_'+exp_control][iv] = fig.add_subplot( grids[loc[0],loc[1]] )

    ## y axis: vertical coordinate
    if interp_P:
        y_axis_rg = P_of_interest
    else:
        pass

    # Plot
    for imp in list_mp: 
        for iv in var_names:
            min_var, max_var = set_range( iv )
            x_axis_rg = np.linspace(min_var, max_var,num_bins)
            X,Y = np.meshgrid(x_axis_rg,y_axis_rg)
            im = ax[imp+'_'+exp_control][iv].pcolor(X,Y,h_count[iv][imp][ida], cmap='gist_heat_r',vmin=0,vmax=20)
            # invert the y-axis
            ax[imp+'_'+exp_control][iv].invert_yaxis()

    # Create a colorbar above the first row of subplots
    cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8]) #fig.add_axes([0.925, 0.52, 0.03, 0.43])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cbar.set_ticks([0, 5, 10, 15, 20])
    cbar.set_ticklabels(['0%', '5%', '10%', '15%', '20%'])

    # axes attributes
    for imp in list_mp:
        for iv in var_names:
            ax[imp+'_'+exp_control][iv].set_ylim([900,100])
            y_ticks = [900,700,500,300,100]
            ax[imp+'_'+exp_control][iv].set_yticks( y_ticks )
            if iv == vars_2by3[0][0] or iv == vars_2by3[1][0]:
                ax[imp+'_'+exp_control][iv].set_yticklabels([str(it) for it in y_ticks])
            else:
                ax[imp+'_'+exp_control][iv].set_yticklabels([])


    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/sys_pdf_Xa_'+ida+'_diff_'+list_mp[0]+'_'+exp_control+'.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    MP = ['WSM6','THO']
    DA = ['CONV','IR','IR+MW']
    # variables of interest
    var_names= ['U','V','W','T','QVAPOR','REFL_10CM']
    # time period
    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708231200','IRMA':'201709040000','JOSE':'201709060000','MARIA':'201709170000'}
    cycles = 25
    lead_t = list(range(0, cycles, 1))

    # vertical interpolation
    interp_P = True
    P_of_interest = list(range( 900,10,-20 )) # 900 to 100 hPa
    interp_H = False
    # domain average
    dom_ave = True
    specify_area = False
    if specify_area:
        radius_threshold = 300 #km
    # histogram
    num_bins = 100

    # Operation
    #!!!!!!!!!!!!!!!!!!!!!!
    calculate_ave = False # default
    #!!!!!!!!!!!!!!!!!!!!!!

    exp_DAdiff = False
    exp_MPdiff = True

    if exp_DAdiff and not exp_MPdiff:
        exp_control = 'CONV'

    if exp_MPdiff and not exp_DAdiff:
        exp_control = 'THO'

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

    # Read and calculate average
    if calculate_ave:
        for ist in Storms:
            for imp in MP:
                for ida in DA:
                    # make dir
                    save_dir = small_dir+ist+'/'+Exper_names[ist][imp][ida]+'/Data_analyze/EnKF/'
                    if specify_area:
                        save_dir = save_dir + 'CircleMean/'
                    else:
                        save_dir = save_dir + 'DomainMean/'
                    savedir_exists = os.path.exists( save_dir )
                    if savedir_exists == False:
                        os.mkdir(save_dir)
                    # calculate
                    print('Calculating the average: '+ist+' '+imp+' '+ida)
                    each3DVar_timeSeries_cal( ist, Exper_names[ist][imp][ida], d_hrs[ist])
    
    # Load pickled data
    dmean_3d = {}
    for ist in Storms:
        dmean_3d[ist] = {}
        for imp in MP:
            dmean_3d[ist][imp] = {}
            for ida in DA:
                dmean_3d[ist][imp][ida] = {}
                for ivar in var_names:
                    dmean_3d[ist][imp][ida][ivar] = Load_domainMean_3DVar_timeSeries( ist, Exper_names[ist][imp][ida], ivar  ) 

    # Histogram of differences between DA experiments
    if exp_DAdiff:

        if exp_control == 'CONV':
            list_da = ['IR','IR+MW']
        elif exp_control == 'IR':
            list_da = ['IR+MW',]
        else:
            raise ValueError("Such condition does not exist!")

        # Make differences at each time: Exp - control
        diff = calc_exp_diff( dmean_3d, MP, list_da )

        # stack storms: data from storms are treated equally
        # nested by Var*MP*DA
        stack_storms = stack_Storms_diff( diff, MP, list_da )

        # Find the range for each var
        #for ivar in var_names:
        #    min_value, max_value = find_value_range(stack_storms[ivar])
        #    print("The range of "+ivar+f" is from {min_value} to {max_value}.")

        # Make histogram
        h_count = {}
        for ivar in var_names:
            print('Bin for '+ ivar)
            h_count[ivar] = {}
            # set up bins and ranges 
            min_var_rg, max_var_rg = set_range( ivar )
            # loop thru MP:
            for imp in MP:
                h_count[ivar][imp] = {}
                for ida in list_da:
                    if interp_P:
                        H_one = np.zeros((len(P_of_interest), num_bins))
                        for il in range(len(P_of_interest)):
                            temp_counts,bin_edges = np.histogram(stack_storms[ivar][imp][ida][il,:], num_bins, (min_var_rg,max_var_rg))
                            percentages = 100 * temp_counts / len(stack_storms[ivar][imp][ida][il,:])
                            H_one[il,:] = percentages
                        h_count[ivar][imp][ida] = H_one
                    else:
                        raise ValueError("The algorithm is not developed yet!")
            
        # plot
        # one DA + different MP
        #for ida in list_da:
        #    plot_sameDA_diffMP( h_count,ida )
        for imp in MP:
            plot_sameMP_diffDA( h_count,imp,list_da )            


    # Histogram of differences between MP experiments
    if exp_MPdiff:

        if exp_control == 'THO':
            list_mp = ['WSM6',]
        elif exp_control == 'WSM6':
            list_mp = ['THO',]
        else:
            raise ValueError("Such condition does not exist!")

        # Make differences at each time: Exp - control
        diff = calc_exp_diff( dmean_3d, list_mp, DA )

        # stack storms: data from storms are treated equally
        # nested by Var*MP*DA
        stack_storms = stack_Storms_diff( diff, list_mp, DA )

        # Find the range for each var
        #for ivar in var_names:
        #    min_value, max_value = find_value_range(stack_storms[ivar])
        #    print("The range of "+ivar+f" is from {min_value} to {max_value}.")

        # Make histogram
        h_count = {}
        for ivar in var_names:
            print('Bin for '+ ivar)
            h_count[ivar] = {}
            # set up bins and ranges 
            min_var_rg, max_var_rg = set_range( ivar )
            # loop thru MP:
            for imp in list_mp:
                h_count[ivar][imp] = {}
                # loop thru DA
                for ida in DA:
                    if interp_P:
                        H_one = np.zeros((len(P_of_interest), num_bins))
                        for il in range(len(P_of_interest)):
                            temp_counts,bin_edges = np.histogram(stack_storms[ivar][imp][ida][il,:], num_bins, (min_var_rg,max_var_rg))
                            percentages = 100 * temp_counts / len(stack_storms[ivar][imp][ida][il,:])
                            H_one[il,:] = percentages
                        h_count[ivar][imp][ida] = H_one
                    else:
                        raise ValueError("The algorithm is not developed yet!")

        # plot: between different MP experiment
        for ida in ['CONV',]:
            plot_sameDA_betweenDiffMP( h_count,ida )




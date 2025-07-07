
from numba import njit, prange
import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
from wrf import getvar,interplevel
import numpy.ma as ma
import scipy as sp
#from scipy import interpolate
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from matplotlib import colors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import pickle
import random

import Util_Vis
import Util_data as UD
import Diagnostics as Diag


## For each obs loc, find the nearest model grid point (good for mercator)
@njit(parallel=True)
def nearest_axis( obs,model ):

    res = np.zeros( (obs.shape[0]), )
    res[:] = np.nan
    for io in prange( obs.shape[0] ):
        for i in range( model.shape[0]-1 ):
            if model[i+1] < obs[io]:
                continue
            elif model[i] > obs[io]:
                continue
            else:
                if model[i] <= obs[io] and model[i+1] >= obs[io]:
                    if abs(model[i]-obs[io]) < abs(model[i+1]-obs[io]):
                        res[io] = i
                    else:
                        res[io] = i+1
    # Make sure every obs has a nearest model grid
    assert res.any() != np.nan
    return res

def Find_nearest_grid( wrf_file,lon_obs,lat_obs ): #lon_obs,lat_obs: lists

    print('------- Search for the nearest model grid for the obs ------')
    
    # Read model lon and lat
    ncdir = nc.Dataset( wrf_file, 'r')
    lon_x1d = ncdir.variables['XLONG'][0,0,:]
    lon_x1d = lon_x1d.filled(np.nan)
    lat_x1d = ncdir.variables['XLAT'][0,:,0]
    lat_x1d = lat_x1d.filled(np.nan)
    lon_x = ncdir.variables['XLONG'][0,:,:].flatten()
    lat_x = ncdir.variables['XLAT'][0,:,:].flatten()

    # Loop thru all obs and search each's left- and bottom-nearest model grid along x and y direction
    # returned is i,j in model domain
    Idx_i = nearest_axis( lon_obs,lon_x1d )
    Idx_j = nearest_axis( lat_obs,lat_x1d )

    # Transform to a 2d meshgrid and find the corresponding idx
    Idx_nearest = []
    for io in range(len(lon_obs)):
        Idx_nearest.append( int(Idx_j[io]*len(lon_x1d)+Idx_i[io]) )

    # check
    check_idx = random.randint(0, len(lon_obs))-1
    print('Checking the '+str(check_idx)+'th obs... lon: '+str(lon_obs[check_idx])+' lat: '+str(lat_obs[check_idx]))
    print('Its nearest model grid -- lon: '+str( lon_x[Idx_nearest[check_idx]] )+' lat: '+str( lat_x[Idx_nearest[check_idx]] ))

    if len(np.unique(Idx_nearest)) != len(lon_obs):
        warnings.warn('The nearest model grids might be repetitive!')
        print('Number of obs is '+str(len(lon_obs))+' while the number of the unique model location is '+str(len(np.unique(Idx_nearest))))

    return Idx_nearest


# ---------------------------------------------------------------------------------------------------------------
#    Object: explain the PSFC horizontal correlation 
# ------------------------------------------------------------------------------------------------------------------

def find_mslp( Storm,wrf_file):

    ncdir = nc.Dataset( wrf_file )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # sea level pressure
    #slp = getvar(ncdir, 'slp')
    slp = UD.compute_slp( ncdir )
    #print( np.amax( abs(slp - slp_getvar.values )) )
    # original SLP
    slp_values = slp #slp.values
    slp_values[slp_values > 1030] = np.nan
    # smoothed SLP
    slp_smt_values = sp.ndimage.gaussian_filter(slp, [1,1]) #[11,11]
    slp_smt_values[slp_smt_values > 1030] = np.nan
    # simulated storm center
    if Storm == 'HARVEY': # Harvey is special with its location near land!

        # eyeball where the storm is
        if DAtime <= '201708221600':
            lat_mask = lat <= 18
            lon_mask = lon <= -91.5
            mask = lat_mask | lon_mask
        elif (DAtime >= '201708221700') & (DAtime <= '201708222300'):
            lat_mask = lat <= 18
            lon_mask =  (lon >= -88) | (lon <= -91.5)
            mask = lat_mask | lon_mask
        else:
            mask = lat <= 18
        slp_masked = ma.masked_array(slp_values, mask=mask)
        minslp = np.nanmin( slp_masked ) 

        slp_smooth_masked = ma.masked_array(slp_smt_values, mask=mask)
        idx = np.nanargmin( slp_smooth_masked )
        lat_minslp = lat.flatten()[idx] 
        lon_minslp = lon.flatten()[idx] 
    else:
        minslp = np.nanmin( slp_values ) 
        idx = np.nanargmin( slp_smt_values )
        lat_minslp = lat.flatten()[idx]       
        lon_minslp = lon.flatten()[idx] 
    return {'mslp':minslp,'lat_minslp':lat_minslp,'lon_minslp':lon_minslp}

# Loop thru all ensemble members &
# Identify the mslp for that member &
# Save mslps to a txt file
def identify_mslp_ens( wrf_dir, key='input'):

    mslp_info = np.zeros( shape=(num_ens,4) )
    #mslp_info = np.nan

    for ie in range(num_ens):
        id_num = ie+1
        wrf_file = wrf_dir+'wrf_enkf_'+key+'_d03_'+f"{id_num:03}" 
        d_mslp = find_mslp( Storm,wrf_file)
        mslp_info[ie,1] = "{0:.3f}".format( d_mslp['lon_minslp'] )
        mslp_info[ie,2] = "{0:.3f}".format(d_mslp['lat_minslp'] )
        mslp_info[ie,3] = "{0:.3f}".format(d_mslp['mslp'] )

    # May save the mslp info
    if If_save:
        header = ['ID','Lon','Lat','mslp']
        des_path = wrf_dir+ DAtime+'_enkf_'+key+'_mslp.txt'
        with open(des_path,'w') as f:
            # Add header 
            f.write('\t'.join( item.rjust(6) for item in header ) + '\n' )
            # Write the record to the file serially
            len_records = np.shape( mslp_info )[0]
            for irow in range( len_records ):
                irecord = mslp_info[irow,:].astype(str) 
                irecord[0] = f"{irow+1:03}"
                #print(irecord)
                f.write('\t'.join( item.rjust(6) for item in irecord ) + '\n')
    print('Save '+des_path)
    
    return None

# Find the slp at the obs location
def find_slp_obsLoc( wrf_file, idx_grid ):

    ncdir = nc.Dataset( wrf_file )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # sea level pressure
    #slp = getvar(ncdir, 'slp')
    slp = UD.compute_slp( ncdir )
    # simulated storm center
    if Storm == 'HARVEY': # Harvey is special with its location near land!
        pass
    else:
        lon_slp = lon.flatten()[idx_grid].item()
        lat_slp = lat.flatten()[idx_grid].item()
        slp_obsLoc = slp.flatten()[idx_grid].item()
    return {'slp_obsLoc':slp_obsLoc,'lat_slp':lat_slp,'lon_slp':lon_slp}

# Loop thru all ensemble members &
# Identify the slp for that member at the obs location &
# Save mslps to a txt file
def identify_slp_ens_obsLoc( obs, wrf_dir, key=None ):

    slp_obsLoc_info = np.zeros( shape=(num_ens,4) )

    # Find the nearest grid point to the obs
    lon_obs = np.array( obs['lon'] )
    lat_obs = np.array( obs['lat'] )
    wrf_file = wrf_dir+'wrf_enkf_'+key+'_d03_001'
    idx_grid = Find_nearest_grid( wrf_file,lon_obs,lat_obs ) 

    # Identify the slp at the obs location
    for ie in range(num_ens):
        id_num = ie+1
        wrf_file = wrf_dir+'wrf_enkf_'+key+'_d03_'+f"{id_num:03}"
        d_slp_obsLoc = find_slp_obsLoc( wrf_file, idx_grid )
        slp_obsLoc_info[ie,1] = "{0:.3f}".format( d_slp_obsLoc['lon_slp'] )
        slp_obsLoc_info[ie,2] = "{0:.3f}".format(d_slp_obsLoc['lat_slp'] )
        slp_obsLoc_info[ie,3] = "{0:.3f}".format(d_slp_obsLoc['slp_obsLoc'] )

    # May save the slp_obsLoc info
    if If_save:
        header = ['ID','Lon','Lat','slp_obsLoc']
        des_path = wrf_dir+ DAtime+'_enkf_'+key+'_slp_obsLoc.txt'
        with open(des_path,'w') as f:
            # Add header 
            f.write('\t'.join( item.rjust(6) for item in header ) + '\n' )
            # Write the record to the file serially
            len_records = np.shape( slp_obsLoc_info )[0]
            for irow in range( len_records ):
                irecord = slp_obsLoc_info[irow,:].astype(str)
                irecord[0] = f"{irow+1:03}"
                #print(irecord)
                f.write('\t'.join( item.rjust(6) for item in irecord ) + '\n')
    print('Save '+des_path)

    return None


def read_mslp_ens( wrf_dir,DAtime,key=None ):

    id_ = []
    lon_mslp = []
    lat_mslp = []
    mslp = []

    with open( wrf_dir+ DAtime+'_enkf_'+key+'_mslp.txt' ) as f:
        next(f)
        all_lines = f.readlines()

    for line in all_lines:
        split_line = line.split()
        id_.append( int(split_line[0]) )
        lon_mslp.append( float(split_line[1]) )
        lat_mslp.append( float(split_line[2]) )
        mslp.append( float(split_line[3]) )
    
    d_mslp = {'id':id_,'mslp':mslp,'lat_mslp':lat_mslp,'lon_mslp':lon_mslp}
    return d_mslp

# Read  HPI info from pre-calculated files
# e.g., HPI_wrf_enkf_input_d03_.201709030000_201709040000.txt 
def Read_mslp_EnsMean( mfile ):

    d_mean = {}
    print('Reading ', mfile)
    with open(mfile) as f:
        next(f)
        all_lines = f.readlines()

    for line in all_lines:
        if 'Success' in line:
            break
        split_line = line.split()
        time = split_line[0]
        lon =  float(split_line[1])
        lat =  float(split_line[2])
        mslp = float(split_line[3])
        d_mean[time] = (lon,lat,mslp) 

    return d_mean


def read_min_PSFC_ens( wrf_dir,DAtime,key=None ):

    id_ = []
    lon_mpsfc = []
    lat_mpsfc = []
    mpsfc = []

    for ie in range(num_ens):
        id_.append( ie+1 )
        wrf_file = wrf_dir+'wrf_enkf_'+key+'_d03_'+f"{ie+1:03}"
        ncdir = nc.Dataset( wrf_file )
        lat = ncdir.variables['XLAT'][0,:,:]
        lon = ncdir.variables['XLONG'][0,:,:]
        PSFC = ncdir.variables['PSFC'][0,:,:]

        mpsfc.append( np.nanmin( PSFC )/100 )
        idx = np.nanargmin( PSFC )
        lat_mpsfc.append( lat.flatten()[idx] )
        lon_mpsfc.append( lon.flatten()[idx] )


    d_mpsfc = {'id':id_,'mslp':mpsfc,'lat_mslp':lat_mpsfc,'lon_mslp':lon_mpsfc}
    return d_mpsfc

def mean_geolocation(coordinates):
    """
    Calculates the mean latitude and longitude of a set of geographic coordinates.

    Parameters:
    coordinates (list of tuples): List of (longitude, latitude) in degrees.

    Returns:
    tuple: (mean_longitude, mean_latitude) in degrees.
    """
    # Convert latitude and longitude from degrees to radians
    rad_coords = [(np.radians(lon), np.radians(lat)) for lon, lat in coordinates]

    # Convert spherical to Cartesian coordinates
    x = [np.cos(lat) * np.cos(lon) for lon, lat in rad_coords]
    y = [np.cos(lat) * np.sin(lon) for lon, lat in rad_coords]
    z = [np.sin(lat) for lon, lat in rad_coords]

    # Compute the mean of each Cartesian component
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    mean_z = np.mean(z)

    # Convert the mean Cartesian coordinates back to spherical coordinates
    mean_lon = np.arctan2(mean_y, mean_x)
    hyp = np.sqrt(mean_x**2 + mean_y**2)
    mean_lat = np.arctan2(mean_z, hyp)

    # Convert the result from radians to degrees
    mean_lon_deg = np.degrees(mean_lon)
    mean_lat_deg = np.degrees(mean_lat)

    return mean_lon_deg, mean_lat_deg

def explain_HroiCorr_Pres( wrf_dir, DAtime, var_name, key, x_mean):
    
    # Read mslp for all input members
    if var_name == 'PSFC':
        d_mslp = read_min_PSFC_ens( wrf_dir,DAtime, key )
    elif var_name == 'slp':
        d_mslp = read_mslp_ens( wrf_dir,DAtime,key )

    # Calculate the mean of geographic coordinates for the whole ensemble   
    ensXb_mslp_locs = list(zip(d_mslp['lon_mslp'],d_mslp['lat_mslp']))
    mean_xb_lon, mean_xb_lat = mean_geolocation( ensXb_mslp_locs )

    # Read correlations between a Tb and a column of model var
    if var_name == 'slp':
        des_path = wrf_dir+DAtime+"_d03_hroi_corr_Obs_" + obs_type +'_model_PSFC.pickle'
    else:
        des_path = wrf_dir+DAtime+"_d03_hroi_corr_Obs_" + obs_type +'_model_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        hori_corr = pickle.load( f )
    print('Shape of hori_corr: '+ str(np.shape(hori_corr)))

    # Locate the TCvitals mslp
    lon_obs = np.array( d_obs[DAtime]['lon'] )
    lat_obs = np.array( d_obs[DAtime]['lat'] )

    # Plot
    # read model attributes
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    # read lon and lat
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = UD.read_wrf_domain( wrf_file )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5,figsize=(6.5,6), dpi=300) #(6.5,6)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    if Storm == 'IRMA':
        min_corr = -0.6
        max_corr = 0.6
    else:
        min_corr = -0.8
        max_corr = 0.8
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)
    #cs = ax.flat[isub].scatter(lon,lat,5,Interp_corr[isub,:],cmap='RdBu_r',edgecolors='none',transform=ccrs.PlateCarree(),)
    cs_corr = ax.scatter(xlon,xlat,5,hori_corr,cmap='RdBu_r',edgecolors='none',vmin=min_corr,vmax=max_corr,transform=ccrs.PlateCarree(),)
    # Mark the observed location: TCvital
    ax.scatter(lon_obs, lat_obs, c='green', s=15, marker='s', edgecolors='green', transform=ccrs.PlateCarree())
    # Mark the mslp for the whole EnKF input ensemble
    if Storm == 'IRMA':
        min_mslp = 950
        max_mslp = 990
    else:
        min_mslp = 1008
        max_mslp = 1012
    cs_mslp = ax.scatter(d_mslp['lon_mslp'],d_mslp['lat_mslp'],15,d_mslp['mslp'],cmap='bone',marker='*',vmin=min_mslp,vmax=max_mslp,transform=ccrs.PlateCarree())
    #cs_mslp = ax.scatter(d_mslp['lon_mslp'],d_mslp['lat_mslp'],10,d_mslp['id'],cmap='jet',marker='*',vmin=1,vmax=60,transform=ccrs.PlateCarree())
    
    # Mark the mean of EnsXb mslp
    ax.scatter(mean_xb_lon, mean_xb_lat,c='#FF0000', s=15, marker='*',transform=ccrs.PlateCarree(),)

    # Mark the mslp in the ensemble mean 
    ax.scatter(x_mean[DAtime][0],x_mean[DAtime][1],15,x_mean[DAtime][2],cmap='bone',marker='o',vmin=min_mslp,vmax=max_mslp,transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    color_ticks1 = np.linspace(min_corr, max_corr, 5, endpoint=True)
    cbar1 = fig.colorbar(cs_corr, cax=cbaxes,fraction=0.046, pad=0.04, extend='both')
    cbar1.set_ticks( color_ticks1 )
    cbar1.ax.tick_params(labelsize=11)

    cax_below = fig.add_axes([0.25, 0.05, 0.7, 0.02])  # [left, bottom, width, height]
    color_ticks2 = np.linspace(min_mslp, max_mslp, 5, endpoint=True)
    cbar2 = fig.colorbar(cs_mslp, cax=cax_below, orientation='horizontal', extend='both')
    cbar2.set_ticks( color_ticks2 )
    cbar2.ax.tick_params(labelsize=10)
    if key == 'input':
        if var_name == 'PSFC':
            fig.text(0.018, 0.05, 'Ens Xb mpsfc (hPa)', fontweight='bold',fontsize=11)
        elif var_name == 'slp':
            fig.text(0.018, 0.05, 'Ens Xb mslp (hPa)', fontweight='bold',fontsize=11)
    else:
        if var_name == 'PSFC':
            fig.text(0.018, 0.05, 'Ens Xa mpsfc (hPa)', fontweight='bold',fontsize=11)
        elif var_name == 'slp':
            fig.text(0.018, 0.05, 'Ens Xa mslp (hPa)', fontweight='bold',fontsize=11)


    fig.text( 0.1,0.93,'Assimilated obs: '+"{0:.2f}".format((d_obs[DAtime]['obs'][0]/100))+' hPa',fontsize=10,color='green',rotation='horizontal',fontweight='bold')

    if var_name == 'PSFC':
        title_name = 'Mean of min PSFC locations in '
    elif var_name == 'slp':
        title_name = 'Mean of mslp locations in '
    if key == 'input':
        title_name = title_name + 'Ens Xb'
    else:
        title_name = title_name + 'Ens Xa'

    fig.text( 0.5,0.93,title_name,fontsize=10,color='#FF0000',rotation='horizontal',fontweight='bold')

    #subplot title
    font = {'size':12,}
    ax.set_title( 'Horizontal Corr: obs ' + obs_type +' & model PSFC', font, fontweight='bold')

    #title for all
    tt_name = Storm+': '+Exper_name
    fig.suptitle(tt_name, fontsize=10, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    gl.top_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.right_labels = False

    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 12}

    # Save the figure
    if key == 'input':
        save_des = plot_dir+'Xb_explain_'+DAtime+'_HroiCorr_obs_' + obs_type +'_model_' +  var_name + '_mslp.png'
    else:
        save_des = plot_dir+'Xa_explain_'+DAtime+'_HroiCorr_obs_' + obs_type +'_model_' +  var_name + '_mslp.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'#'/expanse/lustre/scratch/zuy121/temp_project/Pro2_PSU_MW/' #'/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'#'/expanse/lustre/projects/pen116/zuy121/Pro2_PSU_MW/'  #'/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'MARIA'
    DA = 'CONV'
    MP = 'WSM6'
    fort_v = ['obs_type','lat','lon','obs']
    sensor = 'abi_gr'

    # observation type 
    obs_assimilated = True
    if obs_assimilated:
        obs_type = 'slp' # Radiance

    # model variable
    if Storm == 'HARVEY':
        model_v = [ 'slp',]
    else:
        model_v = [ 'PSFC',]
    
    # time
    start_time_str = '201709160000'
    end_time_str = '201709170000'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    # vertical interpolation if needed
    interp_P = False
    P_range = np.arange( 995,49,-20 )
    interp_H = True
    H_range = list(np.arange(1,21,1))

    # if calculate data
    calculate_ens_data = False
    if calculate_ens_data:
        # at obs location
        at_obs_res = False
    If_save = True

    # plot
    If_plot_corr_snapshot = True
    if If_plot_corr_snapshot:
        key = 'input'

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

    # Read info of obs point of interest
    d_obs = {}
    if obs_assimilated:
        print('------------ Read obs info from enkf diagnostics fort.10000 --------------')
        for DAtime in DAtimes:
            # Read assimilated obs
            file_Diag = big_dir+Storm+'/'+Exper_name+'/run/'+DAtime+'/enkf/d03/fort.10000'
            dt = datetime.strptime(DAtime, "%Y%m%d%H%M")
            DAtime_wrf = dt.strftime("%Y-%m-%d_%H:%M")
            file_hpi = small_dir+'/Obs_input_EnKF/'+Storm+'/HPI/HPI_obs_gts_'+DAtime_wrf+':00.3DVAR'

            if obs_type == 'slp':
                d_obs[DAtime] = Diag.Find_min_slp( file_Diag, fort_v, file_hpi )

    # Calculate more data
    if calculate_ens_data:
        for DAtime in DAtimes:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            print('At '+DAtime)

            for var_name in model_v:
                start_time=time.process_time()
                if var_name == 'PSFC' or var_name == 'slp':
                    if not at_obs_res:
                        print('Finding the mslp for each EnKF member...')
                        identify_mslp_ens( wrf_dir, 'output' )
                    else:
                        print('Finding the slp at the obs location for each EnKF member...')
                        identify_slp_ens_obsLoc( d_obs[DAtime], wrf_dir, 'input' )

                end_time = time.process_time()
                print ('time needed: ', end_time-start_time, ' seconds')


    # Plot the horizontal correlations per snapshot
    if If_plot_corr_snapshot:
        print('------------ Plot the horizontal correlation --------------')

        # Read mslp of the ensemble mean 
        saved_dir = small_dir+'Clean_results/'+Storm+'/'+Exper_name+'/Data_analyze/'
        if Storm == 'HARVEY':
            x_mean = Read_mslp_EnsMean( saved_dir+'HPI_wrf_enkf_'+key+'_d03_mean.201708221200_201708231200.txt' )
        elif Storm == 'JOSE':
            x_mean = Read_mslp_EnsMean( saved_dir+'HPI_wrf_enkf_'+key+'_d03_mean.201709050000_201709060000.txt' )
        elif Storm == 'IRMA':
            x_mean = Read_mslp_EnsMean( saved_dir+'HPI_wrf_enkf_'+key+'_d03_mean.201709030000_201709040000.txt' )
        elif Storm == 'MARIA':
            x_mean = Read_mslp_EnsMean( saved_dir+'HPI_wrf_enkf_'+key+'_d03_mean.201709160000_201709170000.txt' )


        # Loop thru times
        for DAtime in DAtimes:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            print('At '+DAtime)

            for var_name in model_v:
                print('Plot horizontal correlation: '+var_name+'...')
                #var_dim = def_vardim( var_name )

                if var_name == 'PSFC':
                    plot_dir=small_dir+'Clean_results/'+Storm+'/'+Exper_name+'/Vis_analyze/hori_Corr/explain_obs_slp/'
                    explain_HroiCorr_Pres( wrf_dir, DAtime, var_name, key, x_mean)

                if var_name == 'slp':
                    plot_dir=small_dir+'Clean_results/'+Storm+'/'+Exper_name+'/Vis_analyze/hori_Corr/explain_obs_slp/'
                    explain_HroiCorr_Pres( wrf_dir, DAtime, var_name, key, x_mean)

                #if var_dim == '2D':
                #    HroiCorr_snapshot( DAtime,var_name,var_dim )
                #else:
                #    if interp_H and not interp_P:
                #        HroiCorr_snapshot( DAtime,var_name,var_dim,H_range)












import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
from wrf import getvar,interplevel
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

import Util_Vis
import Util_data as UD
import Diagnostics as Diag

def def_vardim( var_name ):
    if var_name == 'PSFC':
        return '2D'
    elif 'Q' in var_name:
        return '3D'

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
    slp_original = slp_values
    # smoothed SLP
    slp_smt_values = slp #sp.ndimage.gaussian_filter(slp, [1,1]) #[11,11]
    slp_smt_values[slp_smt_values > 1030] = np.nan
    slp_smooth = slp_smt_values
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
def identify_mslp_ens( wrf_dir ):

    mslp_info = np.zeros( shape=(num_ens,4) )
    #mslp_info = np.nan

    for ie in range(num_ens):
        id_num = ie+1
        wrf_file = wrf_dir+'wrf_enkf_input_d03_'+f"{id_num:03}" 
        d_mslp = find_mslp( Storm,wrf_file)
        #print(type(id_num))
        mslp_info[ie,1] = "{0:.3f}".format( d_mslp['lon_minslp'] )
        mslp_info[ie,2] = "{0:.3f}".format(d_mslp['lat_minslp'] )
        mslp_info[ie,3] = "{0:.3f}".format(d_mslp['mslp'] )

    # May save the mslp info
    if If_save:
        header = ['ID','Lon','Lat','mslp']
        des_path = wrf_dir+ DAtime+'_enkf_input_mslp.txt'
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

def read_mslp_ens( wrf_dir,DAtime ):

    id_ = []
    lon_mslp = []
    lat_mslp = []
    mslp = []

    with open( wrf_dir+ DAtime+'_enkf_input_mslp.txt' ) as f:
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


def read_min_PSFC_ens( wrf_dir,DAtime ):

    id_ = []
    lon_mpsfc = []
    lat_mpsfc = []
    mpsfc = []

    for ie in range(num_ens):
        id_.append( ie+1 )
        wrf_file = wrf_dir+'wrf_enkf_input_d03_'+f"{ie+1:03}"
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


def explain_HroiCorr_PSFC( wrf_dir,DAtime ):
    
    # Read mslp for all input members
    #d_mslp = read_mslp_ens( wrf_dir,DAtime )

    d_mslp = read_min_PSFC_ens( wrf_dir,DAtime )

    # Read correlations between a Tb and a column of model var
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
    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5,figsize=(6.5,6), dpi=300)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    min_corr = -0.6
    max_corr = 0.6
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)
    #cs = ax.flat[isub].scatter(lon,lat,5,Interp_corr[isub,:],cmap='RdBu_r',edgecolors='none',transform=ccrs.PlateCarree(),)
    cs_corr = ax.scatter(xlon,xlat,5,hori_corr,cmap='RdBu_r',edgecolors='none',vmin=min_corr,vmax=max_corr,transform=ccrs.PlateCarree(),)
    # Mark the observed location: TCvital
    ax.scatter(lon_obs, lat_obs, c='green', s=10, marker='s', edgecolors='green', transform=ccrs.PlateCarree())
    # Mark the mslp for the whole EnKF input ensemble
    #min_mslp = 950
    #max_mslp = 990
    #cs_mslp = ax.scatter(d_mslp['lon_mslp'],d_mslp['lat_mslp'],10,d_mslp['mslp'],cmap='bone',marker='*',vmin=min_mslp,vmax=max_mslp,transform=ccrs.PlateCarree())

    cs_mslp = ax.scatter(d_mslp['lon_mslp'],d_mslp['lat_mslp'],10,d_mslp['id'],cmap='jet',marker='*',vmin=1,vmax=60,transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    color_ticks1 = np.linspace(min_corr, max_corr, 5, endpoint=True)
    cbar1 = fig.colorbar(cs_corr, cax=cbaxes,fraction=0.046, pad=0.04, extend='both')
    cbar1.set_ticks( color_ticks1 )
    cbar1.ax.tick_params(labelsize=11)

    cax_below = fig.add_axes([0.1, 0.05, 0.8, 0.02])  # [left, bottom, width, height]
    color_ticks2 = np.linspace(min_mslp, max_mslp, 5, endpoint=True)
    cbar2 = fig.colorbar(cs_mslp, cax=cax_below, orientation='horizontal', extend='both')
    cbar2.set_ticks( color_ticks2 )
    cbar2.ax.tick_params(labelsize=10)

    #subplot title
    font = {'size':12,}
    ax.set_title( 'Horizontal Corr: obs ' + obs_type +' & model ' +  var_name, font, fontweight='bold')

    #title for all
    title_name = Storm+': '+Exper_name
    fig.suptitle(title_name, fontsize=10, fontweight='bold')

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
    save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/hori_Corr/explain_'+DAtime+'_HroiCorr_obs_' + obs_type +'_model_' +  var_name + '_id.png'
    #save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/hori_Corr/Interp_H_corr_ms_'+DAtime+'_'+var_name+'_'+sensor+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()



if __name__ == '__main__':

    big_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'#'/expanse/lustre/scratch/zuy121/temp_project/Pro2_PSU_MW/' #'/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'#'/expanse/lustre/projects/pen116/zuy121/Pro2_PSU_MW/'  #'/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'CONV'
    MP = 'WSM6'
    fort_v = ['obs_type','lat','lon','obs']
    sensor = 'abi_gr'

    # observation type 
    obs_assimilated = True
    if obs_assimilated:
        obs_type = 'slp' # Radiance

    # model variable
    model_v = [ 'PSFC',]#'QSNOW','QCLOUD','QRAIN','QICE','QGRAUP']

    # time
    start_time_str = '201709030000'
    end_time_str = '201709030000'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    #to_obs_res = False
    #ens_Interp_to_obs = False

    # vertical interpolation if needed
    interp_P = False
    P_range = np.arange( 995,49,-20 )
    interp_H = True
    H_range = list(np.arange(1,21,1))

    calculate_ens_data = False
    If_save = True

    If_plot_corr_snapshot = True

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
            if obs_type == 'slp':
                d_obs[DAtime] = Diag.Find_min_slp( file_Diag, fort_v )

    # Calculate more data
    if calculate_ens_data:
        for DAtime in DAtimes:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            print('At '+DAtime)

            for var_name in model_v:
                start_time=time.process_time()
                if var_name == 'PSFC':
                    print('Finding the mslp for each EnKF input member...')
                    identify_mslp_ens( wrf_dir )

                end_time = time.process_time()
                print ('time needed: ', end_time-start_time, ' seconds')


    # Plot the horizontal correlations per snapshot
    if If_plot_corr_snapshot:
        print('------------ Plot the horizontal correlation --------------')
        for DAtime in DAtimes:
            #Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            print('At '+DAtime)

            for var_name in model_v:
                print('Plot horizontal correlation: '+var_name+'...')
                #var_dim = def_vardim( var_name )

                if var_name == 'PSFC':
                    explain_HroiCorr_PSFC( wrf_dir, DAtime )


                #if var_dim == '2D':
                #    HroiCorr_snapshot( DAtime,var_name,var_dim )
                #else:
                #    if interp_H and not interp_P:
                #        HroiCorr_snapshot( DAtime,var_name,var_dim,H_range)












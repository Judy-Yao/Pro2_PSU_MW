#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

from numba import njit, prange
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

import Util_Vis
import Read_Obspace_IR as ROIR
import Diagnostics as Diag


# ------------------------------------------------------------------------------------------------------
#           Operation: Verify the use of data
# ------------------------------------------------------------------------------------------------------

## Demonstrated the wrf_enkf_input_d03_mean is identical to the mean of wrf_enkf_input_d03_0* (diff is noise)
# One way: use ncea to do an ensemble average (not shown)


# ------------------------------------------------------------------------------------------------------
#           Operation: Perform Operations
# ------------------------------------------------------------------------------------------------------

## Matrix multiplication
# Use Numba to speed up loops (by converting it into machine code)
@njit(parallel=True)
def mat_mult(A, B):
    # check the fact
    assert A.shape[1] == B.shape[0]
    res = np.zeros((A.shape[0], B.shape[1]), )
    for i in prange(A.shape[0]):
        for k in range(A.shape[1]):
            for j in range(B.shape[1]):
                res[i,j] += A[i,k] * B[k,j]
    return res

## Calculate ensemble variance*(N-1) members for 2D variables
@njit(parallel=True)
def var_2D(A):
    # shape of A: num_ens,xmax*ymax
    res = np.zeros( (A.shape[1],), ) 
    for j in prange(A.shape[1]):
        for i in range(A.shape[0]):
            res[j] += A[i,j] * A[i,j]
    return res


## Calculate ensemble variance*(N-1) members for 3D variables
@njit(parallel=True)
def var_3D(A):
    # shape of A: nLevel,num_ens,xmax*ymax
    res = np.zeros( (A.shape[0], A.shape[2]), ) 
    for k in prange(A.shape[2]):
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                res[i,k] += A[i,j,k,] * A[i,j,k]
    return res


## Calculate the covariance*(N-1) between Tbs (at several locations) and the model variale (3D: level,lon*lat)
@njit(parallel=True)
def cross_Tb_to3D( xb_ens,hxb_ens ):
    assert xb_ens.shape[1] == hxb_ens.shape[0]
    res = np.zeros((xb_ens.shape[0], xb_ens.shape[2], hxb_ens.shape[1]), ) # levels,model_var at lon*lat,obs (at lon*lat)
    # for each level of variable: calculate corr at each grid point
    for i in prange( xb_ens.shape[0] ):
        res[i,:,:] = mat_mult( xb_ens[i,:,:].transpose(), hxb_ens)
    return res

## Calculate the correlation between Tbs (at several locations) and the model variale (3D: level,lon*lat)
@njit(parallel=True)
def corr_Tb_to3D( cov_xb_hxb,stddev_xb,stddev_hxb ):
    assert cov_xb_hxb.shape[1] == stddev_xb.shape[1] 
    assert cov_xb_hxb.shape[2] == stddev_hxb.shape[0]
    res = np.zeros( (cov_xb_hxb.shape[0],cov_xb_hxb.shape[1],cov_xb_hxb.shape[2]), )
    for j in prange(cov_xb_hxb.shape[1]):
        for i in range(cov_xb_hxb.shape[0]):
            for k in range(cov_xb_hxb.shape[2]):
                res[i,j,k] = cov_xb_hxb[i,j,k]/stddev_xb[i,k]
                res[i,j,k] = res[i,j,k]/stddev_hxb[k]
    return res

## For each obs loc, find the nearest model grid point (good for mercator)
#@njit(parallel=True)
def nearest_axis( obs,model ):

    res = np.zeros( (obs.shape[0]), )
    res[:] = np.nan
    for io in range( obs.shape[0] ):
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





# Read simulated Tb at obs locations for all members
def read_ens_obspace( ens_Tb_file, sensor ):

    # Number of ensemble members
    num_ens = 60

    ch_obs = []
    lat_obs = []
    lon_obs = []
    Yo_obs = []

    with open( ens_Tb_file ) as f:
        next(f)
        all_lines = f.readlines()
    
    for line in all_lines:
        split_line = line.split()
        ch_obs.append( int(split_line[0]) )
        lat_obs.append( float(split_line[1]) )
        lon_obs.append( float(split_line[2]) )
        Yo_obs.append( float(split_line[3]) )

    hxb_ens = np.zeros( (num_ens,len(Yo_obs)) )
    hxb_ens[:] = np.nan
    il = 0
    for line in all_lines:
        split_line = line.split() #!!! the line didn't exist! Bug!!
        tmp = [float(it) for it in split_line[4:]]
        hxb_ens[:,il] = np.array( tmp )
        il = il+1

    lat_obs = np.array( lat_obs )
    lon_obs = np.array( lon_obs )
    ch_obs = np.array( ch_obs )
    Yo_obs = np.array( Yo_obs )
    hxb_ens = hxb_ens 

    d_obspace = {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'ch_obs':ch_obs, 'Yo_obs':Yo_obs, 'hxb_ens':hxb_ens}
    return d_obspace


# ------------------------------------------------------------------------------------------------------
#           Object: Hxb stddev ; Operation: verify
# ------------------------------------------------------------------------------------------------------

def verify_hxb( var, wrf_dir=None ):

    #file_Diag = wrf_dir+'/run/201709170000/enkf/d03/fort.10000'
    file_Diag = wrf_dir+'/run/'+DAtime+'/enkf/d03/fort.10000'
    d_obs = Diag.Find_IR( file_Diag, fort_v )
    lon_obs = list( d_obs['lon'] )
    lat_obs = list( d_obs['lat'] )

    # Read WRF domain
    #wrf_input_mean = wrf_dir+'/fc/201709160000/wrf_enkf_input_d03_mean'
    wrf_input_mean = wrf_dir+'/fc/'+DAtime+'/wrf_enkf_input_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_input_mean )
    ncdir = nc.Dataset( wrf_input_mean, 'r')
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()

    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, figsize=(6,6), dpi=400)
    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)
    if to_obs_res:
        cs = ax.scatter(lon_obs,lat_obs,3,var,cmap='rainbow',vmin=0,vmax=15,transform=ccrs.PlateCarree())
    else:
        cs = ax.scatter(xlon,xlat,1,var,cmap='rainbow',vmin=0,vmax=15,transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    color_ticks = [0,3,6,9,12,15]
    cbar = fig.colorbar(cs, cax=cbaxes)
    #cbar = fig.colorbar(cs, cax=cbaxes,ticks=color_ticks,fraction=0.046, pad=0.04)
    #bounds_str =  [ str(item) for item in color_ticks ]
    cbar.set_ticks( color_ticks )
    cbar.ax.tick_params(labelsize=12)

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(1):
        gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')

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
        gl.xlabel_style = {'size': 12}
        gl.ylabel_style = {'size': 12}

    #title 
    ax.set_title('Stddev of Hxb',fontsize=12, fontweight='bold')
    fig.suptitle(Storm+': '+Exper_name, fontsize=10, fontweight='bold')

    # Save the figure
    if to_obs_res: 
        des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/stddev_obsRes_Hxb_'+DAtime+'.png'
    else:
        des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/stddev_modelRes_Hxb_'+DAtime+'.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)
    plt.close()
    return None


# ------------------------------------------------------------------------------------------------------
#           Object: ; Operation: Calculate ensemble perturbations/variances of Xb and Hxb
# ------------------------------------------------------------------------------------------------------

# Calculate values on model space
def cal_pert_stddev_modelRes_Hxb( DAtime, sensor, Hx_dir, If_save, wrf_dir):

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    ### ------------------------- Perturbation -------------------------
    print("Read the ensemble and calculate the perturbations for Hxb in model space...")
    start_time=time.process_time()

    # Read the ensemble mean calculated from the other program: Obspace_compare_IR_txt_bin.py 
    meanYb = []
    mean_hxb = Hx_dir + "mean_model_res_d03" + DAtime + '_' +  sensor + '.txt'
    print('Reading the ensemble mean of Hxb: ', mean_hxb,'...')
    with open(mean_hxb) as f:
        next(f)
        all_lines = f.readlines()
    for line in all_lines:
        split_line = line.split()
        meanYb.append( float(split_line[3]) )
    meanYb = np.array( meanYb )
    
    hxb_ens = np.zeros( shape=[num_ens,xmax*ymax] ) # the 0 to last - 1 rows are ensemble perturbations, the last row is mean value
    hxb_ens[:] = np.nan
    #hxb_ens[num_ens,:] = meanYb
    # Read the ensemble of Hxb at model locations
    file_hxb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem0*.bin') )
    for ifile in file_hxb:
        idx = file_hxb.index( ifile )
        tmp_control = np.fromfile( ifile,dtype='<f4') # <: little endian; f: float; 4: 4 bytes
        n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
        tmp_data = tmp_control[:].reshape(n_ch+2,ymax,xmax)
        Yb_x = tmp_data[2,:,:].flatten()
        hxb_ens[idx,:] = Yb_x - meanYb

    # May save the perturbations
    if If_save:
        des_path = Hx_dir+ "Hxb_ensPert_modelRes_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( hxb_ens, f )
        f.close()
        print('Save '+des_path)

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    ### ------------------------- Variance -------------------------
    print('Calculating the ensemble variance of Hxb in model space......' )
    start = time.perf_counter()
    var_hxb = var_2D( hxb_ens[:num_ens,:] )
    end = time.perf_counter()
    print("Elapsed (after compilation) = {}s".format((end - start)))
    var_hxb = var_hxb / ( num_ens - 1 )
    stddev_hxb = np.sqrt( var_hxb )
    #print(stddev_hxb[71190])
    #print(np.std(hxb_ens[:num_ens,71190]))

    print(np.amax( stddev_hxb ))
    verify_hxb( stddev_hxb, wrf_dir )
    #print(np.std(hxb_ens[:num_ens,71190]))

    # May save the standard deviation
    if If_save:
        des_path = Hx_dir+ "Hxb_ensStddev_modelRes_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( stddev_hxb, f )
        f.close()
        print('Save '+des_path)

    return None


# Calculate values in observation space
def cal_pert_stddev_obsRes_Hxb( DAtime, sensor, Hx_dir, If_save, fort_v,  If_interp_before_pert=True, wrf_dir=None):

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    ### ------------------------- Perturbation -------------------------
    print("Read the ensemble and calculate the perturbations for Hxb in obs space...")
    start_time=time.process_time()
    
    # Read the ensemble mean calculated from the other program: Obspace_compare_IR_txt_bin.py 
    latYb = []
    lonYb = []
    meanYb = []
    mean_hxb = Hx_dir + "mean_obs_res_d03" + DAtime + '_' +  sensor + '.txt'
    print('Reading the ensemble mean of Hx: ', mean_hxb,'...')
    with open(mean_hxb) as f:
        next(f)
        all_lines = f.readlines()
    for line in all_lines:
        split_line = line.split()
        latYb.append( float(split_line[0]) )
        lonYb.append( float(split_line[1]) )
        meanYb.append( float(split_line[4]) )
    latYb = np.array( latYb )
    lonYb = np.array( lonYb )
    meanYb = np.array( meanYb )

    hxb_ens_os = np.zeros( shape=[num_ens+1,len(meanYb)] ) # the 0 to last - 1 rows are ensemble perturbations, the last row is mean value
    hxb_ens_os[:] = np.nan
    hxb_ens_os[num_ens,:] = meanYb
    
    # interp model-equivalent ens Tbs to obs location first
    if If_interp_before_pert: 
        # Read the ensemble of Hxb at obs locations
        ens_Tb_file = Hx_dir + "/Hxb_ens_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
        print('Reading the ensemble of Hx: ', ens_Tb_file,'...')
        d_obspace = read_ens_obspace( ens_Tb_file, sensor )
        hxb_ens_ol = d_obspace['hxb_ens']
        print('Shape of hxb_ens_obs: '+str(np.shape(hxb_ens_os)))
        for im in range( num_ens ):
            hxb_ens_os[im,:] = hxb_ens_ol[im,:] - meanYb
    
    # May save the perturbations
    if If_save:
        des_path = Hx_dir+ "Hxb_ensPert_obsRes_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( hxb_ens_os, f )
        f.close()
        print('Save '+des_path)

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    ### ------------------------- Variance -------------------------
    print('Calculating the ensemble variance of Hxb in obs space......' )
    start = time.perf_counter()
    var_hxb = var_2D( hxb_ens_os[:num_ens,:] )
    end = time.perf_counter()
    print("Elapsed (after compilation) = {}s".format((end - start)))
    var_hxb = var_hxb / ( num_ens-1 )
    stddev_hxb = np.sqrt( var_hxb )
    print(np.amax( stddev_hxb ))
    
    # Find the maxium stddev of hxb and its location
    max_stddev = np.amax( stddev_hxb )
    idx_max_obs =  np.where( stddev_hxb==max_stddev )[0][0] 
    lat_obs = []
    lat_obs.append( latYb[idx_max_obs] )
    lat_obs = np.array(lat_obs)
    lon_obs = []
    lon_obs.append( lonYb[idx_max_obs] )
    lon_obs = np.array( lon_obs )

    # verify the standard deviation of hxb
    verify_hxb( stddev_hxb, wrf_dir )

    # May save the standard deviation
    if If_save:
        des_path = Hx_dir+ "Hxb_ensStddev_obsRes_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( stddev_hxb, f )
        f.close()
        print('Save '+des_path)

    d_posi = {'idx':idx_max_obs,'lat_obs':lat_obs,'lon_obs':lon_obs,'meanYb_obs':meanYb[idx_max_obs]}
    print('Pre-calculated mean: '+str(meanYb[idx_max_obs]))

    return d_posi


# For a ensemble of IR at an obs location, compare its difference from the nearest model location.
def compare_obs_model( DAtime, sensor, Hx_dir, wrf_dir, d_obs ):

    # Read the ensemble of IR at obs locations
    ens_Tb_file = Hx_dir + "/Hxb_ens_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    print('Reading the ensemble of Hx: ', ens_Tb_file,'...')
    d_obspace = read_ens_obspace( ens_Tb_file, sensor )
    hxb_ens_os = d_obspace['hxb_ens']

    # Read the ensemble of IR at model locations
    # -Number of ensemble members
    num_ens = 60
    # -Dimension of the domain
    xmax = 297
    ymax = 297

    hxb_ens = np.zeros( shape=[num_ens,xmax*ymax] ) # the 0 to last - 1 rows are ensemble perturbations, the last row is mean value
    hxb_ens[:] = np.nan
    #hxb_ens[num_ens,:] = meanYb
    # Read the ensemble of Hxb at model locations
    file_hxb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem0*.bin') )
    for ifile in file_hxb:
        idx = file_hxb.index( ifile )
        tmp_control = np.fromfile( ifile,dtype='<f4') # <: little endian; f: float; 4: 4 bytes
        n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
        tmp_data = tmp_control[:].reshape(n_ch+2,ymax,xmax)
        if idx == 0:
            lon_x = tmp_data[0,:,:].flatten()
            lat_x = tmp_data[1,:,:].flatten()
        Yb_x = tmp_data[2,:,:].flatten()
        hxb_ens[idx,:] = Yb_x 

    # ------------- Search for the model grid that is closest to the obs location -----------------
    print('------- Search for the nearest model grid for each obs ------')

    # Read model lon and lat
    mean_xb = wrf_dir + '/fc/'+ DAtime + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    lon_x1d = ncdir.variables['XLONG'][0,0,:]
    lat_x1d = ncdir.variables['XLAT'][0,:,0]
    lon_x = ncdir.variables['XLONG'][0,:,:].flatten()
    lat_x = ncdir.variables['XLAT'][0,:,:].flatten()

    # loop thru all obs
    # returned is i,j in model domain
    Idx_i = nearest_axis( d_obs['lon_obs'],lon_x1d )
    Idx_j = nearest_axis( d_obs['lat_obs'],lat_x1d ) 

    # Transform to a 2d meshgrid and find the corresponding idx
    Idx_nearest = []
    for io in range(len(d_obs['lon_obs'])):
        Idx_nearest.append( int(Idx_j[io]*len(lon_x1d)+Idx_i[io]) )

    print('Idx_nearest: '+str(Idx_nearest[0]))
    # compare the difference between var on the obs location and model location
    print('Model location: '+str(lon_x[Idx_nearest[0]])+','+str(lat_x[Idx_nearest[0]]))
    print('Obs location: '+str(d_obs['lon_obs'])+','+str(d_obs['lat_obs']))
    print('Ensemble of IR at this obs location - at the nearest model location...')
    ir_obs = d_obspace['hxb_ens'][:,d_obs['idx']]
    print('Mean calculated from the 60 member at one obs location: '+str(np.mean(ir_obs)))
    ir_model = hxb_ens[:,Idx_nearest[0]]
    diff = ir_obs - ir_model
    print(np.std(ir_obs))
    print(np.std(ir_model))
    sort_diff = sorted(diff,reverse=True)
    print(np.mean(sort_diff))

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'HARVEY'
    Exper_name = 'JerryRun/IR_THO'

    v_interest = [ 'QVAPOR',]
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201708221200'
    end_time_str = '201708221200'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    to_obs_res = True
    ens_Interp_to_obs = False
    If_interp_before_pert = True
    
    If_cal_pert_stddev = True
    If_save = True
    If_plot = False
    # -------------------------------------------------------    

    if not Consecutive_times:
        IR_times = ['201709160600','201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000']#['201709160000','201709161200','201709161800',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]


    # Calculate ensemble perturbations and variances
    if If_cal_pert_stddev:
        print('------------ Calculate the ensemble perturbations --------------')
        for DAtime in IR_times:
            # Hxb
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir =  big_dir+Storm+'/'+Exper_name
            if to_obs_res:
                d_obs = cal_pert_stddev_obsRes_Hxb( DAtime, sensor, Hx_dir, If_save, fort_v, If_interp_before_pert,wrf_dir)
                compare_obs_model( DAtime, sensor, Hx_dir, wrf_dir, d_obs )
            else:
                cal_pert_stddev_modelRes_Hxb( DAtime, sensor, Hx_dir, If_save, wrf_dir)










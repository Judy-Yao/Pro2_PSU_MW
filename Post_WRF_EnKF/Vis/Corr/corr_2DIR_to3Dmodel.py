#!/usr/bin/env python3

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


## Calculate the covariance*(N-1) between Tbs (2D: lon*lat) and the model variale (3D: level,lon*lat)
@njit(parallel=True)
def cross_Tb_to3D( xb_ens,hxb_ens ):
    assert xb_ens.shape[1] == hxb_ens.shape[0]
    res = np.zeros((xb_ens.shape[0], xb_ens.shape[2], hxb_ens.shape[1]), ) # levels,model_var at lon*lat,obs (at lon*lat)
    # for each level of variable: calculate corr at each grid point
    for i in prange( xb_ens.shape[0] ):
        res[i,:,:] = mat_mult( xb_ens[i,:,:].transpose(), hxb_ens)
    return res

## Calculate the correlation between Tbs (2D: lon*lat) and the model variale (3D: level,lon*lat)
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


# ------------------------------------------------------------------------------------------------------
#           Object: ; Operation: Calculate ensemble perturbations/variances of Xb and Hxb
# ------------------------------------------------------------------------------------------------------

def cal_pert_stddev_Hxb( DAtime, sensor, Hx_dir, If_save):

    ### ------------------------- Perturbation -------------------------
    print("Read the ensemble and calculate the perturbations for Hxb...")
    start_time=time.process_time()
    
    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297
   
    hxb_ens = np.zeros( shape=[num_ens+1,xmax*ymax] ) # the 0 to last - 1 rows are ensemble perturbations, the last row is mean value
    # Read the ensemble mean calculated from the other program: Obspace_compare_IR_txt_bin.py 
    meanYb = []
    mean_hxb = Hx_dir + "mean_model_res_d03" + DAtime + '_' +  sensor + '.txt'
    print('Reading the ensemble mean of Hx: ', mean_hxb,'...')
    with open(mean_hxb) as f:
        next(f)
        all_lines = f.readlines()
    for line in all_lines:
        split_line = line.split()
        meanYb.append( float(split_line[3]) )
    hxb_ens[num_ens,:] = np.array( meanYb )
    # Read the ensemble of Hxb
    file_hxb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem0*.bin') )
    for ifile in file_hxb:
        tmp_control = np.fromfile( ifile,dtype='<f4') # <: little endian; f: float; 4: 4 bytes
        n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
        tmp_data = tmp_control[:].reshape(n_ch+2,ymax,xmax)
        Tb_all = np.array(tmp_data[2,:,:].flatten())
        idx = file_hxb.index( ifile )
        hxb_ens[idx,:] = Tb_all - hxb_ens[num_ens,:]

    # May save the perturbations
    if If_save:
        des_path = Hx_dir+ "Hxb_ens_pert_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( hxb_ens, f )
        f.close()
        print('Save '+des_path)

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    ### ------------------------- Variance -------------------------
    print('Calculating the ensemble variance of Hxb......' )
    start = time.perf_counter()
    var_hxb = var_2D( hxb_ens[:num_ens,:] )
    end = time.perf_counter()
    print("Elapsed (after compilation) = {}s".format((end - start)))
    var_hxb = var_hxb / ( num_ens-1 )
    stddev_hxb = np.sqrt( var_hxb )

    # May save the perturbations
    if If_save:
        des_path = Hx_dir+ "Hxb_ens_stddev_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( stddev_hxb, f )
        f.close()
        print('Save '+des_path)

    return None



def cal_pert_stddev_xb( DAtime, wrf_dir, var_name ):
    
    ### ------------------------- Perturbation -------------------------
    print("Read the ensemble and calculate the perturbations for Xb...")
    start_time=time.process_time()
    
    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297
    if 'Q' in var_name:
        nLevel = 42
    
    xb_ens = np.zeros( shape=[num_ens+1,nLevel,xmax*ymax] )
    # read the ensemble mean of xb
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    var = ncdir.variables[var_name][0,:,:,:]
    xb_ens[num_ens,:] = var.reshape(nLevel,xmax*ymax)
    # read the ensemble of xb
    file_xb = sorted( glob.glob(wrf_dir + '/wrf_enkf_input_d03_0*') )
    for ifile in file_xb:
        idx = file_xb.index( ifile )
        ncdir = nc.Dataset( ifile, 'r')
        var = ncdir.variables[var_name][0,:,:,:]
        xb_ens[idx,:] = var.reshape(nLevel,xmax*ymax)
    xb_ens = xb_ens.reshape(nLevel,num_ens+1,xmax*ymax)

    # May save the perturbations
    if If_save:
        des_path = wrf_dir+ "xb_d03_ens_pert_" + DAtime + '_' + var_name + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( xb_ens, f )
        f.close()
        print('Save '+des_path)

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    ### ------------------------- Variance -------------------------
    print('Calculating the ensemble variance of model variable......' )
    start = time.perf_counter()
    var_xb = var_3D( xb_ens[:,:num_ens,:] )
    end = time.perf_counter()
    print("Elapsed (after compilation) = {}s".format((end - start)))
    var_xb = var_xb / ( num_ens-1 )
    stddev_xb = np.sqrt( var_xb )

    # May save the perturbations
    if If_save:
        des_path = wrf_dir+ "xb_d03_ens_stddev_" + DAtime + '_' +  var_name + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( stddev_xb, f )
        f.close()
        print('Save '+des_path)

    return None


# ------------------------------------------------------------------------------------------------------
#           Object: ensemble correlations of Xb and Hxb
# ------------------------------------------------------------------------------------------------------

# Calculate the ensemble correlation
def calculate_corr( DAtime, Hx_dir, sensor, wrf_dir, var_name, idx_xb, idx_hxb):

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297
    if 'Q' in var_name:
        nLevel = 42

    # Read ensemble perturbations of Hxb
    des_path = Hx_dir+ "Hxb_ens_pert_" + DAtime + '_' +  sensor + '.pickle'
    with open( des_path,'rb' ) as f:
        hxb_ens = pickle.load( f )
    print('Shape of hxb_ens: '+ str(np.shape(hxb_ens)))
    # Read ensemble perturbations of xb
    des_path = wrf_dir+ "xb_d03_ens_pert_" + DAtime + '_' + var_name +  '.pickle'
    with open( des_path,'rb' ) as f:
        xb_ens = pickle.load( f )
    print('Shape of xb_ens: '+ str(np.shape(xb_ens)))
    # Read ensemble stand deviation of Hxb
    des_path = Hx_dir+ "Hxb_ens_stddev_" + DAtime + '_' +  sensor + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_hxb = pickle.load( f )
    print('Shape of stddev_hxb: '+ str(np.shape(stddev_hxb)))
    # Read ensemble standard deviation of xb
    des_path = wrf_dir+ "xb_d03_ens_stddev_" + DAtime + '_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_xb = pickle.load( f )
    print('Shape of stddev_xb: '+ str(np.shape(stddev_xb))) 

    # Calculate the covariance between Xb and Hxb
    print('Calculating the covariance between Tbs and ' + var_name + '......' )
    start = time.perf_counter()
    cov_xb_hxb = cross_Tb_to3D( xb_ens[:,:num_ens,idx_xb],hxb_ens[:num_ens,idx_hxb] )
    end = time.perf_counter()
    print("Elapsed (after compilation) of covariance calculation = {}s".format((end - start)))
    cov_xb_hxb = cov_xb_hxb / ( num_ens-1 )

    # Calculate the correlation between Xb and Hxb
    print('Calculating the correlation between Tbs and ' + var_name + '......' )
    start = time.perf_counter()
    corr_xb_hxb = corr_Tb_to3D( cov_xb_hxb,stddev_xb,stddev_hxb[idx_hxb] )
    end = time.perf_counter()
    print("Elapsed (after compilation) = {}s".format((end - start)))

    # sanity check
    assert  0 <= abs(corr_xb_hxb).all() and abs(corr_xb_hxb).all() <= 1
    
    # May save the correlations
    #if If_save:
    #    des_path = Hx_dir+ "Ens_corr_Hxb_xb_" + DAtime + '_' +  sensor + '_' + var_name +  '.pickle'
    #    f = open( des_path, 'wb' )
    #    pickle.dump( corr_xb_hxb, f )
    #    f.close()
    #    print('Save '+des_path)
    return corr_xb_hxb

# ------------------------------------------------------------------------------------------------------
#           Operation: Plot the correlation between the specified Tb and the 3D model (at specified levels)
# ------------------------------------------------------------------------------------------------------
def find_cloudy_allEns( DAtime, sensor, Hx_dir):

    print("Looking for regions that are cloudy for all ensemble members...")
    start_time=time.process_time()

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297
    # Cloudy Tb threshold
    cold_Tb = 226 # !!! after some tests 206
    idx_cloud_ens = [ [] for i in range(num_ens) ] 
    # Collect indices that are cloudy for all ensemble members
    file_hxb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem0*.bin') )
    for ifile in file_hxb:
        tmp_control = np.fromfile( ifile,dtype='<f4') # <: little endian; f: float; 4: 4 bytes
        n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
        tmp_data = tmp_control[:].reshape(n_ch+2,ymax,xmax)
        Tb_all = np.array(tmp_data[2,:,:].flatten())
        idx = file_hxb.index( ifile )
        tmp_idx = np.where(Tb_all<cold_Tb)
        idx_cloud_ens[idx] = tmp_idx[0]
    # Find the common index (or indices)
    common_ini = list(set( idx_cloud_ens[0] ).intersection( idx_cloud_ens[1] ))
    for i in range(1,num_ens-1):
        if i == 1:
            common_idx = common_ini
        common_idx = list(set( common_idx ).intersection( idx_cloud_ens[i+1] ))
    print('Indices of cloudy region:'+str(common_idx))
    return common_idx


def corr_cloud_3Dmodel( DAtime, Hx_dir, sensor, wrf_dir, small_dir, Storm, Exper_name, var_name ):
 
    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    # Calculate the correlation
    idx_xb = np.arange(xmax*ymax)
    cloud_hxb = find_cloudy_allEns( DAtime, sensor, Hx_dir)
    corr_xb_hxb = calculate_corr( DAtime, Hx_dir, sensor, wrf_dir, var_name, idx_xb, cloud_hxb)
    #des_path = Hx_dir+ "Ens_corr_Hxb_xb_" + DAtime + '_' +  sensor + '_' + var_name +  '.pickle'
    #with open( des_path,'rb' ) as f:
    #    corr_xb_hxb = pickle.load( f )
    print('Shape of corr_xb_hxb: '+ str(np.shape(corr_xb_hxb)))

    # Read model attributes
    nLevel = corr_xb_hxb.shape[0] 
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    # pressure levels
    PB = ncdir.variables['PB'][0,:,:,:]
    P = ncdir.variables['P'][0,:,:,:]
    P_hpa = (PB + P)/100 # 0 dimension: bottom to top
    P_hpa = P_hpa.reshape(nLevel,xmax*ymax)
    # lon and lat
    lon = ncdir.variables['XLONG'][0,:,:].flatten()
    lat = ncdir.variables['XLAT'][0,:,:].flatten()
    print([lon[idx] for idx in cloud_hxb])
    print([lat[idx] for idx in cloud_hxb])

    # Specify pressure levels of interest
    P_of_interest = [100,850,500]

    #f_interp = interpolate.interp1d( P_hpa[iz,:].reshape( np.shape, var_mean )
    #    var_interp = f_interp( P_of_interest ) 



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # Configuration
    Storm = 'MARIA'
    Exper_name = 'IR-J_DA+J_WRF+J_init-SP-intel17-THO-24hr-hroi900'
    v_interest = [ 'QVAPOR',]
    sensor = 'abi_gr'
    ch_list = ['8',]
    start_time_str = '201709160000'
    end_time_str = '201709160000'
    Consecutive_times = True
    If_cal_pert_stddev = False
    If_cal_corr = False
    If_save = True
    If_plot = True
    

    if not Consecutive_times:
        IR_times = ['201709050000','201709050600','201709051200',]#'201708221800','201708230000','201708230600','201708231200']
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
            cal_pert_stddev_Hxb( DAtime, sensor, Hx_dir, If_save)
            # Xb
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in v_interest:
                cal_pert_stddev_xb( DAtime, wrf_dir, var_name )

    # Calculate correlations
    if If_cal_corr:
        print('------------ Calculate the correlation between IR Tb and model variables --------------')
        for DAtime in IR_times:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in v_interest:
                calculate_corr( DAtime, Hx_dir, sensor, ch_list, wrf_dir, var_name, If_save )
   
    # Plot the correlations
    if If_plot:
        print('------------ Plot the correlation --------------')
        for DAtime in IR_times:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            print('At '+DAtime)
            for var_name in v_interest:
                print('Plot '+var_name+'...')
                corr_cloud_3Dmodel( DAtime, Hx_dir, sensor, wrf_dir, small_dir, Storm, Exper_name, var_name )










    

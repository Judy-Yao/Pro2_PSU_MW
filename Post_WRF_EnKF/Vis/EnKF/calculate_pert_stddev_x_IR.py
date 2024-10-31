
from numba import njit, prange
import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
import matplotlib
from wrf import getvar,interplevel
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
import Util_data as UD
import Read_Obspace_IR as ROIR
import Diagnostics as Diag
#import matlab.engine


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

def vertical_interp( ncdir,array,levels,interp_H=None,interp_P=None):

    if interp_H and not interp_P:
        #interp_arr = np.zeros( [len(H_of_interest),len(idx_xb)] )
        z = getvar(ncdir, 'z', units='km')
        interp_arr = np.zeros( (len(levels),z.shape[1],z.shape[2]) )
        array =  array.reshape( (z.shape) )
        start_time=time.process_time()
        for ih in levels:
            interp_arr[levels.index(ih),:,:] = interplevel(array, z, ih)
        end_time = time.process_time()
        print ('time needed for the interpolation: ', end_time-start_time, ' seconds')
        print('Min of interpolated_arr: '+str(np.amin( interp_arr )))
        print('Max of interpolated_arr: '+str(np.amax( interp_arr )))
        return interp_arr
    elif interp_P and not interp_H:
        pres = getvar(ncdir, 'pres', units='hPa')
        interp_arr = np.zeros( (len(levels),pres.shape[1],pres.shape[2]) )
        array =  array.reshape( (pres.shape) )
        start_time=time.process_time()
        for ih in levels:
            interp_arr[levels.index(ih),:,:] = interplevel(array, pres, ih)
        end_time = time.process_time()
        print ('time needed for the interpolation: ', end_time-start_time, ' seconds')
        print('Min of interpolated_arr: '+str(np.amin( interp_arr )))
        print('Max of interpolated_arr: '+str(np.amax( interp_arr )))
        return interp_arr
    else:
        pass



# ------------------------------------------------------------------------------------------------------
#           Object: ; Operation: Write or read Tbs at obs locations
# ------------------------------------------------------------------------------------------------------
# Interpolate IR Tbs in model resolution to assimilated obs locations and Write it to a txt file
def interp_simu_to_obs_matlab_ens( d_obs, Hx_dir, sensor, ch_list,  DAtime ):

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    start_time=time.process_time()

    # Read IR_obs attributes
    IR_obs = d_obs['obs']
    lon_obs = d_obs['lon']
    lat_obs = d_obs['lat']
 
    hxb_ens = np.zeros( shape=[num_ens,len(d_obs['obs_type'])] )
    file_hxb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem0*.bin') )
    # Get common variables: lat,lon
    tmp_control = np.fromfile( file_hxb[0],dtype='<f4')
    n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
    tmp_data = tmp_control[:].reshape(n_ch+2,ymax,xmax)
    lon_x =  list(tmp_data[0,:,:].flatten()) #longitude
    lat_x = list(tmp_data[1,:,:].flatten())  #latitude
    ## Interpolate to obs location
    # start a matlab process
    eng = matlab.engine.start_matlab() 
    for ifile in file_hxb:
        print(ifile)
        idx = file_hxb.index( ifile )
        tmp_control = np.fromfile( ifile,dtype='<f4') # <: little endian; f: float; 4: 4 bytes
        n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
        tmp_data = tmp_control[:].reshape(n_ch+2,ymax,xmax)
        Yb_x = list(tmp_data[2,:,:].flatten())
        
        for ich in ch_list:
            if ifile == file_hxb[0]:
                Ch_all = np.full( np.shape(IR_obs), ich )
            if sum(np.isnan(Yb_x)) != 0:
                warnings.warn('NaN value exists in Yb_x!') 
            # interpolate simulated Tbs to obs location
            mYb_obspace = eng.griddata(matlab.double(lon_x), matlab.double(lat_x), matlab.double(Yb_x), matlab.double(lon_obs.tolist()), matlab.double(lat_obs.tolist()) )
            #mYb_obspace = eng.griddata(matlab.double(lon_x), matlab.double(lat_x), matlab.double(Yb_x), matlab.double(lon_obs), matlab.double(lat_obs) )
            Yb_obspace = np.array(mYb_obspace._data)
            hxb_ens[idx,:] = ["{0:.4f}".format(item) for item in Yb_obspace]
    # end the matlab process
    eng.quit()

    # Stack each list into an array
    all_attrs = np.row_stack( (Ch_all,lat_obs,lon_obs,IR_obs,hxb_ens)  )
    all_attrs = all_attrs.transpose()
    print('Shape of all_attrs: '+str(np.shape(all_attrs)))
    # ---- Write to file and save it to the disk ----
    header = ['Ch_num','Lat','Lon','Tb_obs','Hxb']
    file_name = Hx_dir + "/Hxb_ens_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    with open(file_name,'w') as f:
        # Add header 
        f.write('\t'.join( item.rjust(6) for item in header ) + '\n' )
        # Write the record to the file serially
        len_records = np.shape( all_attrs )[0]
        for irow in range( len_records ):
            irecord =  [str(item) for item in all_attrs[irow,:] ]
            f.write('\t'.join( item.rjust(6) for item in irecord ) + '\n')
    print('Save '+file_name)
    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    return None 

# Read simulated Tb at obs locations for all members
def read_ens_obspace( ens_Tb_file, sensor ):

    # --------------------------------------------------------------------------------
    # Data inside ens_Tb_file is like 
    # each row (from left to right): at a location with one ch_num,
    #      ch_num, lat, lon, Tb_obs, Tb_simu_mem1, Tb_simu_mem2, ..., Tb_simu_mem60
    # each column (from top to down): iterate over all channels and all locations 
    # --------------------------------------------------------------------------------

    # Number of ensemble members
    num_ens = 60

    ch_obs = []
    lat_obs = []
    lon_obs = []
    Yo_obs = []

    # Read member-invariate variables:ch_num, lat, lon, Tb_obs
    with open( ens_Tb_file ) as f:
        next(f)
        all_lines = f.readlines()
    
    for line in all_lines:
        split_line = line.split()
        ch_obs.append( int(split_line[0]) )
        lat_obs.append( float(split_line[1]) )
        lon_obs.append( float(split_line[2]) )
        Yo_obs.append( float(split_line[3]) )

    # Read member-variate variable: simulated Tb 
    hxb_ens = np.zeros( (num_ens,len(Yo_obs)) )
    hxb_ens[:] = np.nan
    il = 0
    for line in all_lines:
        split_line = line.split() #!!! the line didn't exist! Bug!!
        tmp = [float(it) for it in split_line[4:]] # for each line: read the content (Tb_ens of a sample) after ch_num, lat, lon, Tb_obs
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
    mean_hxb = Hx_dir + "mean_model_res_d03_" + DAtime + '_' +  sensor + '.txt'
    print('Reading the ensemble mean of Hxb: ', mean_hxb,'...')
    with open(mean_hxb) as f:
        next(f)
        all_lines = f.readlines()
    for line in all_lines:
        split_line = line.split()
        meanYb.append( float(split_line[3]) )
    meanYb = np.array( meanYb )
    
    hxb_ens = np.zeros( shape=[num_ens+1,xmax*ymax] ) # the 0 to last - 1 rows are ensemble perturbations, the last row is mean value
    hxb_ens[:] = np.nan
    hxb_ens[num_ens,:] = meanYb
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
    var_hxb = var_hxb / ( num_ens-1 )
    stddev_hxb = np.sqrt( var_hxb )

    # Check if there are any NaN values using assert
    assert not np.isnan(stddev_hxb).any()

    #verify_hxb( d_obs_ens, DAtime, var, wrf_dir )

    # May save the standard deviation
    if If_save:
        des_path = Hx_dir+ "Hxb_ensStddev_modelRes_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( stddev_hxb, f )
        f.close()
        print('Save '+des_path)

    return None


# Calculate values in observation space
def cal_pert_stddev_obsRes_Hxb( DAtime, sensor, Hx_dir, If_save, fort_v,  wrf_dir=None):

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    ### ------------------------- Perturbation -------------------------
    print("Read the ensemble and calculate the perturbations for Hxb in obs space...")
    start_time=time.process_time()
    
    # Read the ensemble mean calculated from the other program: Obspace_compare_IR_txt_bin.py 
    meanYb = []
    mean_hxb = Hx_dir + "mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    print('Reading the ensemble mean of Hx: ', mean_hxb,'...')
    with open(mean_hxb) as f:
        next(f)
        all_lines = f.readlines()
    for line in all_lines:
        split_line = line.split()
        meanYb.append( float(split_line[4]) )
    meanYb = np.array( meanYb )

    hxb_ens_os_pert = np.zeros( shape=[num_ens+1,len(meanYb)] ) # the 0 to last - 1 rows are ensemble perturbations, the last row is mean value
    hxb_ens_os_pert[:] = np.nan
    hxb_ens_os_pert[num_ens,:] = meanYb
    
    # Read the ensemble of Hxb at obs locations
    ens_Tb_file = Hx_dir + "/Hxb_ens_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    print('Reading the ensemble of Hx: ', ens_Tb_file,'...')
    d_obspace = read_ens_obspace( ens_Tb_file, sensor )
    print('Shape of hxb_ens_obs: '+str(np.shape(d_obspace['hxb_ens'])))

    # Calculate the perturbation of Hxb
    for im in range( num_ens ):
        hxb_ens_os_pert[im,:] = d_obspace['hxb_ens'][im,:] - meanYb
 
    # May save the perturbations
    if If_save:
        des_path = Hx_dir+ "Hxb_ensPert_obsRes_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( hxb_ens_os_pert, f )
        f.close()
        print('Save '+des_path)

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    ### ------------------------- Variance -------------------------
    print('Calculating the ensemble variance of Hxb in obs space......' )
    start = time.perf_counter()
    var_hxb = var_2D( hxb_ens_os_pert[:num_ens,:] )
    end = time.perf_counter()
    print("Elapsed (after compilation) = {}s".format((end - start)))
    var_hxb = var_hxb / ( num_ens-1 )
    stddev_hxb = np.sqrt( var_hxb )

    # Check if there are any NaN values using assert
    assert not np.isnan(stddev_hxb).any()

    # verify the standard deviation of hxb
    verify_hxb( d_obspace, DAtime, stddev_hxb,  wrf_dir )

    # May save the standard deviation
    if If_save:
        des_path = Hx_dir+ "Hxb_ensStddev_obsRes_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( stddev_hxb, f )
        f.close()
        print('Save '+des_path)

    return None


def cal_pert_stddev_xb( DAtime, wrf_dir, var_name, If_save, dim=None ):
  
    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297
    if dim == '3D':
        nLevel = 42
    elif dim == '2D':
        nLevel = 1

    ### ------------------------- Perturbation -------------------------
    print("Read the ensemble and calculate the perturbations for Xb...")
    start_time=time.process_time()
    
    xb_ens_pert = np.zeros( shape=[nLevel,num_ens+1,xmax*ymax] )
    xb_ens_pert[:] = np.nan
    # read the ensemble mean of xb
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    if var_name == 'W':
        var = ncdir.variables[var_name][0,:,:,:]
        var = var.reshape(nLevel+1,xmax*ymax)
        var_half = (var[:-1]+var[1:])/2
        var = np.ma.getdata( var_half )
    elif var_name == 'rt_vo':
        xb_avo = getvar( ncdir, 'avo') # Absolute vorticity, units: 10-5 s-1 
        xb_coriolis_sin = ncdir.variables['F'][0,:,:]/1e-5 #units: 10-5 s-1
        xb_rtvo = np.array(xb_avo - xb_coriolis_sin) #relative vorticity, units: 10-5 s-1 
        var = xb_rtvo.reshape(nLevel,xmax*ymax)
    elif dim == '3D':
        var = ncdir.variables[var_name][0,:,:,:]
        var = var.reshape(nLevel,xmax*ymax)
    elif dim == '2D':
        var = ncdir.variables[var_name][0,:,:]
        var = var.reshape(nLevel,xmax*ymax)

    xb_ens_pert[:,num_ens,:] = var.reshape(nLevel,xmax*ymax)
    # read the ensemble of xb
    file_xb = sorted( glob.glob(wrf_dir + '/wrf_enkf_input_d03_0*') )
    for ifile in file_xb:
        idx = file_xb.index( ifile )
        ncdir = nc.Dataset( ifile, 'r')
        if var_name == 'W':
            var = ncdir.variables[var_name][0,:,:,:]
            var = var.reshape(nLevel+1,xmax*ymax)
            var_half = (var[:-1]+var[1:])/2
            var = np.ma.getdata( var_half )
        elif var_name == 'rt_vo':
            xb_avo = getvar( ncdir, 'avo') # Absolute vorticity, units: 10-5 s-1 
            xb_coriolis_sin = ncdir.variables['F'][0,:,:]/1e-5 #units: 10-5 s-1
            xb_rtvo = np.array(xb_avo - xb_coriolis_sin) #relative vorticity, units: 10-5 s-1 
            var = xb_rtvo.reshape(nLevel,xmax*ymax)
        elif dim == '3D':
            var = ncdir.variables[var_name][0,:,:,:]
            var = var.reshape(nLevel,xmax*ymax)
        elif dim == '2D':
            var = ncdir.variables[var_name][0,:,:]
            var = var.reshape(nLevel,xmax*ymax)

        xb_ens_pert[:,idx,:] = var - xb_ens_pert[:,num_ens,:]
    
    # Check if there are any NaN values using assert
    assert not np.isnan(xb_ens_pert).any()

    # May save the perturbations
    if If_save:
        des_path = wrf_dir+ "xb_d03_"+dim+"_ensPert_" + DAtime + '_' + var_name + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( xb_ens_pert, f )
        f.close()
        print('Save '+des_path)

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    ### ------------------------- Variance -------------------------
    print('Calculating the ensemble variance of model variable......' )
    start = time.perf_counter()
    var_xb = var_3D( xb_ens_pert[:,:num_ens,:] ) # var_xb: nLevel,xmax*ymax
    end = time.perf_counter()
    print("Elapsed (after compilation) = {}s".format((end - start)))
    var_xb = var_xb / ( num_ens-1 )
    stddev_xb = np.sqrt( var_xb )

    # May save the standard deviation
    if If_save:
        des_path = wrf_dir+ "xb_d03_"+dim+"_ensStddev_" + DAtime + '_' +  var_name + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( stddev_xb, f )
        f.close()
        print('Save '+des_path)

    return None


# ------------------------------------------------------------------------------------------------------
#           Object: ensemble correlations of Xb and Hxb
# ------------------------------------------------------------------------------------------------------

# # Calculate the ensemble correlation
# def calculate_corr( DAtime, Hx_dir, sensor, wrf_dir, var_name, idx_xb, idx_hxb):
# 
#     # Read ensemble perturbations of xb
#     des_path = wrf_dir+ "xb_d03_3D_ensPert_" + DAtime + '_' + var_name +  '.pickle'
#     with open( des_path,'rb' ) as f:
#         xb_ens = pickle.load( f )
#     print('Shape of xb_ens: '+ str(np.shape(xb_ens)))
#     # Read ensemble standard deviation of xb
#     des_path = wrf_dir+ "xb_d03_3D_ensStddev_" + DAtime + '_' +  var_name + '.pickle'
#     with open( des_path,'rb' ) as f:
#         stddev_xb = pickle.load( f )
#     print('Shape of stddev_xb: '+ str(np.shape(stddev_xb)))
#     # Read ensemble perturbations of Hxb 
#     if to_obs_res:
#         des_path = Hx_dir+ "Hxb_ensPert_obsRes_" + DAtime + '_' +  sensor + '.pickle'
#     else:
#         des_path = Hx_dir+ "Hxb_ensPert_modelRes_" + DAtime + '_' +  sensor + '.pickle'
#     with open( des_path,'rb' ) as f:
#         hxb_ens = pickle.load( f )
#     print('Shape of hxb_ens: '+ str(np.shape(hxb_ens)))
#     # Read ensemble stand deviation of Hxb
#     if to_obs_res:
#         des_path = Hx_dir+ "Hxb_ensStddev_obsRes_" + DAtime + '_' +  sensor + '.pickle'
#     else:
#         des_path = Hx_dir+ "Hxb_ensStddev_modelRes_" + DAtime + '_' +  sensor + '.pickle'
#     with open( des_path,'rb' ) as f:
#         stddev_hxb = pickle.load( f )
#     print('Shape of stddev_hxb: '+ str(np.shape(stddev_hxb)))
# 
#     # Calculate the covariance between Xb and Hxb
#     print('Calculating the covariance between Tbs and ' + var_name + '......' )
#     start = time.perf_counter()
#     cov_xb_hxb = cross_Tb_to3D( xb_ens[:,:num_ens,idx_xb],hxb_ens[:num_ens,idx_hxb] )
#     end = time.perf_counter()
#     print("Elapsed (after compilation) of covariance calculation = {}s".format((end - start)))
#     cov_xb_hxb = cov_xb_hxb / ( num_ens-1 )
# 
#     # Calculate the correlation between Xb and Hxb
#     print('Calculating the correlation between Tbs and ' + var_name + '......' )
#     start = time.perf_counter()
#     corr_xb_hxb = corr_Tb_to3D( cov_xb_hxb,stddev_xb,stddev_hxb[idx_hxb] )
#     end = time.perf_counter()
#     print("Elapsed (after compilation) = {}s".format((end - start)))
# 
#     # sanity check
#     assert  0 <= abs(corr_xb_hxb).all() and abs(corr_xb_hxb).all() <= 1
#     
#     # May save the correlations
#     if If_save:
#         if to_obs_res:
#             des_path = Hx_dir+ "Hcorr_Hxb_xb_obsRes_" + DAtime + '_' +  sensor + '_' + var_name +  '.pickle'
#         else:
#             des_path = Hx_dir+ "Hcorr_Hxb_xb_modelRes_" + DAtime + '_' +  sensor + '_' + var_name +  '.pickle'
#         f = open( des_path, 'wb' )
#         pickle.dump( corr_xb_hxb, f )
#         f.close()
#         print('Save '+des_path)
#     return corr_xb_hxb

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'CONV'
    MP = 'THO'

    v_interest = [ 'QCLOUD','QRAIN','QICE','QGRAUP'] #[ 'PSFC']
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709030000'
    end_time_str = '201709030000'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    to_obs_res = False
    ens_Interp_to_obs = False
   
    If_cal_pert_stddev = True
    If_save = True
    If_plot = False
    # -------------------------------------------------------    

    Exper_name = UD.generate_one_name( Storm,DA,MP )

    if not Consecutive_times:
        IR_times = ['201709160600','201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000']#['201709160000','201709161200','201709161800',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]


    # Interpolate simulated Tb at model resolution to obs resolution
    if ens_Interp_to_obs:
        print('------- For all members, interpolate Hx in model resolution to obs location ------')
        for DAtime in IR_times:
            # Read assimilated obs 
            file_Diag = big_dir+Storm+'/'+Exper_name+'/run/'+DAtime+'/enkf/d03/fort.10000'
            d_obs = Diag.Find_IR( file_Diag, fort_v )
            # Interpolate
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            interp_simu_to_obs_matlab_ens( d_obs, Hx_dir, sensor, ch_list, DAtime )

    # Calculate ensemble perturbations and variances
    if If_cal_pert_stddev:
        print('------------ Calculate the ensemble perturbations --------------')
        for DAtime in IR_times:
            # Hxb
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir =  big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            if to_obs_res:
                cal_pert_stddev_obsRes_Hxb( DAtime, sensor, Hx_dir, If_save, fort_v, wrf_dir)
            else:
                cal_pert_stddev_modelRes_Hxb( DAtime, sensor, Hx_dir, If_save, wrf_dir)
            # Xb
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in v_interest:
                cal_pert_stddev_xb( DAtime, wrf_dir, var_name, If_save, '2D')

    # Calculate correlations
    #if If_cal_corr:
    #    print('------------ Calculate the correlation between IR Tb and model variables --------------')
    #    for DAtime in IR_times:
    #        Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
    #        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
    #        for var_name in v_interest:
    #            calculate_corr( DAtime, Hx_dir, sensor, ch_list, wrf_dir, var_name, If_save )
   









    

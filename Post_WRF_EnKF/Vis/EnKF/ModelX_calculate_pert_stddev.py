
from numba import njit, prange
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import matplotlib
from wrf import getvar,interplevel
from scipy import interpolate
matplotlib.use("agg")
import time
import pickle

import Util_data as UD


# Define constants
R_D = 287.0
Cpd=7.0*R_D/2.0
P1000MB=100000.0

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
        #print ('time needed for the interpolation: ', end_time-start_time, ' seconds')
        #print('Min of interpolated_arr: '+str(np.amin( interp_arr )))
        #print('Max of interpolated_arr: '+str(np.amax( interp_arr )))
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
#        Operation: set time range 
# ------------------------------------------------------------------------------------------------------

def set_time_for_MakeEns( Storm ):

    if Storm == 'HARVEY':
        return ['201708220000',]
    elif Storm == 'IRMA':
        return ['201709021200',]
    elif Storm == 'JOSE':
        return ['201709041200',]
    elif Storm == 'MARIA':
        return ['201709151200',]

def set_time_for_cyclings( Storm, start_time_str,end_time_str):

    if not Consecutive_times:
        times = ['']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return times


# ------------------------------------------------------------------------------------------------------
#        Operation: read and process variables from WRF files
# ------------------------------------------------------------------------------------------------------
# Read regular variables
def Read_unstagger_WRF( ncdir, var_name ):
   
    # 3D variable
    if var_name == 'U': # m/s
        var = ncdir.variables[var_name][0,:,:,:]
        var_mass = (var[:,:,:-1] + var[:,:,1:])/2
        return var_mass
    elif var_name == 'V': # m/s
        var = ncdir.variables[var_name][0,:,:,:]
        var_mass = (var[:,:-1,:] + var[:,1:,:])/2
        return var_mass
    elif var_name == 'W': # m/s
        var = ncdir.variables[var_name][0,:,:,:]
        var_mass = (var[:-1,:,:] + var[1:,:,:])/2
        return var_mass
    elif var_name == 'T': # K 
        pres = (ncdir.variables['P'][0,:,:,:]+ncdir.variables['PB'][0,:,:,:])
        Theta = ncdir.variables[var_name][0,:,:,:]
        var = (Theta+300)*(pres/P1000MB) ** (R_D/Cpd)
        return var
    elif var_name == 'P': # hp
        var = (ncdir.variables['P'][0,:,:,:]+ncdir.variables['PB'][0,:,:,:])
        return var
    elif 'Q' in var_name: # kg/kg
        var = ncdir.variables[var_name][0,:,:,:] 
        return var
    # 2D variable
    elif var_name == 'PSFC':
        return ncdir.variables[var_name][0,:,:]

# ------------------------------------------------------------------------------------------------------
#           Object: Xb; Operation: Calculate ensemble perturbations/variances of Xb 
# ------------------------------------------------------------------------------------------------------
def Cycling_pert_stddev_xb( DAtime, wrf_dir, var_name, If_save, dim=None ):
  
    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297
    if dim == '2D':
        nLevel = 1
    elif dim == '3D':
        if if_interp:
            if interp_H:
                nLevel = len(H_range)
            else:
                nLevel = len(P_range)
        else:
            nLevel = 42

    ### ------------------------- Perturbation -------------------------
    print("Read the ensemble and calculate the perturbations for Xb...")
    start_time=time.process_time()
    
    xb_ens_pert = np.zeros( shape=[nLevel,num_ens+1,xmax*ymax] )
    xb_ens_pert[:] = np.nan
    # read the ensemble mean of xb
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    var = Read_unstagger_WRF( ncdir, var_name )
    if dim == '3D' and if_interp:
        if interp_H:
            var = vertical_interp( ncdir,var,H_range,interp_H,interp_P)
        else:
            pass
    var = var.reshape(nLevel,xmax*ymax)
    xb_ens_pert[:,num_ens,:] = var.reshape(nLevel,xmax*ymax)
    
    # read the ensemble members of xb
    file_xb = sorted( glob.glob(wrf_dir + '/wrf_enkf_input_d03_0*') )
    for ifile in file_xb:
        idx = file_xb.index( ifile )
        ncdir = nc.Dataset( ifile, 'r')
        var = Read_unstagger_WRF( ncdir, var_name )
        if dim == '3D' and if_interp:
            if interp_H:
                var = vertical_interp( ncdir,var,H_range,interp_H,interp_P)
            else:
                pass
        var = var.reshape(nLevel,xmax*ymax)
        xb_ens_pert[:,idx,:] = var - xb_ens_pert[:,num_ens,:]
    
    # Check if there are any NaN values using assert
    #assert not np.isnan(xb_ens_pert).any()

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


def CV3_pert_stddev_xb( DAtime, wrf_dir, var_name, If_save, dim=None ):

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297
    if dim == '2D':
        nLevel = 1
    elif dim == '3D':
        nLevel = 42

    ### ------------------------- Perturbation -------------------------
    print("Calculate the perturbations for CV3-generated ensemble...")
    start_time=time.process_time()

    xb_ens_pert = np.zeros( shape=[nLevel,num_ens+1,xmax*ymax] )
    xb_ens_pert[:] = np.nan
    # read the GFS field
    mean_xb = wrf_dir + '/wrfinput_d03'
    ncdir = nc.Dataset( mean_xb, 'r')
    var = Read_unstagger_WRF( ncdir, var_name )
    var = var.reshape(nLevel,xmax*ymax)
    xb_ens_pert[:,num_ens,:] = var.reshape(nLevel,xmax*ymax)
    
    # read the CV3-generated ensemble
    file_xb = sorted( glob.glob(wrf_dir + '/wrfinput_d03_0*') )
    for ifile in file_xb:
        idx = file_xb.index( ifile )
        ncdir = nc.Dataset( ifile, 'r')
        var = Read_unstagger_WRF( ncdir, var_name )
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

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR-WSM6Ens'
    MP = 'THO'
    # variables of interest
    v_interest = ['PSFC','QSNOW'] #['U','V','W','T','P','QVAPOR','PSFC'] #[ 'PSFC']
    # Which stage
    MakeEns = False # Ens generated from perturbing GFS field
    if MakeEns:
        times = set_time_for_MakeEns( Storm )
    else:
        start_time_str = '201709030600'
        end_time_str = '201709030600'
        Consecutive_times = True
        times = set_time_for_cyclings( Storm, start_time_str,end_time_str)
    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297
 
    # vertical interpolation if needed
    if_interp = True
    if if_interp:
        interp_P = False
        P_range = np.arange( 995,49,-20 )
        interp_H = True
        H_range = list(np.arange(1,21,1))

    If_save = True
    # -------------------------------------------------------    
    Exper_name = UD.generate_one_name( Storm,DA,MP )

    if not MakeEns:
        print('------------ Calculate the ensemble perturbations for EnKF cyclings--------------')
        for DAtime in times:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in v_interest:
                var_dim = UD.def_vardim( var_name )
                Cycling_pert_stddev_xb( DAtime, wrf_dir, var_name, If_save, var_dim)

    else:
        print('------------ Calculate the ensemble perturbations before 12-hr spinup--------------')
        for tt in times:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+tt+'/'
            for var_name in v_interest:
                print('Working on '+var_name+' ......')
                var_dim = UD.def_vardim( var_name )
                CV3_pert_stddev_xb( tt, wrf_dir, var_name, If_save,var_dim )






    

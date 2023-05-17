#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7


#from numba import njit, prange
import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
import math
#import matplotlib
#from scipy import interpolate
#matplotlib.use("agg")
#import matplotlib.ticker as mticker
#from matplotlib import pyplot as plt
#from matplotlib import colors
#from cartopy import crs as ccrs
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from mpl_toolkits.axes_grid1 import make_axes_locatable
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

#@njit(parallel=True)
def mat_mult(A, B):
    # check the fact
    assert A.shape[1] == B.shape[0]
    res = np.zeros((A.shape[0], B.shape[1]), )
    for i in prange(A.shape[0]):
        for k in range(A.shape[1]):
            for j in range(B.shape[1]):
                res[i,j] += A[i,k] * B[k,j]
    return res

## Calculate the corr between Tbs (2D: lon*lat) and the model variale (3D: level,lon*lat)

#@njit(parallel=True)
def corr_Tb_to3D( xb_ens,hxb_ens ):
    assert xb_ens.shape[1] == hxb_ens.shape[0]
    res = np.zeros((xb_ens.shape[0], xb_ens.shape[2], hxb_ens.shape[1]), )
    for i in prange( xb_ens.shape[0] ):
        res[i,:,:] = mat_mult( xb_ens[i,:,:].transpose(), hxb_ens)
    return res

# ------------------------------------------------------------------------------------------------------
#           Object: ; Operation:
# ------------------------------------------------------------------------------------------------------

def calculate_pert_Hxb( DAtime, sensor, Hx_dir, ch_list, If_save):

    print("Initiate the function to read the ensemble and calculate the perturbations for Hxb...")
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
    return hxb_ens

def calculate_pert_xb( DAtime, wrf_dir, var_name ):
    
    print("Initiate the function to read the ensemble and calculate the perturbations for Xb...")
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
        des_path = wrf_dir+ "xb_d03_ens_pert_" + DAtime + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( xb_ens, f )
        f.close()
        print('Save '+des_path)

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')
    return xb_ens


# ------------------------------------------------------------------------------------------------------
#           Object: Perform Operations
# ------------------------------------------------------------------------------------------------------




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
    If_calculate = True
    If_save = True
    

    if not Consecutive_times:
        IR_times = ['201709050000','201709050600','201709051200',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Calculate correlations
    if If_calculate:
        print('------------ Calculate the correlation between IR Tb and model variables --------------')
        for DAtime in IR_times:
            #Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            # Calculate the corr
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            #calculate_corr( DAtime, sensor, Hx_dir, ch_list, If_save )
            #lala = calculate_pert_Hxb( DAtime, sensor, Hx_dir, ch_list, If_save)
            lala = calculate_pert_xb( DAtime, wrf_dir, 'QVAPOR' ) 
 
    #start = time.perf_counter()
    #corr = corr_Tb_to3D( xb_ens,hxb_ens[:,:100] )
    #corr = mat_mult( xb_ens.transpose(), hxb_ens )
    #end = time.perf_counter()
    #print("Elapsed (after compilation) = {}s".format((end - start)))

    #f = open('./test_data', 'wb')
    #pickle.dump(corr,f)
    #f.close()

    #f = open('./test_data', 'rb')
    #test = pickle.load(f)
    #f.close()

    
    #corr = np.zeros(shape=[xmax*ymax,100])
    #for i in range(100):
    #    corr = np.matmul( np.transpose(xb_ens),hxb_ens[:,i] )
    
    #Lon_2d = Lon_all.reshape(xmax,ymax)
    #print(Lon_2d[0,:])
    #print(type(Lon_2d))
    #end_time = time.process_time()
    #print ('time needed: ', end_time-start_time, ' seconds')
    

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
import matlab.engine


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

    verify_hxb( stddev_hxb, wrf_dir )

    # May save the standard deviation
    if If_save:
        des_path = Hx_dir+ "Hxb_ensStddev_modelRes_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( stddev_hxb, f )
        f.close()
        print('Save '+des_path)

    return None


# Calculate values in observation space
def cal_pert_stddev_obsRes_Hxb( DAtime, sensor, Hx_dir, If_save, fort_v, wrf_dir=None):

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

    hxb_ens_os = np.zeros( shape=[num_ens+1,len(meanYb)] ) # the 0 to last - 1 rows are ensemble perturbations, the last row is mean value
    hxb_ens_os[:] = np.nan
    hxb_ens_os[num_ens,:] = meanYb
    
    # Read the ensemble of Hxb at obs locations
    ens_Tb_file = Hx_dir + "/Hxb_ens_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    print('Reading the ensemble of Hx: ', ens_Tb_file,'...')
    d_obspace = read_ens_obspace( ens_Tb_file, sensor )
    hxb_ens_os = d_obspace['hxb_ens']
    print('Shape of hxb_ens_obs: '+str(np.shape(hxb_ens_os)))
    for im in range( num_ens ):
        hxb_ens_os[im,:] = hxb_ens_os[im,:] - meanYb
    
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

    # Check if there are any NaN values using assert
    assert not np.isnan(stddev_hxb).any()

    # verify the standard deviation of hxb
    verify_hxb( stddev_hxb, wrf_dir )

    # May save the standard deviation
    if If_save:
        des_path = Hx_dir+ "Hxb_ensStddev_obsRes_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( stddev_hxb, f )
        f.close()
        print('Save '+des_path)

    return None


def cal_pert_stddev_xb( DAtime, wrf_dir, var_name, If_save ):
  
    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297
    if 'Q' in var_name:
        nLevel = 42
    elif 'W' in var_name:
        nLevel = 43

    ### ------------------------- Perturbation -------------------------
    print("Read the ensemble and calculate the perturbations for Xb...")
    start_time=time.process_time()
    
    xb_ens = np.zeros( shape=[nLevel,num_ens+1,xmax*ymax] )
    xb_ens[:] = np.nan
    # read the ensemble mean of xb
    mean_xb = wrf_dir + '/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    var = ncdir.variables[var_name][0,:,:,:]
    xb_ens[:,num_ens,:] = var.reshape(nLevel,xmax*ymax)
    # read the ensemble of xb
    file_xb = sorted( glob.glob(wrf_dir + '/wrf_enkf_input_d03_0*') )
    for ifile in file_xb:
        idx = file_xb.index( ifile )
        ncdir = nc.Dataset( ifile, 'r')
        var = ncdir.variables[var_name][0,:,:,:]
        xb_ens[:,idx,:] = var.reshape(nLevel,xmax*ymax) - xb_ens[:,num_ens,:]
    
    # Check if there are any NaN values using assert
    assert not np.isnan(xb_ens).any()

    # May save the perturbations
    if If_save:
        des_path = wrf_dir+ "xb_d03_3D_ensPert_" + DAtime + '_' + var_name + '.pickle'
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

    # May save the standard deviation
    if If_save:
        des_path = wrf_dir+ "xb_d03_3D_ensStddev_" + DAtime + '_' +  var_name + '.pickle'
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

    # Read ensemble perturbations of xb
    des_path = wrf_dir+ "xb_d03_3D_ensPert_" + DAtime + '_' + var_name +  '.pickle'
    with open( des_path,'rb' ) as f:
        xb_ens = pickle.load( f )
    print('Shape of xb_ens: '+ str(np.shape(xb_ens)))
    # Read ensemble standard deviation of xb
    des_path = wrf_dir+ "xb_d03_3D_ensStddev_" + DAtime + '_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_xb = pickle.load( f )
    print('Shape of stddev_xb: '+ str(np.shape(stddev_xb)))
    # Read ensemble perturbations of Hxb 
    if to_obs_res:
        des_path = Hx_dir+ "Hxb_ensPert_obsRes_" + DAtime + '_' +  sensor + '.pickle'
    else:
        des_path = Hx_dir+ "Hxb_ensPert_modelRes_" + DAtime + '_' +  sensor + '.pickle'
    with open( des_path,'rb' ) as f:
        hxb_ens = pickle.load( f )
    print('Shape of hxb_ens: '+ str(np.shape(hxb_ens)))
    # Read ensemble stand deviation of Hxb
    if to_obs_res:
        des_path = Hx_dir+ "Hxb_ensStddev_obsRes_" + DAtime + '_' +  sensor + '.pickle'
    else:
        des_path = Hx_dir+ "Hxb_ensStddev_modelRes_" + DAtime + '_' +  sensor + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_hxb = pickle.load( f )
    print('Shape of stddev_hxb: '+ str(np.shape(stddev_hxb)))

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
    if If_save:
        if to_obs_res:
            des_path = Hx_dir+ "Hcorr_Hxb_xb_obsRes_" + DAtime + '_' +  sensor + '_' + var_name +  '.pickle'
        else:
            des_path = Hx_dir+ "Hcorr_Hxb_xb_modelRes_" + DAtime + '_' +  sensor + '_' + var_name +  '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( corr_xb_hxb, f )
        f.close()
        print('Save '+des_path)
    return corr_xb_hxb


# ------------------------------------------------------------------------------------------------------
#           Operation: find obs of interest
# ------------------------------------------------------------------------------------------------------

# Collect indices that are cloudy for all ensemble members (not very useful for Hxb ens at obs location )
def find_cloudy_allEns( DAtime, sensor, Hx_dir):

    print("Looking for regions that are cloudy for all ensemble members...")
    start_time=time.process_time()

    ens_Tb_file = Hx_dir + "/Hxb_ens_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    print('Reading the ensemble of Hx: ', ens_Tb_file,'...')
    d_obspace = read_ens_obspace( ens_Tb_file, sensor )
    Tb_ens = d_obspace['hxb_ens']
    # Find the Tb threshold when all ens mems are cloudy
    cold_Tb_start = 80
    cold_Tb = cold_Tb_start
    common_idx = None
    while common_idx is None:
        idx_cloud_ens = np.where(Tb_ens<cold_Tb) # technique to save memory: sparisity
        # idea: any member should have at least one sample that meets the first requirement ( its value is less than the threshold)
        if not np.any( idx_cloud_ens ): # no member meets the 1st requirement
            common_idx = None
            cold_Tb = cold_Tb + 1
            continue
        elif len(np.unique(idx_cloud_ens[0])) != num_ens: # not all members meet the 1st requirement
            common_idx = None
            cold_Tb = cold_Tb + 1
            continue
        else: # all members meet the 1st requirement
            idx_mem1 = np.where(idx_cloud_ens[0]==0)[0]
            idx_mem2 = np.where(idx_cloud_ens[0]==1)[0]
            common_ini = list(set( idx_cloud_ens[1][idx_mem1] ).intersection( idx_cloud_ens[1][idx_mem2] ))
            for imem in range(1,num_ens-1):
                    if imem == 1:
                        common_idx = common_ini
                    idx_mem = np.where(idx_cloud_ens[0]==imem)[0]
                    common_idx = list(set( common_idx ).intersection( idx_cloud_ens[1][idx_mem] ))
            if not np.any( common_idx ): # not all members are cloudy at the same location(s)
                common_idx = None
                cold_Tb = cold_Tb + 1
                continue
            else:
                print('The Tb threshold when all members are cloudy at the same locations is: ' +str(cold_Tb)+ ' K!'  )
                print('Indices of cloudy region:'+str(common_idx))
                return common_idx


def find_coldest_obs( DAtime, sensor, Hx_dir ):

    print("Looking for cloudy regions based on observations...")
    idx_cloud = []
    # Read assimilated obs 
    file_Diag = big_dir+Storm+'/'+Exper_name+'/run/'+DAtime+'/enkf/d03/fort.10000'
    d_obs = Diag.Find_IR( file_Diag, fort_v )
    coldest_Tb = np.amin( d_obs['obs'] )
    print('The coldest Tb is '+str(coldest_Tb)+ ' K!')
    idx_cloud.append( d_obs['obs'].index( coldest_Tb ))
    print('Indices of cloudy region:'+str(idx_cloud))
    return idx_cloud,d_obs


# ------------------------------------------------------------------------------------------------------
#           Operation: Plot the correlation between the specified Tb and the 3D model (at specified levels)
# ------------------------------------------------------------------------------------------------------
def plot_corr_cloud_3Dmodel(big_dir, small_dir, Storm, Exper_name, var_name, DAtime, sensor, P_of_interest, d_model_res, cloud_idx ):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # Read Tbs of Hxb
    Tb_file = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/' + "/mean_obs_res_d03" + DAtime + '_' +  sensor + '.txt'
    d_all = ROIR.read_allTb(Tb_file, sensor )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = ROIR.read_TCvitals(small_dir+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(6.5,6.5), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    ### ---Plot Hxb---
    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )
    ax[0,0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax[0,0].coastlines(resolution='10m', color='black',linewidth=0.5)
    xb_Tb = ax[0,0].scatter(d_all['lon_obs'],d_all['lat_obs'],4,c=d_all['meanYb_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree()) 
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax[0,0].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())
    # Colorbar
    #caxes = fig.add_axes([0.12, 0.1, 0.45, 0.02])
    #Tb_bar = fig.colorbar(xb_Tb,ax=ax[0,0],orientation="horizontal", cax=caxes)
    #Tb_bar.ax.tick_params()

    ### ---Plot correlation---
    min_corr = -1
    max_corr = 1
    for isub in range(1,4):
        ax.flat[isub].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.flat[isub].coastlines(resolution='10m', color='black',linewidth=0.5)
        cs = ax.flat[isub].scatter(d_model_res['lon_model'],d_model_res['lat_model'],1.5,c=d_model_res['corr'][isub-1,:],\
                edgecolors='none', cmap='RdBu_r', vmin=min_corr, vmax=max_corr, transform=ccrs.PlateCarree())
        ax.flat[isub].scatter(d_all['lon_obs'][cloud_idx[0]],d_all['lat_obs'][cloud_idx[0]],2.5,c='black',edgecolors='none', transform=ccrs.PlateCarree())

    # Colorbar
    caxes = fig.add_axes([0.2, 0.05, 0.6, 0.02])
    #cb_corr_ticks = np.linspace(min_corr, max_corr, 5, endpoint=True)
    corr_bar = fig.colorbar(cs, ax=ax.flat[2:], orientation="horizontal", cax=caxes)
    corr_bar.ax.tick_params()

    #subplot title
    font = {'size':8,}
    ax.flat[0].set_title('H(xb)', font, fontweight='bold')
    for isub in range(1,4):
        ax.flat[isub].set_title( str(P_of_interest[isub-1])+' hPa', font, fontweight='bold')

    #title for all
    fig.suptitle(Storm+': '+Exper_name+'(corr of Qvapor&Tb)', fontsize=8, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(4):
        gl = ax.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top = False
        gl.xlabels_bottom = True
        if j==0 or j==2:
            gl.ylabels_left = True
            gl.ylabels_right = False
        else:
            gl.ylabels_left = False
            gl.ylabels_right = False

        if j==2 or j==3:
            gl.xlabels_bottom = True
            gl.xlabels_top = False
        else:
            gl.xlabels_bottom = False
            gl.xlabels_top = False


        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 8}
        gl.ylabel_style = {'size': 8}

    # Save figures
    des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Corr/IR/corr_'+DAtime+'_'+var_name+'_'+sensor+'.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)
    plt.close()
    return None


def corr_cloud_3Dmodel( DAtime, Hx_dir, sensor, wrf_dir, small_dir, Storm, Exper_name, var_name ):
 
    # Calculate the correlation
    idx_xb = np.arange(xmax*ymax)
    # Find cloudy region
    #cloud_hxb = find_cloudy_allEns( DAtime, sensor, Hx_dir)
    #cloud_hxb,d_obs = find_coldest_obs( DAtime, sensor, Hx_dir )
    cloud_hxb = [5000,]
    corr_xb_hxb = calculate_corr( DAtime, Hx_dir, sensor, wrf_dir, var_name, idx_xb, cloud_hxb)
    print('Shape of corr_xb_hxb: '+ str(np.shape(corr_xb_hxb)))

    # Average the correlation of the cloudy area
    corr_xb_hxb_cloud = np.mean(corr_xb_hxb, axis=2)
    print('Shape of corr_xb_hxb_cloud: '+ str(np.shape(corr_xb_hxb_cloud)))

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
    print('Locations for cloudy points:')
    for idx in cloud_hxb:
        print(str(lon[idx])+','+str(lat[idx]))

    # Specify pressure levels of interest
    P_of_interest = [100,500,850]
    
    # Calculate the corr at specified levels
    start_time=time.process_time()
    corr_xb_hxb_cloud_P = np.zeros( (len(P_of_interest),corr_xb_hxb.shape[1]),  )
    corr_xb_hxb_cloud_P[:] = np.nan
    for im in range( corr_xb_hxb.shape[1] ):
        f_interp = interpolate.interp1d( P_hpa[:,im], corr_xb_hxb_cloud[:,im])
        corr_xb_hxb_cloud_P[:,im] = f_interp( P_of_interest )
    #corr_xb_hxb_cloud_P = corr_xb_hxb_cloud[:3,:] # test....
    end_time = time.process_time()
    print ('time needed for the interpolation: ', end_time-start_time, ' seconds')
    print('Min of correlation: '+str(np.amin( corr_xb_hxb_cloud_P )))
    print('Max of correlation: '+str(np.amax( corr_xb_hxb_cloud_P )))

    d_model_res = {'lon_model':lon,'lat_model':lat,'corr':corr_xb_hxb_cloud_P}

    # Plot the correlation
    plot_corr_cloud_3Dmodel(big_dir, small_dir, Storm, Exper_name, var_name, DAtime, sensor, P_of_interest, d_model_res, cloud_hxb )



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    Exper_name = 'IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'

    v_interest = [ 'P',]
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709041600'
    end_time_str = '201709041600'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    to_obs_res = True
    ens_Interp_to_obs = False
    
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
            wrf_dir =  big_dir+Storm+'/'+Exper_name
            if to_obs_res:
                cal_pert_stddev_obsRes_Hxb( DAtime, sensor, Hx_dir, If_save, fort_v, wrf_dir)
            else:
                cal_pert_stddev_modelRes_Hxb( DAtime, sensor, Hx_dir, If_save, wrf_dir)
            # Xb
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in v_interest:
                cal_pert_stddev_xb( DAtime, wrf_dir, var_name, If_save )

    # Calculate correlations
    #if If_cal_corr:
    #    print('------------ Calculate the correlation between IR Tb and model variables --------------')
    #    for DAtime in IR_times:
    #        Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
    #        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
    #        for var_name in v_interest:
    #            calculate_corr( DAtime, Hx_dir, sensor, ch_list, wrf_dir, var_name, If_save )
   
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









    

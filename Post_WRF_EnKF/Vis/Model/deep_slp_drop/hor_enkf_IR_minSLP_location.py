
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

from wrf import getvar
from Track_xbxa import read_HPI_model
import Util_data as UD
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

## Calculate ensemble variance*(N-1) members for 1D variables
@njit(parallel=True)
def var_1D(A):
    # shape of A: num_ens
    res = np.zeros( (1), )
    for i in prange(A.shape[0]):
        res += A[i] * A[i]
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
    IR_obs = list(d_obs['obs'])
    lon_obs = list(d_obs['lon'])
    lat_obs = list(d_obs['lat'])
 
    hxb_ens = np.zeros( shape=[num_ens,len(d_obs['obs_type'])] )
    file_hxb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem0*.bin') )
    # Get common variables: lat,lon
    tmp_control = np.fromfile( file_hxb[0],dtype='<f4')
    n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
    tmp_data = tmp_control[:].reshape(n_ch+2,ymax,xmax)
    lon_x =  np.float64(list(tmp_data[0,:,:].flatten())) #longitude; added np.float64 bc eng.griddata returns error without it
    lat_x = np.float64(list(tmp_data[1,:,:].flatten()))  #latitude
    ## Interpolate to obs location
    # start a matlab process
    eng = matlab.engine.start_matlab() 
    for ifile in file_hxb:
        print(ifile)
        idx = file_hxb.index( ifile )
        tmp_control = np.fromfile( ifile,dtype='<f4') # <: little endian; f: float; 4: 4 bytes
        n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
        tmp_data = tmp_control[:].reshape(n_ch+2,ymax,xmax)
        Yb_x = np.float64(list(tmp_data[2,:,:].flatten()))
        
        for ich in ch_list:
            if ifile == file_hxb[0]:
                Ch_all = np.full( np.shape(IR_obs), ich )
            if sum(np.isnan(Yb_x)) != 0:
                warnings.warn('NaN value exists in Yb_x!') 
            # interpolate simulated Tbs to obs location
            mYb_obspace = eng.griddata(matlab.double(lon_x),matlab.double(lat_x),matlab.double(Yb_x),matlab.double(lon_obs),matlab.double(lat_obs))
            Yb_obspace = np.array(mYb_obspace._data)
            hxb_ens[idx,:] = np.around(Yb_obspace, decimals=4) #["{0:.4f}".format(item) for item in Yb_obspace]
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
    cbar.set_ticks( color_ticks )
    cbar.ax.tick_params(labelsize=12)

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(1):
        gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')

        gl.top_labels = False
        gl.bottom_labels = True
        if j==0:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False

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
def cal_pert_stddev_modelRes_Hxb( DAtime, sensor, Hx_dir, wrf_dir):

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

    #verify_hxb( stddev_hxb, wrf_dir )

    # May save the standard deviation
    if If_save:
        des_path = Hx_dir+ "Hxb_ensStddev_modelRes_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( stddev_hxb, f )
        f.close()
        print('Save '+des_path)

    return None

def cal_pert_stddev_obsRes_Hxb( DAtime, sensor, Hx_dir, fort_v, wrf_dir=None):

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
    #verify_hxb( d_obspace, DAtime, stddev_hxb,  wrf_dir )

    # May save the standard deviation
    if If_save:
        des_path = Hx_dir+ "Hxb_ensStddev_obsRes_" + DAtime + '_' +  sensor + '.pickle'
        f = open( des_path, 'wb' )
        pickle.dump( stddev_hxb, f )
        f.close()
        print('Save '+des_path)

    return None


def cal_pert_stddev_SLPatPoint( DAtime, wrf_dir, HPI_models ):
  
    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    ### ------------------------- Perturbation -------------------------
    print("Read the ensemble and calculate the perturbations at the location where min SLP of Xa is...")
     
    xb_ens = np.zeros( shape=[num_ens,] )
    xb_ens[:] = np.nan
    # grab the location where the min SLP of xa is
    idx_dt = IR_times.index( DAtime )
    idx_p = HPI_models['wrf_enkf_output_d03_mean']['idx_location'][idx_dt]
    # read the slp value at this location
    wrf_dir = wrf_dir+'/fc/'+DAtime
    mean_xb = wrf_dir+'/wrf_enkf_input_d03_mean'
    ncdir = nc.Dataset( mean_xb, 'r')
    slp = getvar(ncdir, 'slp')
    data_slp = slp.values
    mean_xb = data_slp.flatten()[idx_p] 
    # calculate the perturbation of slp between the ensemble member and the mean at this location
    file_xb = sorted( glob.glob(wrf_dir + '/wrf_enkf_input_d03_0*') )
    for ifile in file_xb:
        idf = file_xb.index( ifile )
        ncdir = nc.Dataset( ifile, 'r')
        slp = getvar(ncdir, 'slp')
        data_slp = slp.values
        xb_ens[idf] = data_slp.flatten()[idx_p] - mean_xb
        
    # Check if there are any NaN values using assert
    assert not np.isnan(xb_ens).any()

    ### ------------------------- Variance -------------------------
    print('Calculating the ensemble variance of model variable......' )
    start = time.perf_counter()
    var_xb = var_1D( xb_ens )
    end = time.perf_counter()
    print("Elapsed (after compilation) = {}s".format((end - start)))
    var_xb = var_xb / ( num_ens-1 )
    stddev_xb = np.sqrt( var_xb )

    return xb_ens,stddev_xb


# ------------------------------------------------------------------------------------------------------
#           Object: ensemble correlations of Xb and Hxb
# ------------------------------------------------------------------------------------------------------

# Calculate the ensemble correlation
def calculate_corr( DAtime, Hx_dir, sensor, wrf_dir, idx_hxb):

    # Read ensemble perturbations and standard deviations of xb
    xb_ens,stddev_xb = cal_pert_stddev_SLPatPoint( DAtime, wrf_dir, HPI_models ) 
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
    print('Calculating the covariance between Tbs and SLP at this point......' )
    start = time.perf_counter()
    xb_ens_2d = xb_ens.reshape(-1, 1)
    cov_xb_hxb = mat_mult( xb_ens_2d.transpose(),hxb_ens[:num_ens,idx_hxb] )
    end = time.perf_counter()
    print("Elapsed (after compilation) of covariance calculation = {}s".format((end - start)))
    cov_xb_hxb = cov_xb_hxb / ( num_ens-1 )

    # Calculate the correlation between Xb and Hxb
    print('Calculating the correlation between Tbs and SLP at this point......' )
    start = time.perf_counter()
    corr_xb_hxb = np.divide( cov_xb_hxb,stddev_xb )
    corr_xb_hxb = np.divide( corr_xb_hxb, stddev_hxb[idx_hxb] )
    end = time.perf_counter()
    print("Elapsed (after compilation) = {}s".format((end - start)))

    # sanity check
    assert  0 <= abs(corr_xb_hxb).all() and abs(corr_xb_hxb).all() <= 1
    
    return corr_xb_hxb

# Plot the correlation between the 2D Tb and the SLP at the specified location
def plot_hcorr_Tb_SLPatPoint(DAtime,sensor,corr_xb_hxb,d_all,idx_hxb,HPI_models):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, figsize=(6,6), dpi=400)
    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)

    min_corr = -0.5
    max_corr = 0.5
    if to_obs_res:
        if plot_scatter:
            cs = ax.scatter(d_all['lon_obs'][idx_hxb], d_all['lat_obs'][idx_hxb],8,c=corr_xb_hxb,\
                        edgecolors='none', cmap='bwr',vmin=min_corr,vmax=max_corr,transform=ccrs.PlateCarree())
        else:
            bounds = [-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5]#np.linspace(min_corr,max_corr,9)
            cs = ax.tricontourf(d_all['lon_obs'][idx_hxb], d_all['lat_obs'][idx_hxb], corr_xb_hxb, cmap='bwr', \
                        vmin=min_corr, vmax=max_corr, levels=bounds, transform=ccrs.PlateCarree(), extend='both' )   
    else:
        if plot_scatter:
            cs = ax.scatter(d_all['lon_model'][idx_hxb], d_all['lat_model'][idx_hxb],2,c=corr_xb_hxb,\
                        edgecolors='none', cmap='bwr',vmin=min_corr,vmax=max_corr,transform=ccrs.PlateCarree())
        else:
            bounds = [-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5]#np.linspace(min_corr,max_corr,9)
            cs = ax.tricontourf(d_all['lon_model'][idx_hxb], d_all['lat_model'][idx_hxb], corr_xb_hxb, cmap='bwr', \
                        vmin=min_corr, vmax=max_corr, levels=bounds, transform=ccrs.PlateCarree(), extend='both' )

    # Colorbar
    caxes = fig.add_axes([0.91, 0.1, 0.015, 0.8])
    cb_diff_ticks = [-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5]#np.linspace(min_corr, max_corr, 9, endpoint=True)
    cbar = fig.colorbar(cs, ticks=cb_diff_ticks,orientation="vertical",cax=caxes,extend='both')
    cbar.ax.tick_params(labelsize=12)

    # mark the location of min slp of simulation
    idx_t = IR_times.index( DAtime )
    ax.scatter(HPI_models['wrf_enkf_output_d03_mean']['lon'][idx_t], HPI_models['wrf_enkf_output_d03_mean']['lat'][idx_t], s=5, marker='o', facecolor='black',edgecolors='black', transform=ccrs.PlateCarree())

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(1):
        gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=1, color='gray', alpha=0.7, linestyle='--')

        gl.top_labels = False
        gl.bottom_labels = True
        if j==0:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False

        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 12}
        gl.ylabel_style = {'size': 12}

    #title
    title = 'Ens Hcorr of Tb and SLP at min SLP of Xa' 
    ax.set_title(title,fontsize=10, fontweight='bold')
    title_name = Storm+': '+Exper_name+' ('+DAtime+')'+'\nAbs of min SLP increment > '+\
                        str(incre_slp_th)+' hPa'+'\nCircle: center@min_slp of Xa, radius='+str(radius_th)+' km'
    fig.suptitle(title_name, fontsize=9, fontweight='bold')

    # Save figures
    if to_obs_res:
        des_name=plot_dir+DAtime+'_Hcorr_'+sensor+'_SLPatPoint_obsRes.png'
    else:
        des_name=plot_dir+DAtime+'_Hcorr_'+sensor+'_SLPatPoint_modelRes.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)
    plt.close()
    return None


# ------------------------------------------------------------------------------------------------------
#           Object: increment of Xb and Hxb
# ------------------------------------------------------------------------------------------------------

# Calculate the increment K(y-H(Xb))
def calculate_incre( DAtime, Hx_dir, sensor, wrf_dir, idx_hxb):

    # Read ensemble perturbations and standard deviations of xb
    xb_ens,stddev_xb = cal_pert_stddev_SLPatPoint( DAtime, wrf_dir, HPI_models )
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
    print('Calculating the covariance between Tbs and SLP at this point......' )
    start = time.perf_counter()
    xb_ens_2d = xb_ens.reshape(-1, 1)
    cov_xb_hxb = mat_mult( xb_ens_2d.transpose(),hxb_ens[:num_ens,idx_hxb] )
    end = time.perf_counter()
    print("Elapsed (after compilation) of covariance calculation = {}s".format((end - start)))
    cov_xb_hxb = cov_xb_hxb / ( num_ens-1 )

    # Calculate the K 
    sigma_o = 3
    sigma_hxb = stddev_hxb[idx_hxb]
    sigma_obs = sigma_hxb**2 + sigma_o**2
    K = np.divide(cov_xb_hxb,sigma_obs) 

    # Calculate the innovation y-H(Xb) at specified locations
    Tb_file = Hx_dir + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    d_all = ROIR.read_Tb_obsRes(Tb_file, sensor )
    y_hxb = d_all['Yo_obs'][idx_hxb]-d_all['meanYb_obs'][idx_hxb]
    
    # Calculate the increment 
    enkf_ic = K*y_hxb #elementwise multiplication

    return enkf_ic

# Plot the correlation between the 2D Tb and the SLP at the specified location
def plot_enkf_ic_Tb_SLPatPoint(DAtime,sensor,enkf_ic,d_all,idx_hxb,HPI_models):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, figsize=(6,6), dpi=400)
    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)

    min_ic = -10
    max_ic = 10
    cs = ax.scatter(d_all['lon_obs'][idx_hxb], d_all['lat_obs'][idx_hxb],8,c=enkf_ic,\
                edgecolors='none', cmap='bwr',vmin=min_ic,vmax=max_ic,transform=ccrs.PlateCarree())

    # Colorbar
    caxes = fig.add_axes([0.91, 0.1, 0.015, 0.8])
    cb_diff_ticks = np.linspace(min_ic, max_ic, 9, endpoint=True)
    cbar = fig.colorbar(cs, ticks=cb_diff_ticks,orientation="vertical",cax=caxes,extend='both')
    cbar.ax.tick_params(labelsize=12)

    # mark the location of min slp of simulation
    idx_t = IR_times.index( DAtime )
    ax.scatter(HPI_models['wrf_enkf_output_d03_mean']['lon'][idx_t], HPI_models['wrf_enkf_output_d03_mean']['lat'][idx_t], s=5, marker='o', facecolor='black',edgecolors='black', transform=ccrs.PlateCarree())

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(1):
        gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=1, color='gray', alpha=0.7, linestyle='--')

        gl.top_labels = False
        gl.bottom_labels = True
        if j==0:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False

        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 12}
        gl.ylabel_style = {'size': 12}

    #title
    title = 'EnKF increment of SLP at min SLP of Xa'
    ax.set_title(title,fontsize=10, fontweight='bold')
    title_name = Storm+': '+Exper_name+' ('+DAtime+')'+'\nAbs of min SLP increment > '+\
                        str(incre_slp_th)+' hPa'+'\nCircle: center@min_slp of Xa, radius='+str(radius_th)+' km'
    fig.suptitle(title_name, fontsize=9, fontweight='bold')

    # Save figures
    des_name=plot_dir+DAtime+'_EnKF_incre_'+sensor+'_SLPatPoint_obsRes.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)
    plt.close()
    return None

# ------------------------------------------------------------------------------------------------------
#           Object: y-H(Xb)
# ------------------------------------------------------------------------------------------------------
def plot_y_Hxb_SLPatPoint(DAtime,sensor,d_all,idx_hxb,HPI_models):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, figsize=(6,6), dpi=400)
    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)

    min_ic = -10
    max_ic = 10
    cs = ax.scatter(d_all['lon_obs'][idx_hxb], d_all['lat_obs'][idx_hxb],8,c=d_all['Yo_obs'][idx_hxb]-d_all['meanYb_obs'][idx_hxb],\
                edgecolors='none', cmap='bwr',vmin=min_ic,vmax=max_ic,transform=ccrs.PlateCarree())

    # Colorbar
    caxes = fig.add_axes([0.91, 0.1, 0.015, 0.8])
    cb_diff_ticks = [-10,-5,0,5,10] #np.linspace(min_ic, max_ic, 9, endpoint=True)
    cbar = fig.colorbar(cs, ticks=cb_diff_ticks,orientation="vertical",cax=caxes,extend='both')
    cbar.ax.tick_params(labelsize=12)

    # mark the location of min slp of simulation
    idx_t = IR_times.index( DAtime )
    ax.scatter(HPI_models['wrf_enkf_output_d03_mean']['lon'][idx_t], HPI_models['wrf_enkf_output_d03_mean']['lat'][idx_t], s=30, marker='o', facecolor='black',edgecolors='black', transform=ccrs.PlateCarree()) #s=5

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(1):
        gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=1, color='gray', alpha=0.7, linestyle='--')

        gl.top_labels = False
        gl.bottom_labels = True
        if j==0:
            gl.left_labels = True
            gl.right_labels = False
        else:
            gl.left_labels = False
            gl.right_labels = False

        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 12}
        gl.ylabel_style = {'size': 12}

    #title
    title = 'Innovation: Yo - H(Xb)'
    ax.set_title(title,fontsize=10, fontweight='bold')
    title_name = Storm+': '+Exper_name+' ('+DAtime+')'+'\nAbs of min SLP increment > '+\
                        str(incre_slp_th)+' hPa'+'\nCircle: center@min_slp of Xa, radius='+str(radius_th)+' km'
    fig.suptitle(title_name, fontsize=9, fontweight='bold')

    # Save figures
    des_name=plot_dir+DAtime+'_innovation_'+sensor+'_SLPatPoint_obsRes.png'
    plt.savefig( des_name, dpi=300)
    print('Saving the figure: ', des_name)
    plt.close()
    return None





if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    model_resolution = 3000 #m

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR'
    MP = 'WSM6'

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

    deep_slp_incre = True
    incre_slp_th = 10
    plot_circle = True
    radius_th = 200 # km

    to_obs_res = True
    ens_Interp_to_obs = False


    If_save = False
    If_plot = True
    plot_inno = True
    plot_corr = False
    plot_incre = False
    plot_scatter = True
    # -------------------------------------------------------    

    # Create experiment names
    Exper_name = UD.generate_one_name( Storm,DA,MP )
    #Exper_name = 'IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'

    if not Consecutive_times:
        IR_times = []
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Find the value of minimum slp increment
    if deep_slp_incre:
        # ----- Read min slp from model-----------------
        HPI_models = {}
        DAtimes_dir = [big_dir+Storm+'/'+Exper_name+'/fc/'+it for it in IR_times]
        file_kinds = ['wrf_enkf_input_d03_mean','wrf_enkf_output_d03_mean']
        for ifk in file_kinds:
            idx = file_kinds.index( ifk )
            HPI_models[ifk] = read_HPI_model( Storm, Exper_name, ifk, DAtimes_dir )
        incre_slp = np.array(HPI_models['wrf_enkf_output_d03_mean']['min_slp']) - np.array(HPI_models['wrf_enkf_input_d03_mean']['min_slp'])
    
    # Calculate the IR Tbs at obs locations for the whole ensemble
    if ens_Interp_to_obs:
        print('------- For all members, interpolate Hx in model resolution to obs location ------')
        for DAtime in IR_times:
            idx_t = IR_times.index( DAtime )
            if abs(incre_slp[idx_t]) > incre_slp_th:
                print('At '+DAtime)
                # Read assimilated obs 
                file_Diag = big_dir+Storm+'/'+Exper_name+'/run/'+DAtime+'/enkf/d03/fort.10000'
                d_obs = Diag.Find_IR( file_Diag, fort_v )
                # Interpolate
                Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
                interp_simu_to_obs_matlab_ens( d_obs, Hx_dir, sensor, ch_list, DAtime )


    # Calculate ensemble perturbations and variances
    if If_plot and deep_slp_incre:
        start_time=time.process_time()
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/deep_slp_incre/hcorr_Tb_SLPatPoint/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

        for DAtime in IR_times:
            idx_t = IR_times.index( DAtime )
            if abs(incre_slp[idx_t]) > incre_slp_th:
                print('At '+DAtime)
                # Calculate the EnKF update
                Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
                ct_lon = HPI_models['wrf_enkf_output_d03_mean']['lon'][idx_t]
                ct_lat = HPI_models['wrf_enkf_output_d03_mean']['lat'][idx_t]
                if to_obs_res:
                    Tb_file = Hx_dir + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
                    d_all = ROIR.read_Tb_obsRes(Tb_file, sensor )
                    idx_hxb = UD.find_circle_area_ungrid( ct_lon, ct_lat, d_all['lon_obs'], d_all['lat_obs'], radius_th )
                else:
                    Tb_file = Hx_dir + "/mean_model_res_d03_" + DAtime + '_' +  sensor + '.txt'
                    d_all = ROIR.read_Tb_modelRes(Tb_file, sensor )
                    idx_hxb = UD.find_circle_area_ungrid( ct_lon, ct_lat, d_all['lon_model'], d_all['lat_model'], radius_th )
                # Plot y-hxb
                if plot_inno:
                    plot_y_Hxb_SLPatPoint(DAtime,sensor,d_all,idx_hxb,HPI_models)

                if plot_corr or plot_incre:
                    # perturbation and stddev of Hxb
                    Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
                    wrf_dir =  big_dir+Storm+'/'+Exper_name
                    if to_obs_res:
                        cal_pert_stddev_obsRes_Hxb( DAtime, sensor, Hx_dir, fort_v, wrf_dir)
                    else:
                        cal_pert_stddev_modelRes_Hxb( DAtime, sensor, Hx_dir, wrf_dir)
                # Calculate the correlation between 2D Tb and the SLP at the min SLP of Xa
                if plot_corr:
                    corr_xb_hxb = calculate_corr( DAtime, Hx_dir, sensor, wrf_dir, idx_hxb)
                    corr_xb_hxb = np.squeeze( corr_xb_hxb )
                    print('Shape of corr_xb_hxb: '+ str(np.shape(corr_xb_hxb)))
                    plot_hcorr_Tb_SLPatPoint(DAtime,sensor,corr_xb_hxb,d_all,idx_hxb,HPI_models)
                # Calculate the enkf increment
                if plot_incre:
                    enkf_ic = calculate_incre( DAtime, Hx_dir, sensor, wrf_dir, idx_hxb)
                    enkf_ic = np.squeeze( enkf_ic )
                    print('Shape of enkf_ic: '+ str(np.shape(enkf_ic)))
                    plot_enkf_ic_Tb_SLPatPoint(DAtime,sensor,enkf_ic,d_all,idx_hxb,HPI_models)
            else:
                continue
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')


























    

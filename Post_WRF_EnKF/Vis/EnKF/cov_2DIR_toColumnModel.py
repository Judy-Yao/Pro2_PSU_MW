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
import warnings

import Util_data as UD
import Util_Vis
import Read_Obspace_IR as ROIR
import corr_2DIR_to3Dmodel as stat_3D
import corr_2DIR_toColumnModel as stat_Column
import Diagnostics as Diag
#import matlab.engine

# ------------------------------------------------------------------------------------------------------
#           Operation: Perform Operations
# ------------------------------------------------------------------------------------------------------

## Calculate the covariance*(N-1) between a Tb and a column of model variale over a 2D area (lon*lat)
@njit(parallel=True)
def cross_Tb_toCol( xb_ens,hxb_ens ):
    assert xb_ens.shape[1] == hxb_ens.shape[0]
    res = np.zeros( (xb_ens.shape[0],xb_ens.shape[2]),  )# levels, nobs
    # for each level of variable: calculate corr at each grid point
    for n in prange( xb_ens.shape[2] ): # loop thru samples
        for m in range( xb_ens.shape[1] ): # loop thru model levels
            for l in range( xb_ens.shape[0] ): # loop thru ens
                res[l,n] += xb_ens[l,m,n] * hxb_ens[m,n]
    return res

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


def plot_snapshot( lat,lon,Interp_cov,ver_coor ):


    # Read WRF domain
    wrf_file = wrf_dir+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )
    
    # Read Tbs 
    Tb_file = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/' + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    d_all = ROIR.read_Tb_obsRes(Tb_file, sensor )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = UD.read_TCvitals(small_dir+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(2, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(9.75,6.5), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    # Define the colorbar
    max_v = 0.005#0.000003 #max(np.amin(abs(Interp_incre)),np.amax(abs(Interp_incre)))
    min_v = -0.005#-0.000003#0-max_abs
    
    #bounds = [-0.00004,-0.00002,0,0.00002,0.00004]
    #bounds = [-0.00005,-0.00004,-0.00003,-0.00002,-0.00001,0,0.00001,0.00002,0.00003,0.00004,0.00005]

    for isub in range(6):
        ax.flat[isub].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.flat[isub].coastlines(resolution='10m', color='black',linewidth=0.5)
        #cs = ax.flat[isub].contourf(lon.reshape((297,297)),lat.reshape((297,297)),Interp_cov[isub,:].reshape((297,297)),cmap='RdBu_r',levels=bounds,extend='both',transform=ccrs.PlateCarree())
        #cs = ax.flat[isub].scatter(lon,lat,5,Interp_cov[isub,:],cmap='RdBu_r',edgecolors='none',transform=ccrs.PlateCarree(),)
        cs = ax.flat[isub].scatter(lon,lat,5,Interp_cov[isub,:],cmap='RdBu_r',edgecolors='none',vmin=min_v,vmax=max_v,transform=ccrs.PlateCarree(),)
        #cs = ax.flat[isub].scatter(lon,lat,5,Interp_cov[isub,:]*(d_all['Yo_obs']-d_all['meanYb_obs']),cmap='RdBu_r',edgecolors='none',transform=ccrs.PlateCarree(),)
        #cs = ax.flat[isub].scatter(lon,lat,5,Interp_cov[isub,:]*(d_all['Yo_obs']-d_all['meanYb_obs']),cmap='RdBu_r',vmin=min_v,vmax=max_v,edgecolors='none',transform=ccrs.PlateCarree(),)
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
            ax.flat[isub].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    color_ticks = np.linspace(min_v, max_v, 5, endpoint=True)
    cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, )
    cbar.set_ticks( color_ticks )
    cbar.ax.tick_params(labelsize=6)

    #subplot title
    font = {'size':15,}
    for isub in range(6):
        ax.flat[isub].set_title( str(ver_coor[isub])+' KM', font, fontweight='bold')

    #title for all
    fig.suptitle(Storm+':'+Exper_name+'~cov(Tb,'+var_name+')', fontsize=10, fontweight='bold')
    #fig.suptitle(Storm+':'+Exper_name+'~EnKF from Tb to Qvapor', fontsize=10, fontweight='bold')
    #fig.suptitle(Storm+':'+Exper_name+'~sd Qvapor', fontsize=10, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for j in range(6):
        gl = ax.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top = False
        gl.xlabels_bottom = True
        if j==0 or j==3:
            gl.ylabels_left = True
            gl.ylabels_right = False
        else:
            gl.ylabels_left = False
            gl.ylabels_right = False

        if j==3 or j==4 or j==5:
            gl.xlabels_bottom = True
            gl.xlabels_top = False
        else:
            gl.xlabels_bottom = False
            gl.xlabels_top = False

        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10}
        gl.ylabel_style = {'size': 12}

    # Save the figure
    if to_obs_res and interp_H:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/COV/Interp_H_cov_os_'+DAtime+'_'+var_name+'_'+sensor+'.png'
    else:
        save_des = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/COV/Interp_H_cov_ms_'+DAtime+'_'+var_name+'_'+sensor+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


# ------------------------------------------------------------------------------------------------------
#           Object: ensemble variances of columns of Xb and Hxb in 2D; Operation: Calculation
# ------------------------------------------------------------------------------------------------------
def cal_2Dcov_IR_ColVar( DAtime, var_name):

    if 'Q' in var_name:
        nLevel = 42

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

    # Find the location of model grid of interest
    if to_obs_res: # nearest for each obs
        idx_xb = stat_Column.Find_nearest_col( Hx_dir, wrf_dir, DAtime, sensor)
    else: # every model grid point
        idx_xb = np.arange(xmax*ymax)

    # Calculate the covariance between Xb and Hxb (a column of var and a Tb)
    print('Calculating the covariance between Tbs and ' + var_name + '(a column of var and a Tb)......' )
    start = time.perf_counter()
    cov_xb_hxb = cross_Tb_toCol( xb_ens[:,:num_ens,idx_xb],hxb_ens[:num_ens,:] )
    end = time.perf_counter()
    print("Elapsed (after compilation) of covariance calculation = {}s".format((end - start)))
    cov_xb_hxb = cov_xb_hxb / ( num_ens-1 )

    # calculate K
#    demo = stddev_hxb**2 + (3**2)
#    kgain = np.divide(cov_xb_hxb,demo)
#    cov_xb_hxb = kgain

    #tmp = np.divide(cov_xb_hxb,stddev_hxb[idx_xb])
    #tmp = np.divide(tmp,stddev_xb[:,idx_xb])
    #cov_xb_hxb = tmp

    #cov_xb_hxb = stddev_xb[:,idx_xb]

    # Make interpolation
    if interp_P:
        pass
    elif interp_H:
        mean_xa = wrf_dir + '/wrf_enkf_output_d03_mean'
        ncdir = nc.Dataset( mean_xa, 'r')
        #H_of_interest = [3.5,4.0,5.0,7.0,9.0,11.0]
        H_of_interest = [1.0,3.0,5.0,7.0,9.0,11.0]
        Interp_cov_xb_hxb = np.zeros( [len(H_of_interest),len(idx_xb)] )
        PHB = ncdir.variables['PHB'][0,:,:,:]
        PH = ncdir.variables['PH'][0,:,:,:]
        geoHkm = (PHB+PH)/9.8/1000 # in km
        geoHkm = geoHkm.reshape( geoHkm.shape[0],-1)
        geoHkm_half_eta = (geoHkm[:-1]+geoHkm[1:])/2
        geoHkm_half_eta = np.ma.getdata(geoHkm_half_eta)
        xlon = ncdir.variables['XLONG'][0,:,:].flatten()
        xlat = ncdir.variables['XLAT'][0,:,:].flatten()
        start_time=time.process_time()
        for im in range( len(idx_xb) ):
            f_interp = interpolate.interp1d( geoHkm_half_eta[:,im], cov_xb_hxb[:,im])
            Interp_cov_xb_hxb[:,im] = f_interp( H_of_interest )
        end_time = time.process_time()
        print ('time needed for the interpolation: ', end_time-start_time, ' seconds')

        if If_plot_snapshot:
            plot_snapshot( xlat[idx_xb],xlon[idx_xb],Interp_cov_xb_hxb,H_of_interest )
    else:
        pass
    


    return None



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'MARIA'
    Exper_name = 'IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900'

    v_interest = [ 'QVAPOR',]
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709160000'
    end_time_str = '201709160000'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    to_obs_res = True
    ens_Interp_to_obs = False

    interp_P = False
    P_of_interest = list(range( 995,49,-20 ))
    interp_H = True

    If_cal_pert_stddev = False
    If_cal_cov = True
    If_save = True

    If_plot_snapshot = True

   # -------------------------------------------------------  

    if not Consecutive_times:
        DAtimes = ['201709160600','201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Interpolate simulated Tb at model resolution to obs resolution
    if ens_Interp_to_obs:
        print('------- For all members, interpolate Hx in model resolution to obs location ------')
        for DAtime in DAtimes:
            # Read assimilated obs 
            file_Diag = big_dir+Storm+'/'+Exper_name+'/run/'+DAtime+'/enkf/d03/fort.10000'
            d_obs = Diag.Find_IR( file_Diag, fort_v )
            # Interpolate
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            stat_3D.interp_simu_to_obs_matlab_ens( d_obs, Hx_dir, sensor, ch_list,  DAtime )

    # Calculate ensemble perturbations and variances
    if If_cal_pert_stddev:
        print('------------ Calculate the ensemble perturbations --------------')
        for DAtime in DAtimes:
            # Hxb
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir =  big_dir+Storm+'/'+Exper_name
            if to_obs_res:
                print('At obs space...')
                stat_3D.cal_pert_stddev_obsRes_Hxb( DAtime, sensor, Hx_dir, If_save, fort_v, wrf_dir)
            else:
                print('At model space...')
                stat_3D.cal_pert_stddev_modelRes_Hxb( DAtime, sensor, Hx_dir, If_save, wrf_dir)
            # Xb
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in v_interest:
                stat_3D.cal_pert_stddev_xb( DAtime, wrf_dir, var_name, If_save )

    # Calculate correlations between obs and their nearest model columns
    if If_cal_cov:
        print('------------ Calculate the correlation between IR Tbs and columns of model variables --------------')
        for DAtime in DAtimes:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            for var_name in v_interest:
                print('Calculate '+var_name+'...')
                cal_2Dcov_IR_ColVar( DAtime, var_name)



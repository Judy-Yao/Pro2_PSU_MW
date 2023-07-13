#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

from numba import njit, prange
import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
from wrf import getvar, interplevel
import math
import scipy as sp
import scipy.ndimage
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time

import Util_data as UD

# ------------------------------------------------------------------------------------------------------
#           Operation: Perform Operations
# ------------------------------------------------------------------------------------------------------
@njit(parallel=True)
def hydro_mass( full_p, tv, geoHm, q ):
    R = 287.06 
    res = np.zeros( (q.shape[0],q.shape[1]), )
    res[:] = np.nan
    for im in prange( q.shape[1] ):
        for il in range( q.shape[0] ):
            zdiff = geoHm[il+1,im] - geoHm[il,im]   
            res[il,im] = (full_p[il,im]/(R * tv[il,im])) * q[il,im] * zdiff 
    # make sure all values are reasonable
    assert res.any() != np.nan
    return res


def compute_hydro( wrf_files, each_var):

    # Dimension of the domain
    xmax = 297
    ymax = 297
    nLevel = 42
    
    d_hydro = {}

    # read the shared variables: lat/lon
    ncdir = nc.Dataset( wrf_files[0] )
    lat = ncdir.variables['XLAT'][0,:,:].flatten()
    d_hydro['lat'] = lat
    lon = ncdir.variables['XLONG'][0,:,:].flatten()
    d_hydro['lon'] = lon
    # calculate fields needed for mass of hydrometeors
    full_p = np.zeros( [len(wrf_files), nLevel, xmax*ymax] ) 
    tv_k = np.zeros( [len(wrf_files), nLevel, xmax*ymax] )   
    geoHm = np.zeros( [len(wrf_files), nLevel+1, xmax*ymax] )  
    for wrf_file in wrf_files:
        ifile = wrf_files.index(wrf_file)
        ncdir = nc.Dataset( wrf_file, 'r')
        # full pressure
        p = ncdir.variables['P'][0,:,:,:] # perturbation
        pb = ncdir.variables['PB'][0,:,:,:]
        tmp_p = p + pb
        full_p[ifile,:,:] = tmp_p.reshape( tmp_p.shape[0],-1 ) 
        # geopotential height
        ph = ncdir.variables['PH'][0,:,:,:] # perturbation
        phb = ncdir.variables['PHB'][0,:,:,:]
        tmp_geoHm = (ph+phb)/9.8 # in meter
        geoHm[ifile,:,:] = tmp_geoHm.reshape( tmp_geoHm.shape[0],-1 )
        # temperature from potential temperature
        #theta = ncdir.variables['T'][0,:,:,:] # perturbation
        #full_theta = theta + 290 # 290: T_base set in namelist.input for wrf forecast
        # virtual temperature
        tv = getvar(ncdir,'tv',units='K')
        tmp_tv = tv.values
        tv_k[ifile,:,:] = tmp_tv.reshape( tmp_tv.shape[0],-1 ) 
    # calculate the mass of hydrometeor per layer
    for var_name in each_var:
        hydro_layer = np.zeros( [len(wrf_files), nLevel, xmax*ymax] )
        hydro_layer[:] = np.nan 
        for wrf_file in wrf_files:
            ifile = wrf_files.index(wrf_file)
            ncdir = nc.Dataset( wrf_file, 'r')
            ivar = each_var.index(var_name)
            var = ncdir.variables[var_name][0,:,:,:]
            var = var.reshape( var.shape[0],-1 )
            hydro_layer[ifile,:,:] = hydro_mass( full_p[ifile,:,:], tv_k[ifile,:,:], geoHm[ifile,:,:], var )
        d_hydro[var_name] = hydro_layer
 
    return d_hydro

def plot_IC(DAtime, plot_dir, d_hydro, var, IC_xb, IC_xa, ):

    # Read storm center
    dict_btk = UD.read_bestrack(Storm)
    # Find the best-track position
    btk_dt = [it_str for it_str in dict_btk['time'] ]#[datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]
    bool_match = [DAtime == it for it in btk_dt]
    if True in bool_match:
        if_btk_exist = True
        idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
    else:
        if_btk_exist = False

    #------ Domain -------------------
    lat = d_hydro['lat']
    lon = d_hydro['lon']
    lat_min = np.amin( lat )
    lon_min = np.amin( lon )
    lat_max = np.amax( lat )
    lon_max = np.amax( lon )
 
    fig, axs=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(12,5), dpi=400)

    if 'Q' in var:
        if var == 'QCLOUD':
            min_ichydro = 0
            max_ichydro = 3
        elif var == 'QGRAUP' or var == 'QRAIN':
            min_ichydro = 0
            max_ichydro = 10
        elif var == 'QICE':
            min_ichydro = 0
            max_ichydro = 0.5
        elif var == 'QSNOW':
            min_ichydro = 0
            max_ichydro = 5
        else:
            pass
    else:
        if var == 'liquid':
            min_ichydro = 0
            max_ichydro = 3
        elif var == 'ice':
            min_ichydro = 0
            max_ichydro = 3
        elif var == 'all_hydro':
            min_ichydro = 0
            max_ichydro = 3
        else:   
            pass

    # Xb
    axs.flat[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    axs.flat[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    #xb_hydro = axs.flat[0].scatter(lon,lat,1.5,c=IC_xb,edgecolors='none', cmap='ocean_r', transform=ccrs.PlateCarree())
    xb_hydro = axs.flat[0].scatter(lon,lat,1.5,c=IC_xb,edgecolors='none', cmap='ocean_r', vmin=min_ichydro, vmax=max_ichydro,transform=ccrs.PlateCarree())
    # Mark the best track
    if if_btk_exist:
        axs.flat[0].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'white', marker='*',transform=ccrs.PlateCarree())
    # Xa
    axs.flat[1].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    axs.flat[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    #xa_hydro = axs.flat[1].scatter(lon,lat,1.5,c=IC_xa,edgecolors='none', cmap='ocean_r',transform=ccrs.PlateCarree())
    xa_hydro = axs.flat[1].scatter(lon,lat,1.5,c=IC_xa,edgecolors='none', cmap='ocean_r', vmin=min_ichydro, vmax=max_ichydro,transform=ccrs.PlateCarree())
    # Mark the best track
    if if_btk_exist:
        axs.flat[1].scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 20, 'white', marker='*',transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.12, 0.1, 0.5, 0.02])
    ic_bar = fig.colorbar(xa_hydro,ax=axs[0:2],orientation="horizontal", cax=caxes,extend='max')
    ic_bar.ax.tick_params()

    # Xa-Xb (increment)
    if 'Q' in var:
        if var == 'QCLOUD':
            min_incre = -0.5
            max_incre = 0.5
        elif var == 'QGRAUP':
            min_incre = -1
            max_incre = 1
        elif var == 'QICE':
            min_incre = -0.05
            max_incre = 0.05
        elif var == 'QRAIN':
            min_incre = -0.5
            max_incre = 0.5
        elif var == 'QSNOW':
            min_incre = -1
            max_incre = 1
        else:
            pass
    else:
        if var == 'liquid':
            min_incre = -1
            max_incre = 1
        elif var == 'ice':
            min_incre = -1
            max_incre = 1
        elif var == 'all_hydro':
            min_incre = -1
            max_incre = 1
        else:
            pass


    axs.flat[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    axs.flat[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    #incre_hydro = axs.flat[2].scatter(lon,lat,1.5,c=IC_xa - IC_xb,edgecolors='none', cmap='bwr',transform=ccrs.PlateCarree())
    incre_hydro = axs.flat[2].scatter(lon,lat,1.5,c=IC_xa - IC_xb,edgecolors='none', cmap='bwr', vmin=min_incre, vmax=max_incre,transform=ccrs.PlateCarree())
    # Colorbar
    caxes = fig.add_axes([0.65, 0.1, 0.25, 0.02])
    cb_diff_ticks = np.linspace(min_incre, max_incre, 5, endpoint=True)
    #cbar = fig.colorbar(incre_hydro,ax=axs[2:],orientation="horizontal", cax=caxes)
    cbar = fig.colorbar(incre_hydro, ax=axs[2:], ticks=cb_diff_ticks, orientation="horizontal", cax=caxes, extend='both')
    cbar.ax.tick_params()

    #subplot title
    axs.flat[0].set_title('Xb--IC ($\mathregular{kgm^{-2}}$)',fontweight='bold')
    axs.flat[1].set_title('Xa--IC ($\mathregular{kgm^{-2}}$)',fontweight='bold')
    axs.flat[2].set_title('Xa-Xb--IC ($\mathregular{kgm^{-2}}$)', fontweight='bold')    #title for all
    fig.suptitle(Storm+': '+Exper_name+'('+var+' '+DAtime+')', fontsize=12, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2))
    for j in range(3):
        gl = axs.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
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
    
    des = plot_dir+DAtime+'_IC_'+var+'.png'
    plt.savefig( des, dpi=300 )
    print('Saving the figure: ', des)
    plt.close()







if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'HARVEY'
    Exper_name = 'JerryRun/IR_WSM6'
    each_var = ['QCLOUD','QRAIN','QICE','QSNOW','QGRAUP']

    # Time range set up
    start_time_str = '201708221200'
    end_time_str = '201708221200'
    Consecutive_times = True

    #Plot_layer = True
    each_water = False
    Plot_IC = True
    # -------------------------------------------------------    

    if not Consecutive_times:
        DAtimes = ['201709050000',]
        #DAtimes = ['201708230000','201708230600','201708231200','201708231800','201708240000','201708240600','201708241200']
        #DAtimes = ['201709031200','201709031800','201709040000','201709040600','201709041200','201709041800','201709050000']
        #DAtimes = ['201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Plot integrated column of hydrometeors
    if Plot_IC:
        start_time=time.process_time()
        for DAtime in DAtimes:
            print('At '+DAtime)
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'    
            wrf_files = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
            d_hydro = compute_hydro( wrf_files, each_var)
            # ------ Plot -------------------
            plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/Hydro_mass/'
            plotdir_exists = os.path.exists( plot_dir )
            if plotdir_exists == False:
                os.mkdir(plot_dir)
            # --- Condition
            if each_water:
                v_interest = each_var
                for var in v_interest:
                    print('Plotting IC of '+var)
                    IC_xb = np.sum(d_hydro[var][0,:,:],axis=0)
                    IC_xa = np.sum(d_hydro[var][1,:,:],axis=0)   
                    plot_IC(DAtime, plot_dir, d_hydro, var, IC_xb, IC_xa, )       
            else:
                v_interest = ['liquid','ice','all_hydro']
                for var in v_interest:
                    if var == 'liquid':
                        tmp = d_hydro['QCLOUD']+d_hydro['QRAIN']
                    elif var == 'ice':
                        tmp = d_hydro['QICE']+d_hydro['QSNOW']+d_hydro['QGRAUP']
                    elif var == 'all_hydro':
                        tmp = d_hydro['QCLOUD']+d_hydro['QRAIN']+d_hydro['QICE']+d_hydro['QSNOW']+d_hydro['QGRAUP']
                    IC_xb = np.sum( tmp[0,:,:],axis=0 )
                    IC_xa = np.sum( tmp[1,:,:],axis=0 )
                    plot_IC(DAtime, plot_dir, d_hydro, var, IC_xb, IC_xa, )

        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')















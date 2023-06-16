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
from matplotlib import ticker
from matplotlib import colors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import pickle

import Read_Obspace_IR as ROIR

def fmt(x,pos):
    a,b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a,b)



def snapshot_H( xlat,xlon,Interp_diff_xa ):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name[1]+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = ROIR.read_TCvitals(small_dir+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(2, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(9.75,6.5), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']  

    # Define the colorbar
    max_var = 0.001 #max(np.amin(abs(Interp_incre)),np.amax(abs(Interp_incre)))
    min_var = -0.001#0-max_abs

    for isub in range(6):
        ax.flat[isub].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.flat[isub].coastlines(resolution='10m', color='black',linewidth=0.5)
        cs = ax.flat[isub].scatter(xlon,xlat,1.5,c=Interp_diff_xa[isub,:],\
                edgecolors='none', cmap='RdBu_r', vmin=min_var, vmax=max_var, transform=ccrs.PlateCarree())
                #edgecolors='none', cmap='RdBu_r',transform=ccrs.PlateCarree())
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
            ax.flat[isub].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.8])
    #cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, format=ticker.FuncFormatter(fmt))
    cbar = fig.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, )
    cbar.ax.tick_params(labelsize=8)

    #subplot title
    font = {'size':15,}
    for isub in range(6):
        ax.flat[isub].set_title( str(H_of_interest[isub])+' KM', font, fontweight='bold')

    #title for all
    fig.suptitle(Storm+'at 1st cycle: (IR_Xa-conv_Xa)~'+MPS+' '+var_name, fontsize=10, fontweight='bold')

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
    save_des = small_dir+Storm+'/'+Exper_name[0]+'/Vis_analyze/EnKF/ExtraIncreFromIR_'+var_name+'_'+DAtime+'.png'
    plt.savefig( save_des )
    print('Saving '+ save_des)



def extraIncre_from_IR():

    if 'Q' in var_name:
        nLevel = 42

    # Read enkf output from conv-only experiment 
    conv_file = big_dir+Storm+'/'+Exper_name[1]+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    ncdir = nc.Dataset( conv_file, 'r')
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()
    nc_conv = nc.Dataset( conv_file, 'r')
    conv_xa = nc_conv.variables[var_name][0,:,:,:]
    # Read enkf output from IR experiment
    IR_file = big_dir+Storm+'/'+Exper_name[0]+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    nc_IR = nc.Dataset( IR_file, 'r')
    IR_xa = nc_IR.variables[var_name][0,:,:,:]
    # Calculate the difference (IR - conv)
    diff_xa = IR_xa - conv_xa
    diff_xa = diff_xa.reshape( diff_xa.shape[0],-1)

    # Make interpolation
    if interp_P:
        Interp_diff_xa = np.zeros( [len(P_of_interest),xmax*ymax] )
        # ---------- Interpolate to specified pressure levels ----------
        PB = nc_conv.variables['PB'][0,:,:,:]
        P = nc_conv.variables['P'][0,:,:,:]
        P_hpa = (PB + P)/100
        P_hpa = P_hpa.reshape( P_hpa.shape[0],-1)
        start_time=time.process_time()
        for im in range( P_hpa.shape[1] ):
            f_interp = interpolate.interp1d( P_hpa[:,im], diff_xa[:,im])
            Interp_diff_xa[:,im] = f_interp( P_of_interest )
        end_time = time.process_time()
        print ('time needed for the interpolation: ', end_time-start_time, ' seconds')
        if If_plot:
            snapshot_P( xlat,xlon,Interp_diff_xa )

    elif interp_H:
        # ---------- Interpolate to specified geopotential heights ----------
        Interp_diff_xa = np.zeros( [len(H_of_interest),xmax*ymax] )
        PHB = nc_conv.variables['PHB'][0,:,:,:]
        PH = nc_conv.variables['PH'][0,:,:,:]
        geoHkm = (PHB+PH)/9.8/1000 # in km
        geoHkm = geoHkm.reshape( geoHkm.shape[0],-1)
        geoHkm_half_eta = (geoHkm[:-1]+geoHkm[1:])/2
        geoHkm_half_eta = np.ma.getdata(geoHkm_half_eta)
        start_time=time.process_time()
        for im in range( geoHkm_half_eta.shape[1] ):
            f_interp = interpolate.interp1d( geoHkm_half_eta[:,im], diff_xa[:,im])
            Interp_diff_xa[:,im] = f_interp( H_of_interest )
        end_time = time.process_time()
        print ('time needed for the interpolation: ', end_time-start_time, ' seconds')
        if If_plot:
            snapshot_H( xlat,xlon,Interp_diff_xa )
    else:
        pass
        
    return None



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'MARIA'
    Exper_name = ['IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900','J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900',]
    v_interest = [ 'QVAPOR',]
    DAtime = '201709160000'
    MPS = 'WSM6'

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    interp_P = False
    P_of_interest = list(range( 995,49,-20 ))
    interp_H = True
    H_of_interest = [1.0,3.0,5.0,7.0,9.0,11.0,]

    If_plot_snapshot = False
    If_plot = True
    # -------------------------------------------------------   
    for var_name in v_interest:
        print('Plot '+var_name+'...')
        extraIncre_from_IR()
    
# Note: the pressure and geopotential height difference of wrf_enkf_output_d03_mean between IR and conv-only is small with less than 1 hPa and 20 meters.

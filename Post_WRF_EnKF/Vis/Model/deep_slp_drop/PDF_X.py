#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

from numba import njit, prange
import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
from wrf import getvar, interplevel
from scipy import interpolate
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
import warnings

import Util_Vis
import Util_data as UD
import vertical_coor as VC
import Read_Obspace_IR as ROIR
import Diagnostics as Diag

# ------------------------------------------------------------------------------------------------------
#           Operation: Perform Numba Operations
# ------------------------------------------------------------------------------------------------------
#@njit(parallel=True)
def hist2d_numba_seq(tracks, bins, ranges):
    H = np.zeros( (bins[0],bins[1]),dtype=np.uint64 )
    delta = 1 / ((ranges[:, 1] - ranges[:, 0]) / bins)

    for t in range(tracks.shape[1]):
        i = (tracks[0, t] - ranges[0, 0]) * delta[0]
        j = (tracks[1, t] - ranges[1, 0]) * delta[1]
        if 0 <= i < bins[0] and 0 <= j < bins[1]:
            H[int(i), int(j)] += 1

    return H

# Calculate hydro mass at any grid point
@njit(parallel=True)
def hydro_mass( full_p, tv, geoHm, q ):
    R = 287.06
    res = np.zeros( (q.shape[0],q.shape[1]), )
    res[:] = np.nan
    for im in prange( q.shape[1] ):
        for il in range( q.shape[0] ):
            zdiff = geoHm[il+1,im] - geoHm[il,im]
            res[il,im] = (full_p[il,im]/(R * tv[il,im])) * q[il,im] * zdiff #* model_resolution**2
    # make sure all values are reasonable
    assert res.any() != np.nan
    return res

# ------------------------------------------------------------------------------------------------------
#           Operation: Read and process data
# ------------------------------------------------------------------------------------------------------

# Read fields needed for mass of hydrometeors
def Read_for_hydromass( DAtimes ):

    d_forHydro = {}
    for DAtime in DAtimes:
        d_forHydro[DAtime] = {}

    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
        wrf_files = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
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
            tmp_geoHm = (ph+phb)/9.81 # in meter
            geoHm[ifile,:,:] = tmp_geoHm.reshape( tmp_geoHm.shape[0],-1 )
            # virtual temperature
            tv = getvar(ncdir,'tv',units='K')
            tmp_tv = tv.values
            tv_k[ifile,:,:] = tmp_tv.reshape( tmp_tv.shape[0],-1 )

        # Assemble a dictionary
        d_forHydro[DAtime]['full_p'] = full_p
        d_forHydro[DAtime]['geoHm'] = geoHm
        d_forHydro[DAtime]['tv_k'] = tv_k

    return d_forHydro

# Process data 
# 1.interpolate model data to vertical coordinate of interest
# 2.collect over time
def read_modelatObsRes( var_name,DAtimes,d_forHydro=None ):

    d_allx = {}
    for DAtime in DAtimes:
        d_allx[DAtime] = {}

    for DAtime in DAtimes:
        print('Reading '+DAtime+'...')
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
        wrf_files = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
        
        # Find the location of model grid of interest
        if to_obs_res: # nearest for each obs
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            idx_x = VC.Find_nearest_col( Hx_dir, wrf_dir, DAtime, sensor)
        else:
            idx_x = range(xmax*ymax) 

        # impose other conditions
        if limit:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            Tb_file = Hx_dir + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
            d_all = ROIR.read_Tb_obsRes(Tb_file, sensor )
            condi1 = d_all['meanYb_obs'] - d_all['Yo_obs'] > 3 #d_all['Yo_obs'] <= 210
            idx_a = np.where( condi1 )[0]
            condi2 = d_all['Yo_obs'] <= 220 #d_all['meanYb_obs'] <= 210
            idx_b = np.where( condi2 )[0]
            idx_xx = list(set(idx_a)&set(idx_b))

        # Collect data
        for wrf_file in wrf_files:
            ncdir = nc.Dataset( wrf_file, 'r')
            # Read variable of interest
            var_tmp = ncdir.variables[var_name][0,:,:,:]
            if var_name == 'W':
                var_tmp = var_tmp.reshape(nLevel+1,xmax*ymax)
                var_half = (var_tmp[:-1]+var_tmp[1:])
                var_tmp = np.ma.getdata( var_half )
                var = var_tmp[:, idx_x]
                if limit:
                    var = var[:,idx_xx]
            elif 'Q' in var_name:
                var_tmp = ncdir.variables[var_name][0,:,:,:]
                var_tmp = var_tmp.reshape( var_tmp.shape[0],-1 )
                # calculate the mass of hydrometeor at any given point for a file
                ifile = wrf_files.index(wrf_file)
                var = hydro_mass( d_forHydro[DAtime]['full_p'][ifile,:,idx_x], d_forHydro[DAtime]['tv_k'][ifile,:,idx_x], d_forHydro[DAtime]['geoHm'][ifile,:,idx_x], var_tmp[:,idx_x] )*1000 # to gram/m2
                if limit:
                    var = var[:,idx_xx]
            else:
                var_tmp = var_tmp.reshape(nLevel,xmax*ymax)
                var = var_tmp[:, idx_x]
                if limit:
                    var = var[:,idx_xx]

            # Read vertical coordinates
            if interp_P and not interp_H:
                p = ncdir.variables['P'][0,:,:,:] # perturbation
                pb = ncdir.variables['PB'][0,:,:,:]
                tmp_p = p + pb
                tmp_p = tmp_p.reshape( tmp_p.shape[0],-1 )
                ver_coor = tmp_p[:,idx_x]
                if limit:
                    ver_coor = ver_coor[:,idx_xx]
                Level_of_interest = P_like 
            elif not interp_P and interp_H:
                ph = ncdir.variables['PH'][0,:,:,:] # perturbation
                phb = ncdir.variables['PHB'][0,:,:,:]
                tmp_geoHm = (ph+phb)/9.81 # in meter
                tmp_geoHkm = tmp_geoHm.reshape( tmp_geoHm.shape[0],-1 )/1000
                geoHkm_half_eta = (tmp_geoHkm[:-1,:]+tmp_geoHkm[1:,:])/2 # nLevel+1 to nLevel
                geoHkm_half_eta = np.ma.getdata(geoHkm_half_eta)
                ver_coor = geoHkm_half_eta[:,idx_x]
                if limit:
                    ver_coor = ver_coor[:,idx_xx]
                Level_of_interest = H_like

            # Interpolate data to vertical coordinate of interest
            var_interp = np.zeros( [ len(Level_of_interest),np.shape(ver_coor)[1]] )    
            var_interp[:] = np.nan 
            # Interpolate to vertical level of interest
            start_time=time.process_time()
            for im in range( ver_coor.shape[1] ): # number of points
                f_interp = interpolate.interp1d( ver_coor[:,im], var[:,im] )
                var_interp[:,im] = f_interp( Level_of_interest )
            end_time = time.process_time()
            print ('time needed for the interpolation: ', end_time-start_time, ' seconds')

            # Assemble a dictionary
            if 'input' in wrf_file:
                d_allx[DAtime]['input'] = var_interp
            elif 'output' in wrf_file:
                d_allx[DAtime]['output'] = var_interp
            else:
                pass

    return d_allx

# Plot the density of histogram
def plot_hist( var,key,H_all ):
    
    # Use of imshow. For example,
    # a is a 2-by-2 array
    # b = plt.imshow(a)
    # Values of x AND y axes will be -0.5,0,0.5,1,1.5

    # ------ Plot Figure -------------------
    fig, ax=plt.subplots(1, 2,figsize=(12,6), dpi=400)
    
    ## x axis: value range of variable
    #if 'Q' in var:
    #    x_range = list(np.logspace(0,1,num=100))
    #    x_range.insert(0,0)
    #else:
    #    x_range = var_rg
    x_range = var_rg
    x_axis_rg = range(len(x_range ))
    f_xinterp = interpolate.interp1d( x_range, x_axis_rg)
    ## y axis: vertical coordinate
    if interp_P:
        pass
    else:
        st = 0
        end = 21
        step = (end-st)/number_bins
        y_range = np.arange(st,end,step)
    y_axis_rg = range(len(y_range))
    f_yinterp = interpolate.interp1d( y_range, y_axis_rg)

    # Set colorbar
    if key == 'AllTimes':
        vmax = 1e7
    elif key != 'AllTimes' and not limit:
        vmax = 1e4
    elif key is not 'AllTimes' and limit:
        vmax = 1e3

    # Plot
    for i in range(2):
        x = x_axis_rg 
        y = y_axis_rg
        X,Y=np.meshgrid(x,y)

        im = ax[i].pcolor(X,Y,H_all[key][i,:,:], cmap='rainbow_r', vmax=vmax, norm=mcolors.LogNorm())
#        im = ax[i].imshow(H_all[key][i,:,:],cmap='rainbow_r',extent=(0,len(x_range),len(y_range),0),vmax=vmax,norm=mcolors.LogNorm())
#        ax[i].invert_yaxis()

        # Set the x-ticks
        if var == 'W':
            xlabel_like = [-1.5,-1,-0.5,0,0.5,1,1.5,2,0,2.5,3.0]
        elif 'Q' in var:
            xlabel_like = [0,1,2,3,4,5]
        xticks = []
        for it in xlabel_like:
            xticks.append( f_xinterp( it ) )
        ax[i].set_xticks( xticks )
        ax[i].set_xticklabels( [str(it) for it in xlabel_like],fontsize=15 )
        if var == 'W':
            ax[i].set_xlabel(var+' (m s-1)',fontsize=20)
        elif 'Q' in var:
            ax[i].set_xlabel(var+' (gram m-2)',fontsize=20)

        # Plot a line indicating the threshold value
        if var == 'W':
            th_loc = f_xinterp(0)
            ax[i].axvline(x=th_loc,color='black',linestyle='-',linewidth=2)

        # Set the y-ticks
        if interp_P:
            pass
        else:
            ylabel_like = [0.5,5.0,10.0,15.0,20.0]
        yticks = []
        for it in ylabel_like:
            yticks.append( f_yinterp( it ) )
        ax[i].set_yticks( yticks )
        ax[i].set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
        if interp_P:
            pass
        else:
            ax[0].set_ylabel('Height (km)',fontsize=20)
        #ax[i].set_ylim(ymin=f_yinterp(0),ymax=f_yinterp(20)) 


    # Add the color bar
    cbaxes = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    cbar = fig.colorbar(im,cax=cbaxes,fraction=0.046, pad=0.04,extend='max')
    cbar.ax.tick_params(labelsize=15)

    # title
    ax[0].set_title( 'Xb', fontsize = 15, fontweight='bold')
    ax[1].set_title( 'Xa', fontsize = 15, fontweight='bold')
    if to_obs_res and limit:
        title_name = Storm+': '+Exper_name+' '+key+'\n @Obs_location: H(Xb)-obs > 3K (obs<=220K)'
    if to_obs_res and not limit:
        title_name = Storm+': '+Exper_name+' '+key+'\n @Obs_location'
    fig.suptitle(title_name, fontsize=12, fontweight='bold')

    if to_obs_res and limit:
        save_des = plot_dir+key+'_'+var+'overH_PDF_obsRes_limit.png'
    if to_obs_res and not limit:
        save_des = plot_dir+key+'_'+var+'overH_PDF_obsRes.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()




def set_range( var ):

    # range of var values
    if var == 'W':
        min_var_rg = -1.5
        max_var_rg = 3.5
    elif 'Q' in var:
        min_var_rg = 0
        max_var_rg = 5

    return min_var_rg, max_var_rg

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR'
    MP = 'WSM6'

    v_interest = [ 'W','QICE','QCLOUD','QRAIN','QSNOW','QGRAUP','QVAPOR',]
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709030000'
    end_time_str = '201709031200'
    Consecutive_times = True

    # Number of ensemble members
    num_ens = 50
    # Dimension of the domain
    xmax = 297
    ymax = 297
    nLevel = 42 

    # field on model/obs locations 
    to_obs_res = True 

    # number of bins
    number_bins = 100
    # vertical coordinate of interest
    interp_P = False
    P_like = np.linspace(900,10,number_bins)
    interp_H = True 
    H_like = np.linspace(0.5,21,number_bins)
    # limitations
    limit = True

    if_plot = True

    # ------------------------------------------------------- 

    Exper_name = UD.generate_one_name( Storm,DA,MP )#'IR-updateW-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'#UD.generate_one_name( Storm,DA,MP )

    if not Consecutive_times:
        DAtimes = ['201709030000','201709030600','201709031200','201709031800','201709040000','201709040600','201709041200','201709041800','201709050000']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Make plot dir
    plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/PDF_X/'
    plotdir_exists = os.path.exists( plot_dir )
    if plotdir_exists == False:
        os.mkdir(plot_dir)

    # Calculate fields needed for mass of hydrometeors
    d_forHydro = None
    for var in v_interest:
        if 'Q' in var:
             d_forHydro = Read_for_hydromass( DAtimes )
             break
 
    # Loop
    for var in v_interest:
        # Make plot dir
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/PDF_'+var+'/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
        # Read and process model data
        d_allx = read_modelatObsRes( var, DAtimes, d_forHydro)
        # ---- Make bin counts -----
        if interp_P and not interp_H:
            Level_of_interest = P_like
        elif not interp_P and interp_H:
            Level_of_interest = H_like
        # set up bins and ranges 
        min_var_rg, max_var_rg = set_range( var )
        var_rg = np.linspace(min_var_rg, max_var_rg,number_bins+1)

        # counts 
        H_all = {} # np.zeros((2,bins[0], bins[1])) #input/output,num_levels,num_var_values
        
        H_allT = np.zeros((2,len(Level_of_interest), number_bins))
        for outer_key in d_allx:
            print(outer_key)
            H_one = np.zeros((2,len(Level_of_interest), number_bins))
            for inner_key in d_allx[outer_key]:
                    for il in range(len(Level_of_interest)):
                        tmp = np.histogram(d_allx[outer_key][inner_key][il,:], number_bins, (min_var_rg,max_var_rg))
                        if inner_key == 'input':
                            H_one[0,il,:] = tmp[0] 
                        elif inner_key == 'output':
                            H_one[1,il,:] = tmp[0]
            H_all[outer_key] = H_one
            H_allT = H_allT+H_one
        #H_all['AllTimes'] = H_allT
        

        # plot
        if if_plot:
            for key in H_all:
                plot_hist( var,key,H_all )



















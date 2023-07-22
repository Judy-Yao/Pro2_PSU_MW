#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

from numba import njit, prange
import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
from wrf import getvar, ll_to_xy
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
from scipy import interpolate

from Track_xbxa import read_HPI_model
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
            res[il,im] = (full_p[il,im]/(R * tv[il,im])) * q[il,im] * zdiff #* model_resolution**2
    # make sure all values are reasonable
    assert res.any() != np.nan
    return res



# Return: for each var, a 3d field with dimension: number of files, 
# & number of levels, and number of points in the area
def compute_hydro_area( wrf_files, hydros, HPI_models, idx_t):

    # Dimension of the domain
    nLevel = 42

    d_hydro = {}
    
    # Find the min slp in enkf input as the anchor point
    # Find the circled area with the min slp as the center
    output_file = wrf_files[1]
    ncdir = nc.Dataset( output_file, 'r')
    anchor_ij = ll_to_xy(ncdir, HPI_models['wrf_enkf_output_d03_mean']['lat'][idx_t], HPI_models['wrf_enkf_output_d03_mean']['lon'][idx_t])
    idx_area = UD.find_circle_area( output_file, anchor_ij.values[0], anchor_ij.values[1], radius_th, model_resolution/1000) # in km
    #print(len(idx_area))

    # calculate fields needed for mass of hydrometeors
    full_p = np.zeros( [len(wrf_files), nLevel, len(idx_area)] )
    tv_k = np.zeros( [len(wrf_files), nLevel, len(idx_area)] )
    geoHm = np.zeros( [len(wrf_files), nLevel+1, len(idx_area)] ) 
    for wrf_file in wrf_files:
        ifile = wrf_files.index(wrf_file)
        ncdir = nc.Dataset( wrf_file, 'r')
        # full pressure
        p = ncdir.variables['P'][0,:,:,:] # perturbation
        pb = ncdir.variables['PB'][0,:,:,:]
        tmp_p = p + pb
        tmp_p = tmp_p.reshape( tmp_p.shape[0],-1 )
        full_p[ifile,:,:] = tmp_p[:,idx_area]
        # geopotential height
        ph = ncdir.variables['PH'][0,:,:,:] # perturbation
        phb = ncdir.variables['PHB'][0,:,:,:]
        tmp_geoHm = (ph+phb)/9.81 # in meter
        tmp_geoHm = tmp_geoHm.reshape( tmp_geoHm.shape[0],-1 )
        geoHm[ifile,:,:] = tmp_geoHm[:,idx_area]
        # virtual temperature
        tv = getvar(ncdir,'tv',units='K')
        tmp_tv = tv.values
        tmp_tv = tmp_tv.reshape( tmp_tv.shape[0],-1 )
        tv_k[ifile,:,:] = tmp_tv[:,idx_area]
    # calculate the mass of hydrometeor per layer
    for var_name in hydros:
        hydro_layer = np.zeros( [len(wrf_files), nLevel, len(idx_area)] )
        hydro_layer[:] = np.nan
        for wrf_file in wrf_files:
            ifile = wrf_files.index(wrf_file)
            ncdir = nc.Dataset( wrf_file, 'r')
            ivar = hydros.index(var_name)
            var = ncdir.variables[var_name][0,:,:,:]
            var = var.reshape( var.shape[0],-1 )
            hydro_layer[ifile,:,:] = hydro_mass( full_p[ifile,:,:], tv_k[ifile,:,:], geoHm[ifile,:,:], var[:,idx_area] )
        d_hydro[var_name] = hydro_layer

    # Calculate the vertical coordinate
    Height = np.mean( geoHm, axis=0) 
    Height = np.mean( Height, axis=1 )/1000 # convert to km
    Height = (Height[:-1]+Height[1:])/2
    return d_hydro, Height


def plot_all_multiTime(DAtime, plot_dir, d_hydro, Height ):

    fig, ax=plt.subplots(5, 1, sharex='all', figsize=(10,12), dpi=300)

    # Manually set discrete values on x and y axis and interpolate data to these values
    ## x axis: correlation value
    x_range = np.arange(0,3000.5,50) #2100.5
    x_axis_rg = range(len(x_range))
    f_xinterp = interpolate.interp1d( x_range, x_axis_rg)
    ## y axis: model level height
    y_range = np.arange(0,31,1)
    y_axis_rg = range(len(y_range))
    f_yinterp = interpolate.interp1d( y_range, y_axis_rg)
    loc_iny = f_yinterp( Height )

    labels = ['Xb_','Xa_']
    lstyle = ['-','--']
    Color = ['#f032e6','#911eb4','#4363d8','#f58231','#469990']

    for key in d_hydro.keys():
        for ifile in range( len(wrf_files) ):
            for var in hydros:
                idx = hydros.index( var )
                PF_x = np.sum(d_hydro[key][var][ifile,:,:],axis=1)
                loc_inx = f_xinterp( PF_x )
                if key == DAtime:
                    ax.flat[idx].plot( loc_inx,loc_iny,Color[idx],linewidth=3,label=labels[ifile]+'t0',linestyle=lstyle[ifile] ) 
                else:
                    ax.flat[idx].plot( loc_inx,loc_iny,'black',linewidth=3,label=labels[ifile]+'t1',linestyle=lstyle[ifile] )

    ax.flat[2].legend(loc='upper right',fontsize='16')
    #ax.flat[2].legend(bbox_to_anchor=(1.2, 1.0),frameon=True,loc='upper right',fontsize='18')
    # subplot title and labels
    ylabel_like = [0.0,5.0,10.0,15.0,20.0]
    yticks = []
    list_y_range = list(y_range)
    for it in ylabel_like:
        yticks.append( list_y_range.index(it) ) 
    
    for var in hydros:
        idx = hydros.index( var )
        ax.flat[idx].set_title( var, fontsize = 20)
        ax.flat[idx].set_yticks( yticks )
        ax.flat[idx].set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
        ax.flat[idx].set_ylim(ymin=0,ymax=20.5) # cut off data above 25km        

    # a common y label
    fig.text(0.06,0.5,'Height (km)',ha='center',va='center',rotation='vertical',fontsize=20)
    # set X label
    ax.flat[-1].set_xticks( x_axis_rg[::10] )
    ax.flat[-1].set_xticklabels( ['0','500','1000','1500','2000','2500','3000'],fontsize=20 ) 
    ax.flat[-1].set_xlabel('Mass (kg m-2)',fontsize=20)
    ax.flat[-1].set_xlim(xmin=0)

    # Set title
    title_name = Storm+': '+Exper_name+' '+DAtime+'  \nProfile: hydrometeors in the circled area \n(center@min slp of Xa, radius=200km)'
    fig.suptitle(title_name, fontsize=15, fontweight='bold')

    # Save the figure
    save_des = plot_dir+DAtime+'_VP_each_hydro_twotimes_area.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


    return None



# Plot height profile of the sum of hydrometeors on the same plot
def plot_all_oneTime(DAtime, plot_dir, d_hydro, Height):


    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)

    # Manually set discrete values on x and y axis and interpolate data to these values
    ## x axis: correlation value
    x_range = np.arange(0,3000.5,50) #2100.5
    x_axis_rg = range(len(x_range))
    f_xinterp = interpolate.interp1d( x_range, x_axis_rg)
    ## y axis: model level height
    y_range = np.arange(0,31,1)
    y_axis_rg = range(len(y_range))
    f_yinterp = interpolate.interp1d( y_range, y_axis_rg)
    loc_iny = f_yinterp( Height )

    labels = ['Xb_','Xa_']
    lstyle = ['-','--']  
    Color = ['#f032e6','#911eb4','#4363d8','#f58231','#469990']

    for ifile in range( len(wrf_files) ):
        for var in hydros:
            idx = hydros.index( var )
            PF_x = np.sum(d_hydro[var][ifile,:,:],axis=1)
            loc_inx = f_xinterp( PF_x )
            ax.plot( loc_inx,loc_iny,Color[idx],linewidth=2,label=labels[ifile]+var,linestyle=lstyle[ifile] )

    # set lables
    ax.legend(loc='upper right',fontsize='10')
    # set X label
    ax.set_xticks( x_axis_rg[::10] )
    ax.set_xticklabels(  ['0','500','1000','1500','2000','2500','3000'],fontsize=15 ) 
    ax.set_xlabel('Mass (kg m-2)',fontsize=15)
    ax.set_xlim(xmin=0) 
    # set Y label
    ylabel_like = [0.0,5.0,10.0,15.0,20.0]
    yticks = []
    list_y_range = list(y_range)
    for it in ylabel_like:
        yticks.append( list_y_range.index(it) )
    ax.set_yticks( yticks )
    ax.set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
    ax.set_ylabel('Height (km)',fontsize=15)
    ax.set_ylim(ymin=0,ymax=20.5) # cut off data above 25km

    # Set title
    title_name = 'Profile: hydrometeors in the circled area (center=min slp of Xb, radius=200km)'
    ax.set_title( title_name,fontweight="bold",fontsize='12' )
    fig.suptitle(Storm+': '+Exper_name+' '+DAtime , fontsize=10, fontweight='bold')

    # Save the figure
    save_des = plot_dir+DAtime+'_VP_each_hydro_area.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


    return None




if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    model_resolution = 3000 #m
    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR'
    MP = 'THO'
    hydros = ['QCLOUD','QRAIN','QICE','QSNOW','QGRAUP']

    # Time range set up
    start_time_str = '201709030000'
    end_time_str = '201709050000'
    Consecutive_times = True

    deep_slp_incre = True
    incre_slp_th = 20 # threshold of increment, unit:hpa
    radius_th = 200 #km

    each_water = True
    Plot_oneTime = False
    # -------------------------------------------------------    

    # Create experiment names
    Exper_name =  UD.generate_one_name( Storm,DA,MP ) 

    if not Consecutive_times:
        DAtimes = ['201709180000',]#'201709041500']
        #DAtimes = ['201708230000','201708230600','201708231200','201708231800','201708240000','201708240600','201708241200']
        #DAtimes = ['201709031200','201709031800','201709040000','201709040600','201709041200','201709041800','201709050000']
        #DAtimes = ['201709161200','201709161800','201709170000','201709170600','201709171200','201709171800','201709180000']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]


    # Plot integrated column of hydrometeors
    if deep_slp_incre:
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/deep_slp_incre/VProfile_Hydro/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
        # ----- Read min slp from model-----------------
        HPI_models = {}
        DAtimes_dir = [big_dir+Storm+'/'+Exper_name+'/fc/'+it for it in DAtimes]
        file_kinds = ['wrf_enkf_input_d03_mean','wrf_enkf_output_d03_mean']
        for ifk in file_kinds:
            idx = file_kinds.index( ifk )
            HPI_models[ifk] = read_HPI_model( Storm, Exper_name, ifk, DAtimes_dir )
        incre_slp = np.array(HPI_models['wrf_enkf_output_d03_mean']['min_slp']) - np.array(HPI_models['wrf_enkf_input_d03_mean']['min_slp'])

    start_time=time.process_time()
    for DAtime in DAtimes:
        if deep_slp_incre:
            idx_t = DAtimes.index( DAtime )
            print('At '+DAtime)
            if abs(incre_slp[idx_t]) > incre_slp_th:
                print('At '+DAtime)
                wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
                wrf_files = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
                d_hydro,Height = compute_hydro_area( wrf_files, hydros, HPI_models, idx_t)
                # --- Condition
                if each_water:
                    v_interest = hydros
                    if Plot_oneTime:
                        plot_all_oneTime(DAtime, plot_dir, d_hydro, Height )
                    else:
                        d_hydro_times = {}
                        d_hydro_times[DAtime] = d_hydro
                        # Read out the next time
                        DAtime_next = datetime.strptime(DAtime,"%Y%m%d%H%M") + timedelta(hours=1) 
                        DAtime_next = DAtime_next.strftime("%Y%m%d%H%M") 
                        idx_t = DAtimes.index( DAtime_next )
                        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime_next+'/'
                        wrf_files = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
                        d_hydro_next,Height = compute_hydro_area( wrf_files, hydros, HPI_models, idx_t)
                        d_hydro_times[DAtime_next] = d_hydro_next
                        # plot
                        plot_all_multiTime(DAtime, plot_dir, d_hydro_times, Height ) 
                else:
                    pass
                #v_interest = ['liquid','ice','all_hydro']
                #for var in v_interest:
                #    if var == 'liquid':
                #        tmp = d_hydro['QCLOUD']+d_hydro['QRAIN']
                #    elif var == 'ice':
                #        tmp = d_hydro['QPFE']+d_hydro['QSNOW']+d_hydro['QGRAUP']
                #    elif var == 'all_hydro':
                #        tmp = d_hydro['QCLOUD']+d_hydro['QRAIN']+d_hydro['QPFE']+d_hydro['QSNOW']+d_hydro['QGRAUP']
                #    PF_xb = np.sum( tmp[0,:,:],axis=0 )
                #    PF_xa = np.sum( tmp[1,:,:],axis=0 )
                #    plot_PF(DAtime, plot_dir, d_hydro, var, PF_xb, PF_xa, )

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')


import os
import numpy as np
from datetime import datetime, timedelta
import netCDF4 as nc
import matplotlib
from scipy import interpolate
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import math
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
matplotlib.use("agg")
import pickle

import Util_data as UD
import ModelX_calculate_pert_stddev as util


def read_wrf_domain( wrf_file ):

    print('Read domain info from: ' + wrf_file)
    ncdir = nc.Dataset(wrf_file, 'r')

    Lat_x = ncdir.variables['XLAT'][0,:,:].flatten() #latitude: XLAT(time, y, x)
    Lon_x = ncdir.variables['XLONG'][0,:,:].flatten() #longitude: XLONG(time, y, x)

    lat_min = np.min( Lat_x.flatten() )
    lat_max = np.max( Lat_x.flatten() )
    lon_min = np.min( Lon_x.flatten() )
    lon_max = np.max( Lon_x.flatten() )

    d03_list = {'lon':Lon_x,'lat':Lat_x,'lat_min':lat_min, 'lat_max':lat_max, 'lon_min':lon_min, 'lon_max':lon_max}
    return d03_list


# Read in pre-calcualted stddev
def Read_EnsStddev( wrf_dir,xdim,tt,var_name ):
    if xdim == '3D':
        nLevel = 42
    elif xdim == '2D':
        nLevel = 1

    # Read ensemble standard deviation of xb
    des_path = wrf_dir+ 'xb_d03_'+xdim+'_ensStddev_' + tt + '_' +  var_name + '.pickle'
    with open( des_path,'rb' ) as f:
        stddev_xb = pickle.load( f )
    print('Shape of stddev_xb: '+ str(np.shape(stddev_xb)))
    assert not np.isnan(stddev_xb).any()

    return stddev_xb

# Plot the vertical profile of horizontal-averaged stddev
# x axis: range of stddev; y axis: vertical height
def plot_profile_hormean( std, v_3d, tt):

    fig, ax=plt.subplots(2, 3, sharey='all', figsize=(9,7), dpi=300)

    # Get unstaggered values of geo height in KM
    wrf_path = big_dir+Storm+'/'+Exper_name[MP[0]]+'/fc/'+tt+'/'
    if MakeEns:
        wrf_file = wrf_path+'wrfinput_d03'
    else:
        wrf_file = wrf_path+'wrf_enkf_output_d03'
    ncdir = nc.Dataset( wrf_file, 'r')
    ph = ncdir.variables['PH'][0,:,:,:] # perturbation
    phb = ncdir.variables['PHB'][0,:,:,:]
    geoHkm = (ph+phb)/9.81/1000 # in kilo
    geoHkm_Dmean = np.nanmean( geoHkm, axis=1 )
    geoHkm_half_eta = (geoHkm_Dmean[:-1]+geoHkm_Dmean[1:])/2
    Height = np.ma.getdata(geoHkm_half_eta)

    # Make Y axis
    y_range = np.arange(0,31,1)
    y_axis_rg = range(len(y_range))
    f_yinterp = interpolate.interp1d( y_range, y_axis_rg)
    loc_iny = f_yinterp( Height )

    # Customization
    color = {'THO':'blue','WSM6':'red'}
    labels = {'THO':'THO','WSM6':'WSM6'}

    for imp in MP:
        for ivar in v_3d:
            idx = v_3d.index(ivar)
            ax.flat[idx].plot( np.nanmean(std[imp][ivar][tt],axis=1),loc_iny,color[imp],linewidth=3)

    # Subplot title and labels
    ylabel_like = [0.0,5.0,10.0,15.0,20.0]
    yticks = []
    for iy in ylabel_like:
        yticks.append( f_yinterp( iy ) )
    for ivar in v_3d:
        idx = v_3d.index(ivar)
        ax.flat[idx].set_title( ivar, fontsize = 15)
        ax.flat[idx].set_ylim(ymin=0,ymax=20.5) # cut off data above 25km
        ax.flat[idx].set_yticks( yticks )
    ax[0,0].set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
    ax[1,0].set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )

    # a common y label
    fig.text(0.06,0.5,'Height (km)',ha='center',va='center',rotation='vertical',fontsize=20)
    # set X label
    #ax.set_xticks( x_axis_rg[::10] )
    #if Storm == 'IRMA' and MP == 'THO':
    #    ax.set_xticklabels(  ['0','500','1000','1500','2000','2500','3000','3500','4000',],fontsize=15 )
    #else:
    #    ax.set_xticklabels(  ['0','500','1000','1500','2000','2500','3000',],fontsize=15 )
        #ax.set_xticklabels(  ['0','500','1000','1500','2000',],fontsize=15 )
    #ax.set_xlabel('Mass (kg m-2)',fontsize=20)
    #ax.set_xlim(xmin=0)

    # Set title
    #title_name = Storm+': '+Exper_name+' '+DAtime+'  \nProfile: hydrometeors in the circled area \n(center@min slp of Xa, radius=200km)'
    #fig.suptitle(title_name, fontsize=12, fontweight='bold')

    # Save the figure
    if MakeEns:
        save_des = plot_dir+'test.png'
    else:
        save_des = plot_dir+tt+'/VP_'+var+'_twotimes_area.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()

# Plot the spatial distribution of standard deviation
def spatial_dis_2D( std, ivar, tt ):

    fig, ax=plt.subplots(1, 3, sharey='all', figsize=(5,2.5), dpi=300)

    # Get domain range
    wrf_path = big_dir+Storm+'/'+Exper_name[MP[0]]+'/fc/'+tt+'/'
    if MakeEns:
        wrf_file = wrf_path+'wrfinput_d03'
    else:
        wrf_file = wrf_path+'wrf_enkf_output_d03_mean'
    d_wrf_d03 = read_wrf_domain( wrf_file )


    # Plot
    fig, axs=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,2.5), dpi=400)

    # Define the domain
    lon = d_wrf_d03['lon']
    lat = d_wrf_d03['lat']
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    # Set the map
    for i in range(3):
        axs.flat[i].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        axs.flat[i].coastlines(resolution='10m', color='black',linewidth=0.5)

    # Ens 1
    min_std = 0
    max_std = 1.5
    Ens1_hPa = std[MP[0]][ivar][tt] / 100
    la1 = axs.flat[0].scatter(lon,lat,1.5,Ens1_hPa,vmin=min_std, vmax=max_std,cmap='binary',transform=ccrs.PlateCarree())
    caxes = fig.add_axes([0.12, 0.1, 0.25, 0.02])
    obs_bar = fig.colorbar(la1,ax=axs[0],orientation="horizontal", cax=caxes)
    obs_bar.ax.tick_params(labelsize=6)

    Ens2_hPa = std[MP[1]][ivar][tt] / 100
    la2 = axs.flat[1].scatter(lon,lat,1.5,Ens2_hPa,cmap='binary',vmin=min_std,vmax=max_std,transform=ccrs.PlateCarree())
    caxes = fig.add_axes([0.52, 0.1, 0.25, 0.02])
    obs_bar = fig.colorbar(la2,ax=axs[1],orientation="horizontal", cax=caxes)
    obs_bar.ax.tick_params(labelsize=6)
    
    axs.flat[2].scatter(lon,lat,1.5,Ens1_hPa - Ens2_hPa,cmap='seismic',transform=ccrs.PlateCarree())

    #obs_s = axs.flat[0].scatter(d_all['lon_obs'],d_all['lat_obs'],1.5,c=d_all['Yo_obs'],edgecolors='none', cmap=IRcmap, vmin=min_obs, vmax=max_obs,transform=ccrs.PlateCarree())
    #if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
    #    axs.flat[0].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='darkviolet', transform=ccrs.PlateCarree())
    #axs.flat[0].add_patch(patches.Polygon(path,facecolor='none',edgecolor='white',linewidth=0.5 ))
    # Colorbar
    #caxes = fig.add_axes([0.12, 0.1, 0.25, 0.02])
    #obs_bar = fig.colorbar(obs_s,ax=axs[0],orientation="horizontal", cax=caxes)
    #obs_bar.ax.tick_params(labelsize=6)

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(3):
        gl = axs[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.7, linestyle='--')

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
        gl.xlabel_style = {'size': 4}
        gl.ylabel_style = {'size': 6}

    if MakeEns:
        save_des = plot_dir+'test2D.png'
    else:
        save_des = plot_dir+tt+'/VP_'+var+'_twotimes_area.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()





if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'HARVEY'
    DA = 'CONV'
    MP = ['THO','WSM6']
    # variables of interest
    v_3d = [ 'U','V','W','T','P','QVAPOR',] #[ 'PSFC']
    v_2d = ['PSFC']
    v_all = v_3d + v_2d
    # Which stage
    MakeEns = True # Ens generated from perturbing GFS field
    if MakeEns:
        times = util.set_time_for_MakeEns( Storm )
    else:
        start_time_str = '201709160000'
        end_time_str = '201709160000'
        Consecutive_times = True
        times = util.set_time_for_cyclings( Storm, start_time_str,end_time_str )
    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    # operation:
    profile_hor_mean = False
    plot_2D = True

    # -------------------------------------------------------
    Exper_name = {}
    for imp in MP:
        Exper_name[imp] = UD.generate_one_name( Storm,DA,imp )

    # Read in pre-calculated stddev
    std = {}
    for imp in MP:
        std[imp] = {}
        for ivar in v_all:
            var_dim = UD.def_vardim( ivar )
            std[imp][ivar] = {}
            for tt in times:
                wrf_dir = big_dir+Storm+'/'+Exper_name[imp]+'/fc/'+tt+'/'
                xdim = UD.def_vardim(ivar)
                std[imp][ivar][tt] = Read_EnsStddev( wrf_dir,xdim,tt,ivar )

    if not MakeEns:
        #print('------------ Calculate the ensemble perturbations for EnKF cyclings--------------')
        for DAtime in times:
            pass

    else:
        # set plot dir
        for imp in MP:
            plot_dir = small_dir+Storm+'/'+Exper_name[imp]+'/Vis_analyze/CV3/'
            plotdir_exists = os.path.exists( plot_dir )
            if plotdir_exists == False:
                os.mkdir(plot_dir)

        for tt in times:
            
            if profile_hor_mean:
                plot_profile_hormean( std, v_3d, tt )

            if plot_2D:
                for ivar in v_2d:
                    spatial_dis_2D( std, ivar, tt )
    

            #for var_name in v_interest:
            #    print('Working on '+var_name+' ......')
            #var_dim = UD.def_vardim( var_name )




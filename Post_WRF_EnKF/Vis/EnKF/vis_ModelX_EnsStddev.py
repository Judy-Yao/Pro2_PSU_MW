
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
    if os.path.exists( des_path ):
        with open( des_path,'rb' ) as f:
            stddev_xb = pickle.load( f )
        print('Shape of stddev_xb: '+ str(np.shape(stddev_xb)))
        assert not np.isnan(stddev_xb).any()
        return stddev_xb
    else:
        return None

# Plot the vertical profile of horizontal-averaged stddev
# x axis: range of stddev; y axis: vertical height
def plot_profile_hormean( std, v_3d, tt):

    fig, ax=plt.subplots(2, 3, sharey='all', figsize=(9,7), dpi=300)

    # Get unstaggered values of geo height in KM
    wrf_path = big_dir+Storm+'/'+Exper_name[MP[0]][DA[0]]+'/fc/'+tt+'/'
    if MakeEns:
        wrf_file = wrf_path+'wrfinput_d03'
    else:
        wrf_file = wrf_path+'wrf_enkf_output_d03_mean'
    ncdir = nc.Dataset( wrf_file, 'r')
    ph = ncdir.variables['PH'][0,:,:,:] # perturbation
    phb = ncdir.variables['PHB'][0,:,:,:]
    geoHkm = (ph+phb)/9.81/1000 # in kilo
    geoHkm = geoHkm.reshape( geoHkm.shape[0],-1)
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
    linestyles = {'CONV':'-','CONV-WSM6Ens':':'}
    labels = {'THO':'THO','WSM6':'WSM6'}

    for imp in MP:
        for ida in DA:
            for ivar in v_3d:
                idx = v_3d.index(ivar)
                if std[imp][ida][ivar][tt] is None:
                    continue
                dmean_std = np.nanmean(std[imp][ida][ivar][tt],axis=1)
                if ivar == 'Pres':
                    dmean_std = dmean_std/100 # pa to hPa
                if 'Q' in ivar:
                    dmean_std = dmean_std*1000 # kg/kg to g/kg
                # plot
                ax.flat[idx].plot( dmean_std,loc_iny,color[imp],linestyle=linestyles[ida],linewidth=2.5)
       
    lines = ax.flat[-1].get_lines()
    lgd = MP+['THO-WSM6Ens']
    legend0 = ax.flat[-1].legend(lines,lgd, fontsize='15', loc='upper right')
    # Add the first legend manually to the current Axes
    ax.flat[-1].add_artist(legend0)

    # limit
    for ivar in v_3d:
        idx = v_3d.index(ivar)
        if ivar == 'U' or ivar == 'V':
            ax.flat[idx].set_xlim(xmin=1.5,xmax=5)
        elif ivar == 'W':
            ax.flat[idx].set_xlim(xmin=-0.01,xmax=0.3)
        elif ivar == 'Temp' or ivar == 'Pres':
            ax.flat[idx].set_xlim(xmin=-0.01,xmax=1.5)
        elif ivar == 'QVAPOR':
            ax.flat[idx].set_xlim(xmin=-0.01,xmax=2.5)

    # Subplot title and labels
    ylabel_like = [0.0,5.0,10.0,15.0,20.0]
    yticks = []
    for iy in ylabel_like:
        yticks.append( f_yinterp( iy ) )
    for ivar in v_3d:
        idx = v_3d.index(ivar)
        if ivar == 'U' or ivar == 'V' or ivar == 'W':
            title = ivar +' (m/s)'
        elif ivar == 'Temp':
            title = ivar + ' (k)'
        elif ivar == 'Pres':
            title = ivar + ' (hPa)'
        elif ivar == 'QVAPOR':
            title = ivar +' (g/kg)'
        ax.flat[idx].set_title( title, fontsize = 15)
        ax.flat[idx].set_ylim(ymin=0,ymax=20.5) # cut off data above 25km
        ax.flat[idx].set_yticks( yticks )
    ax[0,0].set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
    ax[1,0].set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )

    # a common y label
    fig.text(0.03,0.5,'Height (km)',ha='center',va='center',rotation='vertical',fontsize=20)
    #ax.set_xlabel('Mass (kg m-2)',fontsize=20)
    #ax.set_xlim(xmin=0)

    # Set title
    if MakeEns:
        title_name = Storm+'--GFS perturbed ensemble: profiles of domain-averaged ensemble spread'
    else:
        title_name = Storm+' '+tt+'\nProfiles of domain-averaged ensemble spread'

    fig.suptitle(title_name, fontsize=14, fontweight='bold')

    # Save the figure
    if MakeEns:
        save_des = plot_dir+'ens_spread_profiles.png'
    else:
        save_des = plot_dir+'ens_spread_profiles_'+tt+'.png'
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
    if MakeEns:
        min_std = 0.5
        max_std = 1.0
    else:
        if Storm == 'IRMA':
            min_std = 0.5
            max_std = 10.0
        else:
            min_std = 0.5
            max_std = 1.3

    Ens1_hPa = std[MP[0]][ivar][tt] / 100
    axs.flat[0].scatter(lon,lat,1.5,Ens1_hPa,vmin=min_std, vmax=max_std,cmap='binary',transform=ccrs.PlateCarree())
    #axs.flat[0].scatter(lon,lat,1.5,Ens1_hPa,cmap='binary',transform=ccrs.PlateCarree())
    Ens2_hPa = std[MP[1]][ivar][tt] / 100
    #la1 = axs.flat[1].scatter(lon,lat,1.5,Ens2_hPa,cmap='binary',transform=ccrs.PlateCarree())
    la1 = axs.flat[1].scatter(lon,lat,1.5,Ens2_hPa,cmap='binary',vmin=min_std,vmax=max_std,transform=ccrs.PlateCarree())
    caxes = fig.add_axes([0.14, 0.1, 0.5, 0.02])
    obs_bar = fig.colorbar(la1,ax=axs[0],orientation="horizontal", cax=caxes,extend='both')
    obs_bar.ax.tick_params(labelsize=6)

    if MakeEns:
        min_std_diff = -0.1
        max_std_diff = 0.1
    else:
        min_std_diff = -0.5
        max_std_diff = 0.5
    la2 = axs.flat[2].scatter(lon,lat,1.5,Ens1_hPa-Ens2_hPa,vmin=min_std_diff,vmax=max_std_diff,cmap='bwr',transform=ccrs.PlateCarree())
    #la2 = axs.flat[2].scatter(lon,lat,1.5,Ens1_hPa - Ens2_hPa,cmap='bwr',transform=ccrs.PlateCarree())
    caxes = fig.add_axes([0.65, 0.1, 0.25, 0.02])
    obs_bar = fig.colorbar(la2,ax=axs[1],orientation="horizontal", cax=caxes,extend='both')
    obs_bar.ax.tick_params(labelsize=6)

    # Mark the location of TCvital
    tc_lon_st, tc_lat_st, tc_slp_st = UD.read_TCvitals(small_dir, Storm, tt)
    for i in range(3):
        axs[i].plot(tc_lon_st,tc_lat_st, '*', color='green', markersize=3, transform=ccrs.PlateCarree())

    # subplot title
    axs[0].set_title( MP[0], fontsize = 8)
    axs[1].set_title( MP[1], fontsize = 8)
    axs[2].set_title( MP[0]+' - '+MP[1], fontsize = 8)

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

    # Set title
    if MakeEns:
        title_name = Storm+':GFS perturbed ensemble\nEnsemble spread of PSFC (hPa)'
    else:
        title_name = Storm+' '+tt+'\nEnsemble spread of PSFC (hPa)'
    fig.suptitle(title_name, fontsize=8, fontweight='bold')

    # Save the figure
    if MakeEns:
        save_des = plot_dir+'ens_spread_PSFC.png'
    else:
        save_des = plot_dir+'ens_spread_PSFC_'+tt+'.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = ['CONV','CONV-WSM6Ens']
    MP = ['THO','WSM6',]
    # variables of interest
    v_3d = [ 'U','V','W','Temp','Pres','QVAPOR',] #[ 'PSFC']
    v_2d = ['PSFC']
    v_all = v_3d + v_2d
    # Which stage
    MakeEns = False # Ens generated from perturbing GFS field
    if MakeEns:
        times = util.set_time_for_MakeEns( Storm )
    else:
        EndSpinup = True
        if EndSpinup:
            times = util.set_time_for_EndSpinup( Storm )
        else:
            start_time_str = '201709030000'
            end_time_str = '201709030000'
            Consecutive_times = True
            times = util.set_time_for_cyclings( Storm,Consecutive_times,start_time_str,end_time_str)
    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    # operation:
    profile_hor_mean = True
    plot_2D = False

    # -------------------------------------------------------
    Exper_name = {}
    for imp in MP:
        Exper_name[imp] = {}
        for ida in DA:
            Exper_name[imp][ida] = UD.generate_one_name( Storm,ida,imp )

    # Read in pre-calculated stddev
    std = {}
    for imp in MP:
        std[imp] = {}
        for ida in DA:
            std[imp][ida] = {}
            for ivar in v_all:
                var_dim = UD.def_vardim( ivar )
                std[imp][ida][ivar] = {}
                for tt in times:
                    if Exper_name[imp][ida] is not None:
                        wrf_dir = big_dir+Storm+'/'+Exper_name[imp][ida]+'/fc/'+tt+'/'
                        xdim = UD.def_vardim(ivar)
                        std[imp][ida][ivar][tt] = Read_EnsStddev( wrf_dir,xdim,tt,ivar )
                    else:
                        std[imp][ida][ivar][tt] = None

    if (not MakeEns) and (not EndSpinup):
        # set plot dir
        pass

    elif MakeEns or ((not MakeEns) and (EndSpinup)):
        # set plot dir
        for imp in MP:
            plot_dir = small_dir+Storm+'/Vis_analyze/CV3/'
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




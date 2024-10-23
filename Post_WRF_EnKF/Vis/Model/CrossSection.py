
import pyproj
import os # functions for interacting with the operating system
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
from scipy import interpolate
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
import math
import warnings
from matplotlib.cm import get_cmap
from wrf import getvar, to_np, interp2dxy, ll_to_xy, xy_to_ll

import Util_wrfpython as UWP
import Util_data as UD


def distance_to_degrees_lat(distance_km):
    """
    Converts distance in kilometers to degrees of latitude.
    """
    return distance_km / 111.32

def distance_to_degrees_lon(distance_km, latitude):
    """
    Converts distance in kilometers to degrees of longitude at a given latitude.
    """
    return distance_km / (111.32 * math.cos(math.radians(latitude)))

def var_cross( var_name,wrf_files ):

    d_var = {}
    d_var['input'] = {}
    d_var['output'] = {}
    
    # ----- Create a cross line and Identify the start and end points -------
    ncdir = nc.Dataset(wrf_files[1], 'r') # enkf output
    lat = ncdir.variables['XLAT'][:].flatten()
    lon = ncdir.variables['XLONG'][:].flatten()
    # Find the storm center
    # sea level pressure
    slp = getvar(ncdir, 'slp')
    # original SLP
    slp_values = slp.values
    slp_values[slp_values > 1030] = np.nan
    slp_smt_values = sp.ndimage.gaussian_filter(slp, [11,11]) #[11,11]
    minslp = np.nanmin( slp_values )
    idx = np.nanargmin( slp_smt_values )
    # identify the geolocation of the pivot point
    lat_c = lat[idx]
    lon_c = lon[idx]
    center = [lat_c,lon_c]
    # Create the start point and end point in lon/lat for the cross section
    if along_lon: # cross section is along the longitude centered at one latitude
        if cut_segment:
            st_lat = lat_c # distance_to_degrees_lat(distance_km)
            ed_lat = lat_c
            lon_RtoD = distance_to_degrees_lon(radius, lat_c)
            st_lon = lon_c - lon_RtoD
            ed_lon = lon_c + lon_RtoD
        else:
            st_lat = lat_c
            ed_lat = lat_c
            st_lon = np.amin(lon) 
            ed_lon = np.amax(lon)
    else:
        if cut_segment:
            st_lon = lon_c # distance_to_degrees_lat(distance_km)
            ed_lon = lon_c
            lat_RtoD = distance_to_degrees_lat( radius )
            st_lat = lat_c - lat_RtoD
            ed_lat = lat_c + lat_RtoD
        else:
            st_lon = lon_c
            ed_lon = lon_c
            st_lat = np.amin(lat)
            ed_lat = np.amax(lat)
    
    # Create the  start point and end point in i/j
    st_ij = ll_to_xy(ncdir, st_lat, st_lon) # What ll_to_xy returns is not the xy coordinate itself but the grid index starting from 0. 
    st_i = st_ij.values[0]
    st_j = st_ij.values[1]
    ed_ij = ll_to_xy(ncdir, ed_lat, ed_lon)
    ed_i = ed_ij.values[0]
    ed_j = ed_ij.values[1]
    st_p = (st_i,st_j)
    ed_p = (ed_i,ed_j) 
    #xlat = ncdir.variables['XLAT'][0,:,:].flatten()  #[0,:,0]
    #xlon = ncdir.variables['XLONG'][0,:,:].flatten()

    # ----- Read variable values and Obtain values along the cross section -----
    # loop thro enkf input and output
    for wrf_file in wrf_files:
        print('Reading '+wrf_file)
        ncdir = nc.Dataset(wrf_file, 'r')
    
        # vertical coordinate
        if use_pressure:
            PB = ncdir.variables['PB'][0,:,:,:]
            P = ncdir.variables['P'][0,:,:,:]
            ver_coor_tmp = (PB + P)/100 # in hPa
        else:
            PHB = ncdir.variables['PHB'][0,:,:,:]
            PH = ncdir.variables['PH'][0,:,:,:]
            geoHkm = (PHB+PH)/9.8/1000 # in km
            geoHkm_half = (geoHkm[:-1,:,:]+geoHkm[1:,:,:])/2
            ver_coor_tmp = np.ma.getdata(geoHkm_half)

        # Get a line in ij space connecting the start and end point 
        xy_line = UWP.xy_edit(ver_coor_tmp, start_point=st_p, end_point=ed_p) # return: [[i,j]]
        ver_coor = interp2dxy(ver_coor_tmp, xy_line)
        ver_coor = np.mean(ver_coor,axis=1)

        # Variable of interest
        if var_name == 'hroi_wind':
            u_tmp = ncdir.variables['U'][0,:,:,:]
            u_var = (u_tmp[:,:,:-1]+u_tmp[:,:,1:])/2
            v_tmp = ncdir.variables['V'][0,:,:,:]
            v_var = (v_tmp[:,:-1,:]+v_tmp[:,1:,:])/2
            var = (u_var ** 2 + v_var ** 2) ** 0.5
            var_cross = interp2dxy(var, xy_line)
        elif var_name == 'W':
            tmp = ncdir.variables['W'][0,:,:,:]
            var = (tmp[:-1,:,:]+tmp[1:,:,:])/2
            var_cross = interp2dxy(var, xy_line)
        elif var_name == 'RH':
            # full pressure
            p = ncdir.variables['P'][0,:,:,:] # perturbation
            pb = ncdir.variables['PB'][0,:,:,:]
            full_p = (p + pb).filled(np.nan)
            itp_p = interp2dxy(full_p, xy_line)
            itp_p_np = itp_p.values
            # full T in kelvin
            theta = ncdir.variables['T'][0,:,:,:] # theta perturbation
            full_theta = (theta + 300).filled(np.nan)
            itp_theta = interp2dxy(full_theta, xy_line)
            itp_theta_np = itp_theta.values
            # temperature in Kelvin
            tmp_tk = UD.compute_tk( itp_p_np, itp_theta_np )
            itp_tk = tmp_tk.reshape(itp_p.shape)
            # qvapor
            qvapor = ncdir.variables['QVAPOR'][0,:,:,:]
            qvapor_filled = qvapor.filled(np.nan)  # Replaces masked values with NaN
            itp_qv = interp2dxy(qvapor_filled, xy_line)
            # relative humidity
            var = UD.njit_compute_rh( itp_p_np,itp_tk,itp_qv.values )
            var_cross = var.reshape(itp_p.shape)
        elif var_name == 'QVAPOR':
            var = getvar(ncdir, 'QVAPOR')
            var_cross = interp2dxy(var, xy_line)
        elif var_name == 'PT':
            tmp = ncdir.variables['T'][0,:,:,:] # potential temperature 
            var = tmp + 300 # potential temperature  
            var_cross = interp2dxy(var, xy_line)

        # Transform from ij space to lat/lon space
        ll_pair = [] #lat/lon
        for it in xy_line:
            lat_lon = xy_to_ll(ncdir, it[0], it[1])
            ll_pair.append( to_np(lat_lon) )
    
        if 'input' in wrf_file:
            d_var['input']= {var_name: var_cross, 'ver_coor':ver_coor, 'll_pair':ll_pair, 'center':center}
        elif 'output' in wrf_file:
            d_var['output']= {var_name: var_cross, 'ver_coor':ver_coor, 'll_pair':ll_pair, 'center':center}
    
    # Plot
    if if_plot:
        plot_var_cross( wrf_files,var_name,d_var )


def plot_var_cross( wrf_files,var_name,d_var ):

    # Create the figure
    fig, ax=plt.subplots(1, 2, sharex='all', linewidth=0.5, figsize=(12,6), dpi=400)

    # Make the contour plot
    for wrf_file in wrf_files:
        i = wrf_files.index(wrf_file)
        # x axis
        center = d_var['input']['center'] #center = [lat_c,lon_c]
        if 'input' in wrf_file:
            ll_pair = d_var['input']['ll_pair']
            ver_coor = d_var['input']['ver_coor']
        else:
            ll_pair = d_var['output']['ll_pair']
            ver_coor = d_var['output']['ver_coor']
        
        x_axis_rg = range(len(ll_pair))

        dis_center = []
        for it in ll_pair:
            dis_center.append( (it[0]-center[0])**2+(it[1]-center[1])**2)
        idx_center = np.argmin( dis_center )

        # y axis: model vertical coordinate
        if use_pressure:
            y_range = np.arange(900,50,50)
        else:
            y_range = np.arange(0,31,1)
        y_axis_rg = range(len(y_range))
        f_yinterp = interpolate.interp1d( y_range, y_axis_rg)
        yv = f_yinterp( ver_coor )
        # Make a mesh grid
        xcoor, ycoor = np.meshgrid( x_axis_rg, yv )

        # Set colorbar
        if 'radial_wind' in var_name:
            bounds = np.linspace(-12,12,9)
            cmap = 'PiYG_r'
        elif 'tangential_wind' in var_name:
            color_intervals = list(np.linspace(0,50,6))
            color_intervals.insert(0,-10.0)
            exist_cmap = plt.cm.jet
            colors = exist_cmap(np.linspace(0,1,len(color_intervals)))
            cmap = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals))
            bounds = color_intervals
        elif 'hroi_wind' in var_name:
            bounds = np.linspace(10,70,7)
            cmap = 'jet'
        elif 'W' in var_name:
            bounds = [-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8]
            cmap = 'PiYG'
        elif 'Q' in var_name:
            bounds = 10.**np.arange(-6,0,1)
            cmap = 'ocean_r'
        elif 'RH' in var_name:
            bounds = [0,10, 20, 30, 40, 60, 70, 80, 90,100]
            cmap = mcolors.LinearSegmentedColormap.from_list(
                'custom_cmap',
                [
                    (0.0, '#8b6c5c'),  # Darker brown for values below 10%
                    (1/7, '#bca89f'),  # Dark brown at 10%
                    (2/7, '#d8cbc4'),  # Light brown at 40%
                    (3/7, '#ffffff'),  # White at 40-60%
                    (4/7, '#cce7c9'),  # Light green at 60%
                    (5/7, '#8bca84'),  # Dark green at 90%
                    (6/7, '#5bb450'),   # Darker green for values above 90%
                    (1.0, '#276221')
                ]
            )
            norm = mcolors.BoundaryNorm(bounds, cmap.N)
        elif 'PT' in var_name:
            bounds = [290,300,310,320,330]
            cmap = 'jet'

        # Plot
        if 'input' in wrf_file:
            value = d_var['input'][var_name]
        else:
            value = d_var['output'][var_name] 
        
        if 'RH' in var_name:
            xc = ax[i].contourf( xcoor,ycoor,to_np(value),cmap=cmap,norm=norm,levels=bounds)
            # Add the color bar
            cbaxes = fig.add_axes([0.92, 0.1, 0.03, 0.8])
            cbar = fig.colorbar(xc,cax=cbaxes,fraction=0.046, pad=0.04)
            cbar.ax.tick_params(labelsize=15)
        else:
            xc = ax[i].contourf( xcoor,ycoor,to_np(value),cmap=cmap,levels=bounds,extend='both')
            #xc = ax[i].contourf( xcoor,ycoor,to_np(value),cmap=cmap,)
            # Add the color bar
            cbaxes = fig.add_axes([0.92, 0.1, 0.03, 0.8])
            cbar = fig.colorbar(xc,cax=cbaxes,fraction=0.046, pad=0.04)
            cbar.ax.tick_params(labelsize=15)

        ax[i].axvline(x=idx_center,color='k',linestyle='-',linewidth=2)

        # Set the x-ticks to use latitude and longitude labels.
        x_ticks = x_axis_rg
        fmt="{:.1f}, {:.1f}"
        x_labels = [fmt.format(pair[0],pair[1]) for pair in np.array(ll_pair)]
        ax[i].set_xticks(x_ticks[::40])
        ax[i].set_xticklabels(x_labels[::40], rotation=15, fontsize=12)

        f_yinterp = interpolate.interp1d( y_range, y_axis_rg)
        y_ticks = f_yinterp( ver_coor )
    
        # Set the y-ticks
        if use_pressure:
            pass
        else: 
            ylabel_like = [0.0,5.0,10.0,15.0,20.0]
        yticks = []
        list_y_range = list(y_range)
        for it in ylabel_like:
            yticks.append( f_yinterp( it ) )
        ax[i].set_yticks( yticks )
        ax[i].set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
        ax[i].set_ylim(ymin=0,ymax=20) # cut off data above 25km   


    if 'hroi_wind' in var_name:
        xb_title = 'Xb:'+var_name+' (m/s)'
        xa_title = 'Xa:'+var_name+' (m/s)'
    elif 'W' in var_name:
        xb_title = 'Xb: vertical velocity(m/s)'
        xa_title = 'Xa: vertical velocity (m/s)'
    elif 'Q' in var_name:
        xb_title = 'Xb:'+var_name+' (kg/kg)'
        xa_title = 'Xa:'+var_name+' (kg/kg)'
    elif 'RH' in var_name:
        xb_title = 'Xb:'+var_name+' (%)'
        xa_title = 'Xa:'+var_name+' (%)'    
    elif 'PT' in var_name:
        xb_title = 'Xb: Potential Temp (K)'
        xa_title = 'Xa: Potential Temp (K)'

    ax[0].set_title( xb_title, fontsize = 15, fontweight='bold')
    ax[1].set_title( xa_title, fontsize = 15, fontweight='bold')

    # Set the x-axis and  y-axis labels
    fig.text(0.04,0.5,'Height (km)',ha='center',va='center',rotation='vertical',fontsize=20)
    fig.text(0.5,0.02,"Latitude, Longitude (degree)",ha='center',va='center',rotation='horizontal', fontsize=13)

    # title
    title_name = Storm+': '+Exper_name+' '+DAtime
    fig.suptitle(title_name, fontsize=15, fontweight='bold')

    # Save figures
    if along_lon:
        figure_des=plot_dir+DAtime+'_'+var_name+'_CrossSection_alongLon.png'
    else:
        figure_des=plot_dir+DAtime+'_'+var_name+'_CrossSection_alongLat.png'
    plt.savefig(figure_des, dpi=400)
    print('Saving the figure: ', figure_des)

    plt.close()


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    MP = 'WSM6'
    DA = 'CONV'
    v_interest = ['RH'] #['hroi_wind','RH']

    start_time_str = '201709030000'
    end_time_str = '201709031200'
    Consecutive_times = True

    # model dimension
    xmax = 297
    ymax = 297
    nLevel = 42

    # Specify cross section location
    cut_segment = False
    radius = 200 #km
    along_lon = False

    use_pressure = False
    if_plot = True
    # ------------------------------------------------------- 
    
    Exper_name = UD.generate_one_name( Storm,DA,MP )

    # Identify DA times in the period of interest
    if not Consecutive_times:
        DAtimes = ['201709140000',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # create plot_dir
    if if_plot:
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/CrossSection/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

    # Plot cross section of a variable    
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
        wrf_files = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']

    # Loop through each DAtime/analysis
    for DAtime in DAtimes:
        wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime
        print('Reading WRF background and analysis: ', wrf_dir)
        DAtime_dt = datetime.strptime( DAtime, '%Y%m%d%H%M' )

        for var_name in v_interest:
            # ------ Plot -------------------
            plot_dir = plot_dir+var_name
            plotdir_exists = os.path.exists( plot_dir )
            if plotdir_exists == False:
                os.mkdir(plot_dir)

                var_cross( var_name,wrf_files ) 

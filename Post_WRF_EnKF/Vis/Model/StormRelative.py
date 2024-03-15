#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import pyproj
import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
#from wrf import getvar, interplevel
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
import matlab.engine
import pickle
from matplotlib.colors import LinearSegmentedColormap

import Util_data as UD

def test_verify_coordinate_conversion( wrf_file ):

    d_wrf_d03 = UD.read_wrf_domain( wrf_file )
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    fig, ax=plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(12,6), dpi=400)

    for i in range(2):
        ax[i].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i].coastlines (resolution='10m', color='black', linewidth=1)

    min_wind = -10
    max_wind = 10
    cs =ax[0].scatter(xlon,xlat,2,c=u10,cmap='RdBu_r',vmin=min_wind,vmax=max_wind,transform=ccrs.PlateCarree())
    ax[1].scatter(lon_polar,lat_polar,2,c=u_out,cmap='RdBu_r',vmin=min_wind,vmax=max_wind,transform=ccrs.PlateCarree())

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for j in range(2):
        gl = ax[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
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
    
    plt.savefig( 'test.png', dpi=300 )
    print('Saving the figure!')
    plt.close()



def azimuth_range_to_lat_lon(azimuths, radius, center_lon, center_lat, geod=None):
    """Convert azimuth and range locations in a polar coordinate system to lat/lon coordinates.

    Pole refers to the origin of the coordinate system.

    Parameters
    ----------
    azimuths : array-like
        array of azimuths defining the grid. If not a `pint.Quantity`,
        assumed to be in degrees.
    radius : array-like
        array of range distances from the pole. Typically in meters.
    center_lat : float
        The latitude of the pole in decimal degrees
    center_lon : float
        The longitude of the pole in decimal degrees
    geod : `pyproj.Geod` or ``None``
        PyProj Geod to use for forward azimuth and distance calculations. If ``None``, use a
        default spherical ellipsoid.

    Returns
    -------
    lon, lat : 2D arrays of longitudes and latitudes corresponding to original locations

    Notes
    -----
    Credit to Brian Blaylock for the original implementation.

    """
    g = pyproj.Geod(ellps='sphere') if geod is None else geod
    try:  # convert range units to meters
        radius = radius.m_as('meters')
    except AttributeError:  # no units associated
        warnings.warn('Range values are not a Pint-Quantity, assuming values are in meters.')
    try:  # convert azimuth units to degrees
        azimuths = azimuths.m_as('degrees')
    except AttributeError:  # no units associated
        warnings.warn(
            'Azimuth values are not a Pint-Quantity, assuming values are in degrees.'
        )
    rng2d, az2d = np.meshgrid(radius, azimuths)
    lats = np.full(az2d.shape, center_lat)
    lons = np.full(az2d.shape, center_lon)
    lon, lat, _ = g.fwd(lons, lats, az2d, rng2d)

    return lon, lat


def geodetic_to_polar( wrf_file ):

    ncdir = nc.Dataset(wrf_file, 'r')
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()  
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()    
    # Read vertical coordinate 
    if use_pressure:
        PB = ncdir.variables['PB'][0,:,:,:]
        P = ncdir.variables['P'][0,:,:,:]
        P_hpa = (PB + P)/100 # in hPa
        tmp_coor = P_hpa.reshape( P_hpa.shape[0],-1)
    else:
        PHB = ncdir.variables['PHB'][0,:,:,:]
        PH = ncdir.variables['PH'][0,:,:,:]
        geoHkm = (PHB+PH)/9.8/1000 # in km
        tmp_coor = geoHkm.reshape( geoHkm.shape[0],-1)

    # Find the polar origin (usually the same as the storm center)
    if Storm == 'MARIA' and DAtime == '201709160000':
        lon_c = -45.8
        lat_c = 11.6
    else:
        slp = UD.compute_slp( ncdir ) # from enkf output
        min_slp = np.min( slp )
        slp_smooth = sp.ndimage.filters.gaussian_filter(slp, [11, 11] )
        idx = np.nanargmin( slp_smooth )
        lat_c = ncdir.variables['XLAT'][:].flatten()[idx]
        lon_c = ncdir.variables['XLONG'][:].flatten()[idx]

    # Set up a numerical polar coordinate defined by azimuth angles and radius
    # Obtain the lon lat in this numerical polor coordinate
    lon_polar,lat_polar = azimuth_range_to_lat_lon(azimuths, radius*1000, lon_c, lat_c) # radius: convert from km to meter

    # Select grids in the geodetic coordinate that not exceed the polar coordinate
    idx_i = list(set( np.where(xlon>=lon_polar.min())[0] )&set( np.where(xlon<=lon_polar.max())[0] ))
    idx_j = list(set( np.where(xlat>=lat_polar.min())[0] )&set( np.where(xlat<=lat_polar.max())[0] ))
    idx_x = list(set(idx_i)&set(idx_j))
    xlon = xlon[idx_x] 
    xlat = xlat[idx_x] 

    # Prepare the grid in different coordinates for future interpolation
    grid_geo = np.concatenate([xlon.reshape((len(idx_x),1)),xlat.reshape((len(idx_x),1))],axis=1)
    grid_polar = np.concatenate([lon_polar.reshape(-1,1),lat_polar.reshape(-1,1)],axis=1)

    # Prepare vertical coordinate
    if use_pressure:
        tmp_coor = tmp_coor[:,idx_x]
        ver_coor = np.mean( tmp_coor,axis=1 )
    else:
       tmp_coor = tmp_coor[:,idx_x]
       tmp_coor = np.mean( tmp_coor, axis=1 )
       ver_coor = (tmp_coor[:-1]+tmp_coor[1:])/2
       ver_coor = np.ma.getdata(ver_coor) 

    d_coor = {'geo_index':idx_x,'geo_coor':grid_geo,'polar_coor':grid_polar,'ver_coor':ver_coor}
    return d_coor

def var_geodetic_to_polar( var_name,d_coor,wrf_files ):

    d_var = {}
    d_var['input'] = {}
    d_var['output'] = {}

    # coordinate info
    idx_x = d_coor['geo_index']
    grid_geo = d_coor['geo_coor']
    grid_polar = d_coor['polar_coor']

    # get lon lat list
    lon_geo = []
    lat_geo = []
    for i in grid_geo:
        lon_geo.append(i[0])
        lat_geo.append(i[1])
   
    lon_polar = []
    lat_polar = []
    for j in grid_polar:
        lon_polar.append(j[0])
        lat_polar.append(j[1]) 

    # Start a matlab process
    eng = matlab.engine.start_matlab()

    # loop thro enkf input and output
    for wrf_file in wrf_files:
        print('Reading '+wrf_file)
        ncdir = nc.Dataset(wrf_file, 'r') 

        if 'hor_wind' in var_name:
            # Read u and v wind
            u = ncdir.variables['U'][0,:,:,:]
            uwnd = (u[:,:,:-1]+u[:,:,1:])/2
            uwnd = uwnd.reshape( uwnd.shape[0],-1)
            u_w = uwnd[:,idx_x]
            v = ncdir.variables['V'][0,:,:,:]
            vwnd = (v[:,:-1,:]+v[:,1:,:])/2
            vwnd = vwnd.reshape( vwnd.shape[0],-1)
            v_w = vwnd[:,idx_x]
            # Interpolate var in geodetic coordinate to polar coordinate
            u_polar = np.zeros((nLevel,len(azimuths),len(radius)))
            v_polar = np.zeros((nLevel,len(azimuths),len(radius)))
            for il in range(nLevel):
                Mu_out = eng.griddata(matlab.double(lon_geo), matlab.double(lat_geo), matlab.double(u_w[il,:].tolist()), matlab.double(lon_polar), matlab.double(lat_polar) )
                u_out = np.array(Mu_out._data)
                u_polar[il,:,:] = u_out.reshape((len(azimuths),len(radius)))
                Mv_out = eng.griddata(matlab.double(lon_geo), matlab.double(lat_geo), matlab.double(v_w[il,:].tolist()), matlab.double(lon_polar), matlab.double(lat_polar) )
                v_out = np.array(Mv_out._data)
                #v_out = interpolate.griddata(grid_geo,v_w[il,:],grid_polar,method='linear')
                v_polar[il,:,:] = v_out.reshape((len(azimuths),len(radius)))
            # Calculate tangential wind and radial wind
            vt = np.zeros((nLevel,len(azimuths),len(radius)))
            vr = np.zeros((nLevel,len(azimuths),len(radius)))

            for il in range(nLevel):
                for ik in range(len(azimuths)):
                    for ir in range(len(radius)):
                        #vt[il,ik,ir] = v_polar[il,ik,ir]*np.cos(azimuths[ik]*np.pi/180)-u_polar[il,ik,ir]*np.sin(azimuths[ik]*np.pi/180)
                        vt[il,ik,ir] = 0-(u_polar[il,ik,ir]*np.cos(azimuths[ik]*np.pi/180)-v_polar[il,ik,ir]*np.sin(azimuths[ik]*np.pi/180))
                        #vr[il,ik,ir] = u_polar[il,ik,ir]*np.cos(azimuths[ik]*np.pi/180)+v_polar[il,ik,ir]*np.sin(azimuths[ik]*np.pi/180)
                        vr[il,ik,ir] = u_polar[il,ik,ir]*np.sin(azimuths[ik]*np.pi/180)+v_polar[il,ik,ir]*np.cos(azimuths[ik]*np.pi/180)
            # Assemble the dictionary 
            if 'input' in wrf_file:
                if Plot_azmean:
                    d_var['input']= {'radial_wind': np.mean(vr,axis=1), 'tangential_wind': np.mean(vt,axis=1)}
            else:
                if Plot_azmean:
                    d_var['output']= {'radial_wind': np.mean(vr,axis=1), 'tangential_wind': np.mean(vt,axis=1)} 

        else:
            # Read other variables
            tmp = ncdir.variables[var_name][0,:,:,:]
            if 'W' == var_name:
                tmp = (tmp[:-1,:,:]+tmp[1:,:,:])/2
                tmp = tmp.reshape( tmp.shape[0],-1)
            else:
                tmp = tmp.reshape( tmp.shape[0],-1)
            var = tmp[:,idx_x]
            # Interpolate var in geodetic coordinate to polar coordinate
            var_polar = np.zeros((nLevel,len(azimuths),len(radius)))
            for il in range(nLevel):
                Mvar_out = eng.griddata(matlab.double(lon_geo), matlab.double(lat_geo), matlab.double(var[il,:].tolist()), matlab.double(lon_polar), matlab.double(lat_polar) )
                var_out = np.array(Mvar_out._data)
                #var_out = interpolate.griddata(grid_geo,var[il,:],grid_polar,method='linear') # cubic
                var_polar[il,:,:] = var_out.reshape((len(azimuths),len(radius)))
            # Assemble the dictionary 
            if 'input' in wrf_file:
                if Plot_azmean:
                    d_var['input']= {var_name: np.mean(var_polar,axis=1)}
            else:
                if Plot_azmean:
                    d_var['output'] = {var_name: np.mean(var_polar,axis=1)}
    
    # End the matlab process
    eng.quit()

    # Pickle the main dictionary to a file
    if if_save_data:
        # Metadata
        d_data = {'d_var':d_var,'d_coor':d_coor}
        # Save data
        save_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/Model/Azmean_'+var_name+'_'+DAtime+'.pickle'
        with open(save_des, 'wb') as handle:
            pickle.dump(d_data, handle) #, protocol=pickle.HIGHEST_PROTOCOL)
        print( 'Saving the data: ', save_des )

    return None

# ------------------------------------------------------------------------------------------------------
#           Operation: Plot
# ------------------------------------------------------------------------------------------------------
def plot_azmean( var_name,d_var,d_coor  ):

    fig, ax=plt.subplots(1, 3, sharex='all', sharey='all', linewidth=0.5, figsize=(21,10), dpi=400)

    # x axis: radius 
    xv = np.linspace(0,400,100)
    x_axis_rg = range(len(xv))
    f_xinterp = interpolate.interp1d( xv, x_axis_rg)
    # y axis: model vertical coordinate
    if use_pressure:
        y_range = np.arange(900,50,50)
    else:
        y_range = np.arange(0,31,1)
    y_axis_rg = range(len(y_range))
    f_yinterp = interpolate.interp1d( y_range, y_axis_rg)
    yv = f_yinterp( d_coor['ver_coor'] ) 
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
    elif 'hor_wind' in var_name:
        # tangential_wind
        color_intervals = list(np.linspace(0,50,6))
        color_intervals.insert(0,-10.0)
        exist_cmap = plt.cm.jet
        colors = exist_cmap(np.linspace(0,1,len(color_intervals)))
        cmap = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals))
        tan_bounds = color_intervals

        tan_diff_bounds = [0,2,4,6,8,10]
    elif 'W' in var_name:
        bounds = [-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
        cmap = 'jet'
    elif 'Q' in var_name:
        bounds = 10.**np.arange(-6,0,1)
        cmap = 'ocean_r'

    # Plot contourf
    if 'Q' in var_name:
        # Make value==0 as np.nan
        var_input = d_var['input'][var_name]
        #mask1 = var_input <= 1e-10
        #mask2 = var_input == np.nan
        #union_mask = np.logical_or(mask1,mask2)
        #var_input[union_mask] = np.nan
        var_output = d_var['output'][var_name]
        #mask1 = var_output <= 1e-10
        #mask2 = var_output == np.nan
        #union_mask = np.logical_or(mask1,mask2)
        #var_output[union_mask] = np.nan

        ax[0].scatter( xcoor, ycoor, var_input, cmap=cmap,norm=mcolors.LogNorm())
        xa = ax[1].scatter( xcoor, ycoor, var_output, cmap=cmap,norm=mcolors.LogNorm())
        #ax[0].contourf( xcoor, ycoor, var_input, cmap=cmap,levels=bounds,norm=mcolors.LogNorm(),extend='both')
        #xa = ax[1].contourf( xcoor, ycoor, var_output, cmap=cmap,levels=bounds,norm=mcolors.LogNorm(),extend='both')
    elif 'hor_wind' in var_name:
        # customize contour 
        wind_level = [-10,-5,0,5,10]
        colors = [(0, 0, 1), (1, 1, 1), (1, 0, 0)]  # Blue, White, Red
        cmap_name = 'blue_black_red'
        wind_level_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=len(wind_level))

        incre_level = [-1,0,1]
        colors = [(0, 0, 1), (0, 0, 0), (1, 0, 0)]  # Blue, White, Red
        cmap_name = 'blue_b_red'
        incre_level_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=len(incre_level))

        # prior
        ax[0].contourf( xcoor, ycoor, d_var['input']['tangential_wind'], cmap=cmap,levels=tan_bounds,extend='both')
        wind_contour = ax[0].contour(xcoor,ycoor,d_var['input']['radial_wind'],levels=wind_level,cmap=wind_level_cmap,linewidths=2)
        plt.clabel(wind_contour,wind_level,inline=True, fmt="%i", use_clabeltext=True, fontsize=18)
        # posterior
        wind_contourf = ax[1].contourf( xcoor, ycoor, d_var['output']['tangential_wind'], cmap=cmap,levels=tan_bounds,extend='both')
        wind_contour = ax[1].contour(xcoor,ycoor,d_var['output']['radial_wind'],levels=wind_level,cmap=wind_level_cmap,linewidths=2)
        plt.clabel(wind_contour,wind_level,inline=True, fmt="%i", use_clabeltext=True, fontsize=18)
        # analysis increment
        incre_tangential = d_var['output']['tangential_wind']-d_var['input']['tangential_wind']
        incre_radial = d_var['output']['radial_wind']-d_var['input']['radial_wind']
        incre_contourf = ax[2].contourf( xcoor, ycoor, incre_tangential, cmap='ocean_r',levels=tan_diff_bounds,extend='both')
        incre_contour = ax[2].contour(xcoor,ycoor,incre_radial,levels=incre_level,cmap=incre_level_cmap,linewidths=2)
        plt.clabel(incre_contour,incre_level,inline=True, fmt="%i", use_clabeltext=True, fontsize=18)
    else:
        ax[0].contourf( xcoor, ycoor, d_var['input'][var_name], cmap=cmap,levels=bounds,extend='both')
        xa = ax[1].contourf( xcoor, ycoor, d_var['output'][var_name], cmap=cmap,levels=bounds,extend='both')

    # Colorbar
    caxes = fig.add_axes([0.12, 0.03, 0.5, 0.02])
    wind_bar = fig.colorbar(wind_contourf,ax=ax[0:2],orientation="horizontal", cax=caxes)
    wind_bar.ax.tick_params(labelsize=18)

    caxes = fig.add_axes([0.67, 0.03, 0.235, 0.02])
    cbar = fig.colorbar(incre_contourf,ax=ax[2],orientation="horizontal",cax=caxes,extend='both')
    cbar.ax.tick_params(labelsize=18)

    # subplot title and labels
    if use_pressure:
        pass
    else:
        ylabel_like = [0.0,5.0,10.0,15.0,20.0]
        yticks = []
        list_y_range = list(y_range)
        for it in ylabel_like:
            yticks.append( list_y_range.index(it) )
        fig.text(0.06,0.5,'Height (km)',ha='center',va='center',rotation='vertical',fontsize=25)

    for idx in range(3):
        ax.flat[idx].set_ylim(ymin=f_yinterp(0),ymax=f_yinterp(20.5)) # cut off data above 25km 
        ax.flat[idx].set_yticks( yticks )
        ax.flat[idx].set_yticklabels( [str(it) for it in ylabel_like],fontsize=20 )
   
    if 'wind' in var_name:
        xb_title = 'Xb:'+var_name+' (m/s)'
        xa_title = 'Xa:'+var_name+' (m/s)'
    elif 'W' in var_name:
        xb_title = 'Xb: vertical velocity(m/s)'
        xa_title = 'Xa: vertical velocity (m/s)'
    elif 'Q' in var_name:
        xb_title = 'Xb:'+var_name+' (kg/kg)'
        xa_title = 'Xa:'+var_name+' (kg/kg)'

    ax[0].set_title( xb_title, fontsize = 20, fontweight='bold')
    ax[1].set_title( xa_title, fontsize = 20, fontweight='bold')
    ax[2].set_title( 'Xa-Xb', fontsize = 20, fontweight='bold')

    # set X label
    xlabel_like = np.linspace(0,300,7)
    xticks = []
    for it in xlabel_like:
        xticks.append(  f_xinterp( it ) )

    for idx in range(3):
        ax.flat[idx].set_xticks( xticks )
        ax.flat[idx].set_xticklabels(  [str(it) for it in xlabel_like],fontsize=18 )
        ax.flat[idx].set_xlabel('Radius (KM)',fontsize=15)
        ax.flat[idx].set_xlim(xmin=f_xinterp(0),xmax=f_xinterp(300))

    # title
    title_name = Storm+': '+Exper_name+' '+DAtime
    fig.suptitle(title_name, fontsize=15, fontweight='bold')

    # Save figures
    figure_des=plot_dir+DAtime+'_'+var_name+'AZmean.png'
    plt.savefig(figure_des, dpi=400)
    print('Saving the figure: ', figure_des)


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'MARIA'
    MP = 'WSM6'
    DA = 'IR'
    v_interest = ['hor_wind',] #hor_wind

    start_time_str = '201709160000'
    end_time_str = '201709160000'
    Consecutive_times = False

    # model dimension
    xmax = 297
    ymax = 297
    nLevel = 42

    use_pressure = False
    azimuths = np.linspace(0,360,73) #units: degree
    radius = np.linspace(0,300,100) # units: km
    if_calculate = True
    if_save_data = True
    Plot_azmean = True

    # -------------------------------------------------------  
    Exper_name = UD.generate_one_name( Storm,DA,MP )

    # Identify DA times in the period of interest
    if not Consecutive_times:
        DAtimes = ['201709160000',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # create dir to save data
    if if_save_data:
        save_dir = small_dir+Storm+'/'+Exper_name+'/Data_analyze/Model/' 
        savedir_exists = os.path.exists( save_dir )
        if savedir_exists == False:
            os.mkdir(save_dir)

    # create plot_dir
    if Plot_azmean:
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/Az_mean/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

    # if calculate the coordinate transformation
    if if_calculate:
        for DAtime in DAtimes:
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            wrf_files = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
            # convert coordinate
            d_coor = geodetic_to_polar( wrf_files[1] ) #enkf output

            start_time=time.process_time()
            for var_name in v_interest:
                print('Convert '+var_name+' from geodetic to polar coordinate...')
                var_geodetic_to_polar( var_name,d_coor,wrf_files )
            end_time = time.process_time()
            print ('time needed: ', end_time-start_time, ' seconds')

    # if Plot
    if Plot_azmean:
        for DAtime in DAtimes:
            for var_name in v_interest:
                # Read transformed data
                saved_des = small_dir+Storm+'/'+Exper_name+'/Data_analyze/Model/Azmean_'+var_name+'_'+DAtime+'.pickle'
                with open(saved_des, 'rb') as handle:
                    data = pickle.load( handle )
                # Plot
                d_var = data['d_var']
                d_coor = data['d_coor']
                if 'hor_wind' in var_name:
                    plot_azmean( 'hor_wind',d_var,d_coor )
                    #plot_azmean( 'radial_wind',d_var,d_coor )
                    #plot_azmean( 'tangential_wind',d_var,d_coor )
                else:
                    plot_azmean( var_name,d_var,d_coor)
 


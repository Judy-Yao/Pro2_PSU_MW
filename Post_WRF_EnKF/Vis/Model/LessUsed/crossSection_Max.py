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
from matplotlib.cm import get_cmap
from wrf import to_np, getvar, CoordPair, vertcross
import Util_data as UD

def var_cross( var_name,wrf_files ):

    d_var = {}
    d_var['input'] = {}
    d_var['output'] = {}
    
    # ----- Create a cross line -------
    # Find the storm center
    ncdir = nc.Dataset(wrf_files[1], 'r') # enkf output
    slp = UD.compute_slp( ncdir ) # from enkf output
    min_slp = np.min( slp )
    slp_smooth = sp.ndimage.filters.gaussian_filter(slp, [11, 11] )
    idx = np.nanargmin( slp_smooth )
    lat_c = ncdir.variables['XLAT'][:].flatten()[idx]
    lon_c = ncdir.variables['XLONG'][:].flatten()[idx]
    # Create the start point and end point for the cross section
    st_p = CoordPair(lat=lat_c, lon=lon_c-2)
    ed_p = CoordPair(lat=lat_c, lon=lon_c+2)
 
    # loop thro enkf input and output
    for wrf_file in wrf_files:
        print('Reading '+wrf_file)
        ncdir = nc.Dataset(wrf_file, 'r')
    
        # vertical coordinate
        if use_pressure:
            ver_coor = getvar(ncdir,'pres',units="hPa")
        else:
            ver_coor = getvar(ncdir,'z',units="km")

        # Variable of interest
        if var_name == 'U':
            var = getvar(ncdir,'ua',units="m s-1")
        elif var_name == 'V':
            var = getvar(ncdir,'va',units="m s-1")
        elif var_name == 'W':
            var = getvar(ncdir,'wa',units="m s-1")

        # Compute the vertical cross-section interpolation.
        var_cross = vertcross(var,ver_coor,wrfin=ncdir,start_point=st_p,end_point=ed_p,latlon=True,meta=True)
        if 'input' in wrf_file:
            d_var['input']= {var_name: var_cross, 'ver_coor':ver_coor}
        elif 'output' in wrf_file:
            d_var['output']= {var_name: var_cross, 'ver_coor':ver_coor}
    
    # Plot
    if if_plot:
        plot_var_cross( var_name,d_var )



def plot_var_cross( var_name,d_var ):

    # Create the figure
    fig, ax=plt.subplots(1, 2, sharex='all', sharey='all', linewidth=0.5, figsize=(12,6), dpi=400)

    # Set colorbar
    if 'hroi_wind' in var_name:
        bounds = np.linspace(10,70,7)
        cmap = 'jet'
    elif 'W' in var_name:
        bounds = [-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8]
        cmap = 'PiYG'

    ax[0].contourf( to_np(d_var['input'][var_name]),cmap=cmap,levels=bounds,extend='both') # to_np: return the numpy array (value)
    xc = ax[1].contourf( to_np(d_var['output'][var_name]),cmap=cmap,levels=bounds,extend='both')

    # Add the color bar
    cbaxes = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    cbar = fig.colorbar(xc,cax=cbaxes,fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=15) 

    # Set the x-ticks to use latitude and longitude labels.
    coord_pairs = to_np(d_var['output'][var_name].coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str(fmt="{:.1f}, {:.1f}") for pair in to_np(coord_pairs)]
    for i in range(2):
        ax[i].set_xticks(x_ticks[::20])
        ax[i].set_xticklabels(x_labels[::20], rotation=15, fontsize=12)

    # Set the y-ticks to be height.
    vert_val = to_np(d_var['output'][var_name].coords["vertical"])
    vert_vals = [ "{0:.1f}".format(item) for item in vert_val ] 
    v_ticks = np.arange(vert_val.shape[0])
    ax[0].set_yticks(v_ticks[::20])
    ax[0].set_yticklabels(vert_vals[::20], fontsize=15)

    if 'U' in var_name:
        xb_title = 'Xb:'+var_name+' (m/s)'
        xa_title = 'Xa:'+var_name+' (m/s)'
    elif 'W' in var_name:
        xb_title = 'Xb: vertical velocity(m/s)'
        xa_title = 'Xa: vertical velocity (m/s)'
    elif 'Q' in var_name:
        xb_title = 'Xb:'+var_name+' (kg/kg)'
        xa_title = 'Xa:'+var_name+' (kg/kg)'

    ax[0].set_title( xb_title, fontsize = 15, fontweight='bold')
    ax[1].set_title( xa_title, fontsize = 15, fontweight='bold')

    # Set the x-axis and  y-axis labels
    fig.text(0.02,0.5,'Height (km)',ha='center',va='center',rotation='vertical',fontsize=20)
    fig.text(0.5,0.02,"Latitude, Longitude",ha='center',va='center',rotation='horizontal', fontsize=15)

    # title
    title_name = Storm+': '+Exper_name+' '+DAtime
    fig.suptitle(title_name, fontsize=15, fontweight='bold')

    # Save figures
    figure_des=plot_dir+DAtime+'_'+var_name+'_CrossSection.png'
    plt.savefig(figure_des, dpi=400)
    print('Saving the figure: ', figure_des)


if __name__ == '__main__':

    big_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/' #'/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'MARIA'
    MP = 'WSM6'
    DA = 'IR'
    v_interest = ['W',]

    start_time_str = '201709180000'
    end_time_str = '201709180000'
    Consecutive_times = True

    # model dimension
    xmax = 297
    ymax = 297
    nLevel = 42

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
        start_time=time.process_time()
        for var_name in v_interest:
            var_cross( var_name,wrf_files ) 
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')


import os
import glob
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math
from datetime import datetime, timedelta
import time
from netCDF4 import Dataset
from wrf import getvar
import scipy as sp
import scipy.ndimage

import Util_Vis
#import Util_data as UD

Req = 6378137.
Rpol = 6356752.31414
H = 42164160.
ramda_o = -89.5*np.pi/180.

# ------------------------------------------------------------------------
#                    Object: Tbs of Forecasts
# ------------------------------------------------------------------------

# Read crtm calcuated IR data from one binary file
def read_simu_IR_one(Hxb_file, ch_list):

    xmax = 297
    ymax = 297
    print('Hxb_file:' + Hxb_file)

    Hxb_data = np.fromfile(Hxb_file,dtype='<f4') # <: little endian; f: float; 4: 4 bytes
    n_ch = len(Hxb_data)/(xmax*ymax) - 2
    n_ch = int(n_ch)
    if n_ch != len(ch_list):
        print('Error!! # of channels in data is '+str(n_ch))
    Hxb_sim = Hxb_data[:].reshape(n_ch+2,ymax,xmax)

    dict_simu_Tb = {}
    dict_simu_Tb['Lon_x'] = Hxb_sim[0,:,:]
    dict_simu_Tb['Lat_x'] = Hxb_sim[1,:,:]
    dict_simu_Tb['Ch_x'] = ch_list[0]
    dict_simu_Tb['Yb_x'] = Hxb_sim[2,:,:]

    return dict_simu_Tb


def plot_Tb_mspace( Storm, Exper_name, fctime):

    # Read simulated data
    Hxb = big_dir+Storm+'/'+Exper_name+'/'+'CRTM_wrfinput_d03_'+fctime+'.bin'#'CRTM_wrfout_d03_'+fctime+'.bin'
    d_simu = read_simu_IR_one( Hxb, ch_list )

    # Read min slp
    wrf_file =  big_dir+Storm+'/'+Exper_name+'/wrfinput_d03'#'/wrfout_d03_'+fctime
    with Dataset( wrf_file ) as ncid:
        # minimum sea level pressure
        slp = getvar(ncid, 'slp')
        min_slp = np.min( slp )
        # location of the minimum slp
        slp_smooth = sp.ndimage.filters.gaussian_filter(slp, [11, 11] )
        idx = np.nanargmin( slp_smooth )
        lat_storm = ncid.variables['XLAT'][:].flatten()[idx]
        lon_storm = ncid.variables['XLONG'][:].flatten()[idx]

    # ------------------ Plot -----------------------
    f, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,5), dpi=400)

    # Define the domain
    lat_min = np.amin(d_simu['Lat_x'])
    lat_max = np.amax(d_simu['Lat_x'])
    lon_min = np.amin(d_simu['Lon_x'])
    lon_max = np.amax(d_simu['Lon_x'])
  
    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax.scatter(d_simu['Lon_x'], d_simu['Lat_x'],1.4,c=d_simu['Yb_x'],\
                edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())

    # Mark the slp 
    #ax.scatter(lon_storm,lat_storm,20,'blue',marker='*',transform=ccrs.PlateCarree())

    # Colorbar
    caxes = f.add_axes([0.125, 0.05, 0.775, 0.02])
    cbar_ticks = list(range(min_T,max_T,15))
    cbar = f.colorbar(cs, orientation="horizontal",cax=caxes,ticks=cbar_ticks)
    cbar.ax.tick_params(labelsize=10)

    #title for all
    #f.suptitle(Storm+':'+Exper_name+'\n@'+fctime,fontsize=11, fontweight='bold')
    f.suptitle(Storm+' IC:Ref;FC:Ref \n@'+fctime,fontsize=11, fontweight='bold')

    #subplot title
    font = {'size':8,}
    ax.set_title('Forecast--min slp: '+str("{0:.3f}".format(min_slp.values))+' hPa',font,fontweight='bold',fontsize=13)

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for j in range(1):
        gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
       
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
        gl.xlabel_style = {'size': 10}
        gl.ylabel_style = {'size': 10}

    des_name = big_dir+Storm+'/wrfinput_d03_201709151200_GFS_025deg.png'#+Exper_name+'/'+MP+'_'+FCtime+'.png'
    plt.savefig( des_name, dpi=200)
    print('Saving the figure: ', des_name)


# ------------------------------------------------------------------------
#                    Object: Observed Tbs
# ------------------------------------------------------------------------
def read_D_multipleE( Exper_names,FCtime ):

    lat_mins = []
    lat_maxs = []
    lon_mins = []
    lon_maxs = []

    # find two counterparts in expers
    fctime = FCtime[:4]+'-'+FCtime[4:6]+'-'+FCtime[6:8]+'_'+FCtime[8:10]+':00:00'
    wrf_files = []
    for ie in Exper_names:
        #wrf_files.append( big_dir+Storm+'/'+ie+'/wrfout_d03_'+fctime )
        wrf_files.append( big_dir+Storm+'/'+ie+'/wrfinput_d03')

    # Loop thru files
    for ifile in wrf_files:
        ncdir = nc.Dataset(ifile, 'r')
        lat_x = ncdir.variables['XLAT'][0,:,:] 
        lon_x = ncdir.variables['XLONG'][0,:,:] 
        lat_mins.append( np.min( lat_x.flatten() ))
        lat_maxs.append( np.max( lat_x.flatten() ))
        lon_mins.append( np.min( lon_x.flatten() ))
        lon_maxs.append( np.max( lon_x.flatten() ))

    d_wrf_d = {'lat_min':min(lat_mins), 'lat_max':max(lat_maxs), 'lon_min':min(lon_mins), 'lon_max':max(lon_maxs)}
    return d_wrf_d


def read_GOES16(filename, lonlat=True):

    print('Processing '+filename)
    Tb_dict = {}
    ncdir = nc.Dataset( filename )
    # cloud and moisture imagery TB
    tb = ncdir.variables['CMI'][:]
    Tb_dict['Yo'] = tb[::-1,:].flatten()
    if lonlat:
        # GOES fixed grid projection x-coordinate (units: rad)
        x = ncdir.variables['x'][:]
        # GOES fixed grid projection y-coordinate (units: rad) 
        y = ncdir.variables['y'][:]
        xx, yy = np.meshgrid(x,y)
        a = np.sin(xx)**2 + np.cos(xx)**2 * (np.cos(yy)**2 + Req**2/Rpol**2 * np.sin(yy)**2)
        b = - 2.* H * np.cos(xx) * np.cos(yy)
        c = H**2 - Req**2
        Rs = (-b - np.sqrt(b**2 - 4.*a*c)) / (2.*a)
        Sx = Rs * np.cos(xx) * np.cos(yy)
        Sy = -Rs * np.sin(xx)
        Sz = Rs * np.cos(xx) * np.sin(yy)
        lats = np.arctan( Req**2/Rpol**2 * Sz / np.sqrt((H - Sx)**2 + Sy**2))
        lons = ramda_o - np.arctan(Sy / (H - Sx))
        Tb_dict['lon'] = lons[::-1,:].flatten()/np.pi*180.
        Tb_dict['lat'] = lats[::-1,:].flatten()/np.pi*180.

    return Tb_dict


def plot_Tb_obs( FCtime ):

    # convert from date to the day of the year
    FCtime_dt = datetime.strptime(FCtime[:8],"%Y%m%d")
    doy = FCtime_dt.timetuple().tm_yday # day of the year
    print(doy)
    # Read observed Tb file
    #Yo_dir = small_dir+Storm+'/'+Exper_name+'/obs_Tb/'
    Yo_dir = small_dir+Storm+'/Obs_y/'
    print(str(doy))
    Yo_file = glob.glob(Yo_dir+'OR_ABI-L2-CMIPF-M3C08_G16_s'+FCtime[:4]+str(doy)+FCtime[8:10]+'*')
    d_tb = read_GOES16( Yo_file[0], lonlat=True )

    # Find the domain range
    domain = read_D_multipleE( Exper_names, FCtime )
   
    # Find the area 
    con1 = np.where((d_tb['lat']<=domain['lat_max'])&(d_tb['lat'] >= domain['lat_min']))[0]
    con2 = np.where((d_tb['lon'] <= domain['lon_max'])&(d_tb['lon'] >= domain['lon_min']))[0]
    idx = list(set(con1).intersection(set(con2)))
   
    # ------------------ Plot -----------------------
    f, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,5), dpi=400)

    # Define the domain
    lat_min = domain['lat_min']
    lat_max = domain['lat_max']
    lon_min = domain['lon_min']
    lon_max = domain['lon_max']

    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax.scatter(d_tb['lon'][idx], d_tb['lat'][idx],1,c=d_tb['Yo'][idx],\
                edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())

    # Colorbar
    caxes = f.add_axes([0.125, 0.05, 0.775, 0.02])
    cbar_ticks = list(range(min_T,max_T,15)) 
    cbar = f.colorbar(cs, orientation="horizontal",cax=caxes,ticks=cbar_ticks)
    cbar.ax.tick_params(labelsize=10)

    #title for all
    f.suptitle(Storm+': Observed Tb \n@'+FCtime,fontsize=11, fontweight='bold')

    #subplot title
    font = {'size':8,}
    #ax.set_title('Forecast--min slp: '+str("{0:.3f}".format(min_slp.values))+' hPa',font,fontweight='bold',fontsize=13)

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
    gl.top_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.right_labels = False
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}

    #des_name = small_dir+Storm+'/free_run_MPs/obs_Tb/'+'obs_'+FCtime+'.png'
    des_name = 'obs_'+FCtime+'.png'
    plt.savefig( des_name, dpi=200)
    print('Saving the figure: ', des_name)


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'MARIA'
    Exper_name = 'IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900/fc/201709151200'#'IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi300/wrf_df/201708240600'
    #Exper_names = ['IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900/fc/201709151200',]
    #Exper_names = ['IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi300/wrf_df/201708240600/','JerryRun/IR_THO/wrf_df/201708240600/','JerryRun/IR_WSM6/wrf_df/201708240600/']
    #['IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900/wrf_df/201709041800','IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900/wrf_df/201709041800','IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900_origin/wrf_df/201709041800']
    MP = 'WSM6'
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']
    num_ens = 60

    start_time_str = '201709151200'
    end_time_str = '201709151200'
    Consecutive_times = True

    If_plot_obs = False
    If_plot = True
    # -------------------------------------------------------   

    if not Consecutive_times:
        IR_times = ['201709050000','201709050600','201709051200','201709051800','201709060000',]
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 3))]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
   
    # Plot Tbs of forecast
    if If_plot:
        for FCtime in IR_times:
            print('FCtime: '+ FCtime)
            fctime = FCtime[:4]+'-'+FCtime[4:6]+'-'+FCtime[6:8]+'_'+FCtime[8:10]+':00:00'
            plot_Tb_mspace( Storm, Exper_name, fctime,) 
        
    # Plot Tbs of observation
    if If_plot_obs:
        for FCtime in IR_times:
            print('FCtime: '+ FCtime)
            plot_Tb_obs( FCtime ) 



#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os
import glob
import numpy as np
import Util_Vis
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from global_land_mask import globe
import math
from datetime import datetime, timedelta
import matlab.engine



# read thinned IR data from SO files
def read_obs(obs_file):

    with open(obs_file) as tmp:
        obs_all = tmp.readlines()

    Ss_obs = []
    Ch_obs = []
    Lat_obs = []
    Lon_obs = []
    Yo_obs = []
    for line in obs_all:
        line_split = line.split()
      
        Ss_obs.append(line_split[1])
        Ch_obs.append(line_split[2])
        Lat_obs.append(float(line_split[3]))
        Lon_obs.append(float(line_split[4]))
        Yo_obs.append(float(line_split[5]))

    uni_ss_obs = list(set(Ss_obs))
    uni_ch_obs = list(set(Ch_obs))

    print( 'Sensor: '+ uni_ss_obs[0] )
    print( 'Channel number: '+ uni_ch_obs[0] )

    Lat_obs = np.array(Lat_obs)
    Lon_obs = np.array(Lon_obs)
    Yo_obs = np.array(Yo_obs)
    # Note: If you want to build up your matrix one column at a time,
    # you might be best off to keep it in a list until it is finished, and only then convert it into an array.    

    dict_obs_all = {'Lat_obs': Lat_obs, 'Lon_obs': Lon_obs, 'Yo_obs': Yo_obs}
    return dict_obs_all

def read_TCvitals(tc_file, DAtime):

    with open(tc_file) as tmp:
        tc_all = tmp.readlines()

    tc_lat = []
    tc_lon = []
    for line in tc_all:
        line_split = line.split()
        tc_time = line_split[3]+line_split[4]

        if tc_time == DAtime:
            print('Time from TCvitals:', tc_time)
            # Read latitude
            if 'N' in line_split[5]:
                tc_lat.append(float(line_split[5].replace('N',''))/10)
            else:
                tc_lat.append( 0-float(line_split[5].replace('S',''))/10)
            # Read longitude
            if 'W' in line_split[6]:
                tc_lon.append(0-float(line_split[6].replace('W',''))/10)
            else:
                tc_lon.append(float(line_split[6].replace('E',''))/10)
            
            break

    return tc_lon, tc_lat


# Read crtm calcuated IR data from binary file
def read_simu_Tb(Hxb_file, Hxa_file, ch_list):
    
    xmax = 297
    ymax = 297
    print('Hxb_file:' + Hxb_file[0])
    print('Hxa_file:' + Hxa_file[0])
    print('xmax, ymax: '+str(xmax)+' '+str(ymax))

    Hxb_data = np.fromfile(Hxb_file[0],dtype='<f4') # <: little endian; f: float; 4: 4 bytes
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

    Hxa_data = np.fromfile(Hxa_file[0],dtype='<f4')
    n_ch = len(Hxa_data)/(xmax*ymax) - 2
    n_ch = int(n_ch)
    if n_ch != len(ch_list):
        print('Error!! # of channels in data is '+str(n_ch))
    Hxa_sim = Hxa_data[:].reshape(n_ch+2,ymax,xmax)
    dict_simu_Tb['Ya_x'] = Hxa_sim[2,:,:] 
    #for rec in range(n_ch):
    #    dict_simu_Tb[ch_list[rec]] = sim[rec+2,:,:]
    
    return dict_simu_Tb


def plot_Tb(Storm, Exper_name, DAtime):

    # Read Obs data
    obs_file_name = 'radiance_d03_' + DAtime + '_so'
    d_obs = read_obs('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/Obs_y/IR/'+obs_file_name)
    
    # Read simulated data
    ch_list = ['8',]
    Hx_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime
    Hxb_file = [Hx_dir+'/wrf_enkf_input_d03_mean_'+DAtime+'_tb_g16_crtm.bin',]
    Hxa_file = [Hx_dir+'/wrf_enkf_output_d03_mean_'+DAtime+'_tb_g16_crtm.bin',]

    d_simu = read_simu_Tb( Hxb_file, Hxa_file, ch_list )
    # interpolate simulated Tbs to obs physical space
    mYb_obspace = eng.griddata(matlab.double(d_simu['Lon_x'].tolist()), matlab.double(d_simu['Lat_x'].tolist()), matlab.double(d_simu['Yb_x'].tolist()), matlab.double(d_obs['Lon_obs'].tolist()), matlab.double(d_obs['Lat_obs'].tolist()))
    Yb_obspace = np.array(mYb_obspace._data)
    mYa_obspace = eng.griddata(matlab.double(d_simu['Lon_x'].tolist()), matlab.double(d_simu['Lat_x'].tolist()), matlab.double(d_simu['Ya_x'].tolist()), matlab.double(d_obs['Lon_obs'].tolist()), matlab.double(d_obs['Lat_obs'].tolist()))
    Ya_obspace = np.array(mYa_obspace._data)

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = read_TCvitals('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # ------------------ Plot -----------------------
    f, ax=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,2.5), dpi=400)

    # Define the domain
    lat_min = np.amin(d_obs['Lat_obs'].flatten())
    lat_max = np.amax(d_obs['Lat_obs'].flatten())
    lon_min = np.amin(d_obs['Lon_obs'].flatten())
    lon_max = np.amax(d_obs['Lon_obs'].flatten())
  
    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )

    ax[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    ax[0].scatter(d_obs['Lon_obs'],d_obs['Lat_obs'],1.5,c=d_obs['Yo_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())

    ax[1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    ax[1].scatter(d_obs['Lon_obs'], d_obs['Lat_obs'],1.5,c=Yb_obspace,\
                edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax[1].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())

    ax[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax[2].scatter(d_obs['Lon_obs'], d_obs['Lat_obs'],1.5,c=Ya_obspace,\
                edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax[2].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())

    # Colorbar
    caxes = f.add_axes([0.2, 0.97, 0.6, 0.02])
    cbar = f.colorbar(cs, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)

    #subplot title
    font = {'size':8,}
    ax[0].set_title('Yo', font, fontweight='bold')
    ax[1].set_title('H(Xb)', font, fontweight='bold')
    ax[2].set_title('H(Xa)', font, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2))

    for j in range(3):
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
        gl.xlabel_style = {'size': 4}
        gl.ylabel_style = {'size': 6}

    plt.savefig('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IR/'+DAtime+'_g16_ch8_Obspace.png', dpi=300)
 




if __name__ == '__main__':
    Storm = 'MARIA'
    Exper_name = 'newWRF_IR_only'

    start_time_str = '201709171700' 
    end_time_str = '201709171800'
    
    time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
    time_diff_hour = time_diff.total_seconds() / 3600
    time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
    IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    eng = matlab.engine.start_matlab() # start a new matlab process
    for DAtime in IR_times:
        print('DAtime: '+ DAtime)
        plot_Tb( Storm, Exper_name, DAtime)
       
    eng.quit()
    
















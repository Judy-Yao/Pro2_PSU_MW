#!/usr/bin/env python3

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

def read_obs(obs_file):

    with open(obs_file) as tmp:
        obs_all = tmp.readlines()

    Ch_obs = []
    Lat_obs = []
    Lon_obs = []
    Yo_obs = []
    for line in obs_all:
        line_split = line.split() 
        Ch_obs.append(int(line_split[2]))
        Lat_obs.append(float(line_split[3]))
        Lon_obs.append(float(line_split[4]))
        Yo_obs.append(float(line_split[5]))

    Ch_obs = np.array(Ch_obs)
    Lat_obs = np.array(Lat_obs)
    Lon_obs = np.array(Lon_obs)
    Yo_obs = np.array(Yo_obs)
    # Note: If you want to build up your matrix one column at a time,
    # you might be best off to keep it in a list until it is finished, and only then convert it into an array.    

    dict_obs_all = {'Ch_obs': Ch_obs, 'Lat_obs': Lat_obs, 'Lon_obs': Lon_obs, 'Yo_obs': Yo_obs} 
    return dict_obs_all


def getSensor_Ch(obs_file):

    # Read the content inside the microwave obs to a list
    with open(obs_file) as f:
        all_lines = f.readlines() 

    # Declare an empty list 
    sensorCh_rep = [] 
    sensor_rep = []

    # For each string (elements in this list), separate it into words
    for line in all_lines:
        split_all_lines = line.split() 
        sensor_rep.append(split_all_lines[1])
        # Combine sensor name and channel number together 
        sensorCh_rep.append(split_all_lines[1] + ' ' + split_all_lines[2])

    # Find the uqiue sensor / sensor-and-channel combination
    sensor_uni_set = set(sensor_rep)
    sensor_uni = list(sensor_uni_set)

    sensorCh_uni_set = set(sensorCh_rep)
    sensorCh_uni = list(sensorCh_uni_set)

    # For each combination, separate it into sensor name and channel number
    sensor_Ch = []
    for record in sensorCh_uni: 
        sensor_Ch.append(record.split())

    Ch_perSS = []
    # Loop through each unique sensor
    for iss in sensor_uni:
        for ir in sensor_Ch:      
            Ch_perSS.append([]) 
            if sensor_Ch[sensor_Ch.index(ir)][0] == iss:
                Ch_perSS[sensor_uni.index(iss)].append(sensor_Ch[sensor_Ch.index(ir)][1])
            else:
                continue
            
    # Build a dictionary: sensor = channel1, channel2, ...
    dict_ss_ch = {}
    for iss in sensor_uni:
        dict_ss_ch[sensor_uni[sensor_uni.index(iss)]] = Ch_perSS[sensor_uni.index(iss)]

    return dict_ss_ch 


def read_simu_Tb(Hxb_file, Hxa_file):

    ncdir = nc.Dataset(Hxb_file, 'r')

    Lat_x = ncdir.variables['XLAT'][0,:,:] #latitude: XLAT(time, y, x)
    Lon_x = ncdir.variables['XLONG'][0,:,:] #longitude: XLONG(time, y, x)
    Ch_x = ncdir.variables['ch'][:] # int ch(ch)
    Yb_x  = ncdir.variables['Brightness_Temperature'][0,:,:,:] #BackgroundTB:float Brightness_Temperature(time, ch, y, x)    

    ncdir = nc.Dataset(Hxa_file, 'r') 
    Ya_x = ncdir.variables['Brightness_Temperature'][0,:,:,:]
    # Comment: in matlab, data-read order is reversed on the data-storage order.
    
    dict_simu_Tb = {'Ch_x': Ch_x, 'Lat_x': Lat_x, 'Lon_x': Lon_x, 'Yb_x': Yb_x, 'Ya_x': Ya_x}
    return dict_simu_Tb

    
def plot_Tb(Storm, Exper_name, DAtime, sensor):
   
    ch_num = [int(ich) for ich in dict_ss_ch[sensor]]
    # Define the low frequency for each sensor
    d_lowf = {'amsr2':7, 'gmi_gpm_lf':3, 'ssmi': 13}
    # categorize ssmi and ssmis as one kind 'ssmi' since they share the same channel number set up
    if 'ssmi' in sensor:
        sensor_short = 'ssmi'
    else:
        sensor_short = sensor

    # Read data
    obs_file_name = 'microwave_d03_' + DAtime + '_so'
    d_obs = read_obs('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/Obs_y/MW/'+obs_file_name)
   
    Hx_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Obs_Hx/MW/'+DAtime
    Hxb_file = glob.glob( Hx_dir+'/wrf_enkf_input_d03_mean*nc' )
    Hxa_file = glob.glob( Hx_dir+'/wrf_enkf_output_d03_mean*nc')  
    d_simu = read_simu_Tb( Hxb_file[0], Hxa_file[0] )

    # ------------------ Plot -----------------------
    f, ax=plt.subplots(2, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  dpi=400)
    
    # Customize colormap
    max_T=300
    min_T=80
    min_Jet=150
    MWJet = Util_Vis.newJet(max_T=300, min_T=80, min_Jet=150)

    #scatter_size = [2.5, 2.5]
    # Define the domain
    lat_min = np.amin(d_simu['Lat_x'].flatten())
    lat_max = np.amax(d_simu['Lat_x'].flatten())
    lon_min = np.amin(d_simu['Lon_x'].flatten())
    lon_max = np.amax(d_simu['Lon_x'].flatten())

    # Loop over low frequency and high frequency if available
    for irow in range(2):
        # make sure the idex is not out of range    
        if len(ch_num) == 2:
            i = irow
        else:
            if ch_num[0] == d_lowf[sensor_short]:
                i = 0
            else:
                i =  1

        # ---------- obs ---------------------
        if len(ch_num) == 2:
            ch_idx = d_obs['Ch_obs'] == ch_num[i]
        else:
            ch_idx = d_obs['Ch_obs'] == ch_num

        Lat_obs_ch = d_obs['Lat_obs'][ch_idx] 
        Lon_obs_ch = d_obs['Lon_obs'][ch_idx]
        Yo_obs_ch = d_obs['Yo_obs'][ch_idx]

        if d_obs['Ch_obs'][ch_idx].any() == d_lowf[sensor_short]:
            is_ocean = globe.is_ocean(Lat_obs_ch, Lon_obs_ch)
            mask_x = is_ocean
        else:
            mask_x = np.full((np.size(Lon_obs_ch), ), True)
   
        ax[i,0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i,0].coastlines(resolution='10m', color='black',linewidth=0.5)
        ax[i,0].scatter(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],2.5,c=Yo_obs_ch[mask_x],edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
   
        # ------------ simulated Tb ----------------------
        if len(ch_num) == 2:
            for ich in d_simu['Ch_x'].tolist():
                if ich == ch_num[i]:
                    Ch_idx = d_simu['Ch_x'].tolist().index(ich)
        else:
            for ich in d_simu['Ch_x'].tolist():
                if ich == ch_num:
                    Ch_idx = d_simu['Ch_x'].tolist().index(ich)      

        Yb_x_ch = d_simu['Yb_x'][Ch_idx,:,:] 
        Ya_x_ch = d_simu['Ya_x'][Ch_idx,:,:]
   
        if d_simu['Ch_x'][Ch_idx] == d_lowf[sensor_short]:
            is_ocean = globe.is_ocean(d_simu['Lat_x'].flatten(), d_simu['Lon_x'].flatten())  
            mask_x = is_ocean
        else:
            mask_x = np.full((np.size(d_simu['Lon_x'].flatten()), ), True)
      
        ax[i,1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax[i,1].coastlines(resolution='10m', color='black',linewidth=0.5)
        ax[i,1].scatter(d_simu['Lon_x'].flatten()[mask_x], d_simu['Lat_x'].flatten()[mask_x],2.5,c=Yb_x_ch.flatten()[mask_x],\
            edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree()) 

        ax[i,2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax[i,2].coastlines(resolution='10m', color='black',linewidth=0.5)
        cs = ax[i,2].scatter(d_simu['Lon_x'].flatten()[mask_x], d_simu['Lat_x'].flatten()[mask_x],2.5,c=Ya_x_ch.flatten()[mask_x],\
            edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())

#plt.xlim([-91,-84])
#plt.ylim([15,24])

    # Colorbar
    caxes = f.add_axes([0.2, 0.97, 0.6, 0.02])
    cbar = f.colorbar(cs, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)
    #plt.text( 0.8, 0.7, 'Brightness Temperature (K)', fontsize=6, transform=transAxes)

    #subplot title
    font = {'size':8,}
    ax[0,0].set_title('Yo', font, fontweight='bold')
    ax[0,1].set_title('HXb', font, fontweight='bold')
    ax[0,2].set_title('HXa (IR+MW)', font, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2)) 
    for i in range(2):
        for j in range(3):
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                        linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
            # Control if ticks are added to a certain side
            if i==0:
                gl.bottom_labels = False
                gl.top_labels = False
            else:
                gl.bottom_labels = True
                gl.top_labels = False        

            if j==0:
                gl.left_labels = True
                gl.right_labels = False
            else:
                gl.left_labels = False
                gl.right_labels = False  
        # Add ticks
        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 4}
        gl.ylabel_style = {'size': 5} 
    
    plt.savefig('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Vis_analyze/Tb/MW/'+DAtime+'_'+sensor+'.png', dpi=300)
    

if __name__ == '__main__':
    Storm = 'HARVEY'
    Exper_name = 'MW_THO'
    MW_times = ['201708221200','201708221300']
    for DAtime in MW_times:
        obs_file_name = 'microwave_d03_' + DAtime + '_so'
        dict_ss_ch = getSensor_Ch( '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/HARVEY/Obs_y/MW/'+obs_file_name )
        
        for sensor in dict_ss_ch:
            plot_Tb( Storm, Exper_name, DAtime, sensor )



    # Path
    #F_Obs = 'microwave_d03_201708221200_so.txt'
    #F_Hxb_mean = 'wrf_enkf_input_d03_mean.tb.ssmis_f17.crtm.conv.txt'
    #F_Hxa_mean = 'wrf_enkf_output_d03_mean.tb.ssmis_f17.crtm.conv.txt' 

    #SSMIS(F_Obs, F_Hxb_mean, F_Hxa_mean)

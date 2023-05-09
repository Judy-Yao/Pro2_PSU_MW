#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import pickle
import netCDF4 as nc
import Diagnostics as Diag
import math
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import Util_Vis

def read_wrf_domain( wrf_file ):

    print('Read domain info from: ' + wrf_file)
    ncdir = nc.Dataset(wrf_file, 'r')

    Lat_x = ncdir.variables['XLAT'][0,:,:] #latitude: XLAT(time, y, x)
    Lon_x = ncdir.variables['XLONG'][0,:,:] #longitude: XLONG(time, y, x)

    lat_min = np.min( Lat_x.flatten() )
    lat_max = np.max( Lat_x.flatten() )
    lon_min = np.min( Lon_x.flatten() )
    lon_max = np.max( Lon_x.flatten() )

    d03_list = {'lat_min':lat_min, 'lat_max':lat_max, 'lon_min':lon_min, 'lon_max':lon_max}
    return d03_list


def getSensor_Ch( MW_diag ):

    # Declare an empty list 
    sensorCh_rep = []
    sensor_rep = []

    # For each string (elements in this list), separate it into words
    for irec in range(len(MW_diag[0])):
        sensor_rep.append( MW_diag[0][irec])
        # Combine sensor name and channel number together 
        sensorCh_rep.append( MW_diag[0][irec] + ' ' + MW_diag[1][irec] )

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

    # Special treatment to gmi_gpm sensor
    gmi_ch = [key_sensor for key_sensor in dict_ss_ch if 'gmi_gpm' in key_sensor]
    if len( gmi_ch ) >= 1:
        rem_key = ['gmi_gpm_lf','gmi_gpm_hf']
        add_ss = 'gmi_gpm'
        add_ch_num = [dict_ss_ch[ikey][0] for ikey in rem_key]
        [dict_ss_ch.pop(ikey) for ikey in rem_key]
        dict_ss_ch[add_ss] = add_ch_num
        return dict_ss_ch
    else:
        return dict_ss_ch



def Plot_MW( Storm, DAtime, Exper, MW_diag_Expers ):

    sensor = 'ssmis_f17'
    dict_ss_ch = getSensor_Ch( MW_diag_Expers[Exper[0]] )
    ch_num = [int(ich) for ich in dict_ss_ch[sensor]]
    print( ch_num )
    # Define the low and high frequency for each sensor
    d_lowf = {'atms_npp':0, 'amsr2_gcom-w1':7, 'gmi_gpm':3, 'mhs_n19':0, 'mhs_n18':0, 'mhs_metop-a':0, 'mhs_metop-b':0, 'saphir_meghat':0, 'ssmis_f16': 13, 'ssmis_f17': 13, 'ssmis_f18': 13, 'ssmi_f15':1}

    # Read WRF domain
    wrf_file = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/JerryRun/MW_THO/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = read_wrf_domain( wrf_file )

    # ------------------ Plot -----------------------
    f, ax=plt.subplots(2, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(6,4), dpi=400)

    # Customize colormap
    max_T=300
    min_T=80
    min_Jet=150
    MWJet = Util_Vis.newJet(max_T=300, min_T=80, min_Jet=150)

    #scatter_size = [2.5, 2.5]
    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    # Loop over low frequency and high frequency if available
    for input_it in range(2):
        # make sure the idex is not out of range    
        if len(ch_num) == 2:
            if ch_num[input_it] == d_lowf[sensor]:
                i = 0
            else:
                i = 1
        else:
            if ch_num[0] == d_lowf[sensor]:
                i = 0
            else:
                i =  1
        
        ch_obs = [ int(it) for it in MW_diag_Expers[Exper[0]][1] ]
        if len(ch_num) == 2:
            ch_idx_lg = [ int(it) == ch_num[input_it]  for it in MW_diag_Expers[Exper[0]][1] ]
        else:
            ch_idx_lg = [ int(it) == ch_num[0] for it in MW_diag_Expers[Exper[0]][1] ]
        ch_idx = np.where( ch_idx_lg)[0] 

        Lat_obs_ch = [ MW_diag_Expers[Exper[0]][2][idx] for idx in ch_idx ]
        Lon_obs_ch = [ MW_diag_Expers[Exper[0]][3][idx] for idx in ch_idx ]
        Yo_obs_ch =[  MW_diag_Expers[Exper[0]][4][idx] for idx in ch_idx ] 

        # Obs
        ax[i,0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i,0].coastlines(resolution='10m', color='black',linewidth=0.5)
        cs = ax[i,0].scatter(Lon_obs_ch, Lat_obs_ch,2.5,c=Yo_obs_ch,\
                 edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())

        Yb_0 = [  MW_diag_Expers[Exper[0]][5][idx] for idx in ch_idx ]
        Yb_1 = [ MW_diag_Expers[Exper[1]][5][idx] for idx in ch_idx ]
        Yb_diff = [Yb_0[ir] - Yb_1[ir] for ir in range(len(Yb_0)) ]
        idx_b = np.where([it != 0 for it in Yb_diff])[0]
        Lon_Yb = [Lon_obs_ch[it] for it in idx_b ]
        Lat_Yb = [Lat_obs_ch[it] for it in idx_b ]
        #Yb_idxb = [Yb_diff[it] for it in idx_b ]
        # Hxb
        ax[i,1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax[i,1].coastlines(resolution='10m', color='black',linewidth=0.5)
        ax[i,1].scatter(Lon_Yb, Lat_Yb,2.5,'red',\
        #ax[i,1].scatter(Lon_obs_ch, Lat_obs_ch,2.5,Yb_diff,\
                  edgecolors='none', cmap='RdBu_r', vmin=-0.5, vmax=0.5, transform=ccrs.PlateCarree())

        Ya_0 = [ MW_diag_Expers[Exper[0]][6][idx] for idx in ch_idx ]
        Ya_1 = [ MW_diag_Expers[Exper[1]][6][idx] for idx in ch_idx ]
        Ya_diff = [Ya_0[ir] - Ya_1[ir] for ir in range(len(Ya_0)) ]
        idx_a = np.where([it != 0 for it in Ya_diff])[0]
        Lon_Ya = [Lon_obs_ch[it] for it in idx_a ]
        Lat_Ya = [Lat_obs_ch[it] for it in idx_a ]
        Ya_idxa = [Ya_diff[it] for it in idx_a ]
        # Hxa
        ax[i,2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax[i,2].coastlines(resolution='10m', color='black',linewidth=0.5)
        ax[i,2].scatter(Lon_Ya, Lat_Ya,2.5,'red',\
        #ax[i,2].scatter(Lon_obs_ch, Lat_obs_ch,2.5,Ya_diff,\
                edgecolors='none', cmap='RdBu_r', vmin=-0.5, vmax=0.5, transform=ccrs.PlateCarree())


    # Colorbar
    caxes = f.add_axes([0.2, 0.06, 0.6, 0.02])
    #cbar = f.colorbar(cs, orientation="horizontal", cax=caxes, ticks=[-0.5, 0, 0.5])
    cbar = f.colorbar(cs, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)
    #plt.text( 0.8, 0.7, 'Brightness Temperature (K)', fontsize=6, transform=transAxes)

    #subplot title
    font = {'size':8,}
    ax[0,0].set_title('Yo Assimilated', font, fontweight='bold')
    ax[0,1].set_title('H(Xb) Diverged', font, fontweight='bold')
    ax[0,2].set_title('H(Xa) Diverged', font, fontweight='bold')

    f.suptitle(Exper[0]+' v.s. '+Exper[1] , fontsize=10)
    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for i in range(2):
        for j in range(3):
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
            # Control if ticks are added to a certain side
            if i==0:
                gl.xlabels_bottom = False
                gl.xlabels_top = False
            else:
                gl.xlabels_bottom = True
                gl.xlabels_top = False

            if j==0:
                gl.ylabels_left = True
                gl.ylabels_right = False
            else:
                gl.ylabels_left = False
                gl.ylabels_right = False
            # Add ticks
            gl.ylocator = mticker.FixedLocator(lat_ticks)
            gl.xlocator = mticker.FixedLocator(lon_ticks)
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 4}
            gl.ylabel_style = {'size': 6}
    
    save_path = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Exper[1]+'/Vis_analyze/Enkf/label_microwave_201708221200.png'
    plt.savefig(save_path, dpi=300)
    print(save_path)


def Plot_IR(Storm, DAtime, Exper, IR_diag_Expers):

    # Read WRF domain
    wrf_file = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/JerryRun/MW_THO/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = read_wrf_domain( wrf_file )
    # ------------------ Plot -----------------------
    f, ax=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,2.5), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    #Define Tb threshold
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )

    Lat_obs_ch = IR_diag_Expers[Exper[0]][0]
    Lon_obs_ch = IR_diag_Expers[Exper[0]][1]
    Yo_obs_ch =  IR_diag_Expers[Exper[0]][2]
    ax[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax[0].scatter(Lon_obs_ch,Lat_obs_ch,1.5,c=Yo_obs_ch,edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())

    Yb_0 = IR_diag_Expers[Exper[0]][3]
    Yb_1 = IR_diag_Expers[Exper[1]][3]
    Yb_diff = [Yb_0[ir] - Yb_1[ir] for ir in range(len(Yb_0)) ]
    idx_b = list(np.where([it != 0 for it in Yb_diff])[0])
    Lon_Yb = [Lon_obs_ch[it] for it in idx_b ]
    Lat_Yb = [Lat_obs_ch[it] for it in idx_b ]
    ax[1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    ax[1].scatter(Lon_Yb, Lat_Yb,1.5,'red',\
                edgecolors='none', cmap='RdBu_r', vmin=-10, vmax=10, transform=ccrs.PlateCarree())

    Ya_0 = IR_diag_Expers[Exper[0]][4]
    Ya_1 = IR_diag_Expers[Exper[1]][4]
    Ya_diff = [Ya_0[ir] - Ya_1[ir] for ir in range(len(Ya_0)) ]
    idx_a = np.where([it != 0 for it in Ya_diff])[0]
    Lon_Ya = [Lon_obs_ch[it] for it in idx_a ]
    Lat_Ya = [Lat_obs_ch[it] for it in idx_a ]
    ax[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    ax[2].scatter(Lon_Ya, Lat_Ya,1.5,'red',\
    #ax[2].scatter(Lon_obs_ch, Lat_obs_ch,1.5,Ya_diff,\
                edgecolors='none', cmap='RdBu_r', vmin=-0.5, vmax=0.5, transform=ccrs.PlateCarree())

    # Colorbar
    caxes = f.add_axes([0.2, 0.1, 0.6, 0.02])
    #cbar = f.colorbar(cs, orientation="horizontal", cax=caxes, ticks=[-0.5, 0, 0.5])
    cbar = f.colorbar(cs, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)

    #subplot title
    font = {'size':8,}
    ax[0].set_title('Yo Assimilated', font, fontweight='bold')
    ax[1].set_title('H(Xb) Diverged', font, fontweight='bold')
    ax[2].set_title('H(Xa) Diverged', font, fontweight='bold')

    f.suptitle(Exper[0]+' v.s. '+Exper[1] , fontsize=10)
    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

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

    save_path = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Exper[1]+'/Vis_analyze/Enkf/label_radiance_201708221200.png'
    plt.savefig(save_path, dpi=300)
    print(save_path)

def Diagnose(big_dir, Storm, Exper, time, v_interest, Check_IR, Check_MW):
   
    if Check_MW:
        # Read in MW diagnostics with labelled sensor info for each experiment
        MW_diag_Expers = {}
        for iExper in Exper:
            file_Diag = big_dir+Storm+'/'+iExper+'/run/'+time+'/enkf/d03/fort.10000'
            MW_so = big_dir+Storm+'/'+iExper+'/run/'+time+'/enkf/d03//microwave_201708221200_so' 
            MW_diag_Expers[iExper] = Diag.label_mw_obs( file_Diag, MW_so, v_interest )
    
        Plot_MW( Storm, time, Exper, MW_diag_Expers )
  
    if Check_IR:
        # Read in MW diagnostics with labelled sensor info for each experiment
        IR_diag_Expers = {}
        for iExper in Exper:
            file_Diag = big_dir+Storm+'/'+iExper+'/run/'+time+'/enkf/d03/fort.10000'
            IR_diag_Expers[iExper] = Diag.Find_IR( file_Diag, v_interest )
        
        print(IR_diag_Expers[Exper[0]][3] == IR_diag_Expers[Exper[1]][3])

        Plot_IR( Storm, time, Exper, IR_diag_Expers )
  

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'

    # configuration
    Storm = 'HARVEY'
    Exper_name = ['JerryRun/MW_THO','IR+MW-J_DA+J_WRF+J_init-SP-intel17']
    time = '201708221200'
    #filename = '/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/IR+MW-J_DA+J_WRF+J_init-SP-intel19/run/201708221200/enkf/d03/fort.10000'
    v_interest = ['obs_type','lat','lon','height','i_grid','j_grid','k_grid','Hroi','Vroi','obs','obs_error_variance','prior_mean','posterior_mean','prior_spread','posterior_spread'] #['obs_type','lat','lon','obs','prior_mean','posterior_mean']
    Check_MW = True
    Check_IR = True

    Diagnose(big_dir, Storm, Exper_name, time, v_interest, Check_IR, Check_MW)

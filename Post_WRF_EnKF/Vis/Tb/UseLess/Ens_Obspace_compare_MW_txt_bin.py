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
import time

import Read_Obspace_MW as ROMW
import Util_data as UD
import Util_Vis
import Diagnostics as Diag

def RMSE(simu, obs):
    return np.sqrt( ((simu - obs) ** 2).mean() )

def Bias(simu, obs):
    return  np.sum((simu - obs),0)/np.size(obs,0)

def read_memTb( Hx_dir, sensor, DAtime, mem ):

    # number of columns of each record
    ncol = 8
    d_MW = {}
    mem_id = '{:03d}'.format(mem+1)
    input_file = Hx_dir+"/input_mem"+mem_id+'_d03_'+DAtime+'.tb.'+sensor+'.crtm.conv.txt'
    output_file = Hx_dir+"/output_mem"+mem_id+'_d03_'+DAtime+'.tb.'+sensor+'.crtm.conv.txt'

    tmp_control = np.fromfile( input_file, sep=' ' ) # read input file
    lat_obs = [float(item) for item in tmp_control[0::ncol]]
    lon_obs = [float(item) for item in tmp_control[1::ncol]]
    ch_obs = [ int(item) for item in tmp_control[2::ncol] ]
    Yo_obs = [float(item) for item in tmp_control[3::ncol]]
    Yb_obs = [float(item) for item in tmp_control[4::ncol]]
    tmp_control = np.fromfile( output_file, sep=' ' ) # read output file
    Ya_obs = [float(item) for item in tmp_control[4::ncol]]

    lat_obs = np.array( lat_obs )
    lon_obs = np.array( lon_obs )
    ch_obs = np.array( ch_obs )
    Yo_obs = np.array( Yo_obs )
    Yb_obs = np.array( Yb_obs )
    Ya_obs = np.array( Ya_obs )

    d_MW =  {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'ch_obs':ch_obs, 'Yo_obs':Yo_obs, 'Yb_obs':Yb_obs, 'Ya_obs':Ya_obs}

    return d_MW


def plot_Tb(Storm, Exper_name, DAtime, sensor, mem):
   
    ch_num = [int(ich) for ich in dict_ss_ch[sensor]]
    # Define the low and high frequency for each sensor
    d_lowf = {'atms_npp':0, 'amsr2_gcom-w1':7, 'gmi_gpm':3, 'mhs_n19':0, 'mhs_n18':0, 'mhs_metop-a':0, 'mhs_metop-b':0, 'saphir_meghat':0, 'ssmis_f16': 13, 'ssmis_f17': 13, 'ssmis_f18': 13, 'ssmi_f15':1}

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = UD.read_wrf_domain( wrf_file )

    # Read data
    Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/MW/'+DAtime
    d_all = read_memTb(Hx_dir, sensor, DAtime, mem)

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(Storm, DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

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
    # 1st row: lf; 2nd row: hf
    for input_it in range(2):
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
        
        # Filter Tb at this channel
        if len(ch_num) == 2:
            ch_idx = d_all['ch_obs'] == ch_num[input_it]
        else:
            ch_idx = d_all['ch_obs'] == ch_num[0]

        Lat_obs_ch = d_all['lat_obs'][ch_idx] 
        Lon_obs_ch = d_all['lon_obs'][ch_idx]
        Yo_obs_ch = d_all['Yo_obs'][ch_idx]
        Yb_obspace = d_all['Yb_obs'][ch_idx]
        Ya_obspace = d_all['Ya_obs'][ch_idx] 

        # Find the obs over land
        if d_all['ch_obs'][ch_idx][0] == d_lowf[sensor]:
            is_ocean = globe.is_ocean(Lat_obs_ch, Lon_obs_ch)
            mask_x = is_ocean
        else:
            mask_x = np.full((np.size(Lon_obs_ch), ), True)

        # Obs
        ax[i,0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i,0].coastlines(resolution='10m', color='black',linewidth=0.5)
        ax[i,0].scatter(Lon_obs_ch[mask_x], Lat_obs_ch[mask_x],2.5,c=Yo_obs_ch[mask_x],\
                 edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())

        # Hxb
        ax[i,1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax[i,1].coastlines(resolution='10m', color='black',linewidth=0.5)
        ax[i,1].scatter(Lon_obs_ch[mask_x], Lat_obs_ch[mask_x],2.5,c=Yb_obspace[mask_x],\
                 edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree()) 
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
            ax[i,1].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())         
        # Hxa
        ax[i,2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax[i,2].coastlines(resolution='10m', color='black',linewidth=0.5)
        cs = ax[i,2].scatter(Lon_obs_ch[mask_x], Lat_obs_ch[mask_x],2.5,c=Ya_obspace[mask_x],\
                edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
            ax[i,2].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
    
    # Colorbar
    caxes = f.add_axes([0.2, 0.05, 0.6, 0.02])
    cbar = f.colorbar(cs, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)
    #plt.text( 0.8, 0.7, 'Brightness Temperature (K)', fontsize=6, transform=transAxes)
    
    #subplot title
    font = {'size':8,}
    ax[0,0].set_title('Yo', font, fontweight='bold')
    ax[0,1].set_title('H(Xb)', font, fontweight='bold')
    ax[0,2].set_title('H(Xa)', font, fontweight='bold')

    f.suptitle(Storm+': '+Exper_name+' Wrong AZ '+'{:03d}'.format(mem+1), fontsize=8, fontweight='bold')
    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2)) 
    for i in range(2):
        for j in range(3):
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
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
            gl.xlabel_style = {'size': 5}
            gl.ylabel_style = {'size': 6} 
   
    figure_des=small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/MW/Obspace/Ens/'+'{:03d}'.format(mem+1)+'_'+DAtime+'_'+sensor+'_Obspace_wrongAZ.png'
    plt.savefig(figure_des, dpi=400)
    print('Saving the figure: ', figure_des) 
  

def plot_Tb_diff(Storm, Exper_name, DAtime, sensor, dict_ss_len):

    ch_num = [int(ich) for ich in dict_ss_ch[sensor]]
    # Define the low and high frequency for each sensor
    d_lowf = {'atms_npp':0, 'amsr2_gcom-w1':7, 'gmi_gpm':3, 'mhs_n19':0, 'mhs_n18':0, 'mhs_metop-a':0, 'mhs_metop-b':0, 'saphir_meghat':0, 'ssmis_f16': 13, 'ssmis_f17': 13, 'ssmis_f18': 13, 'ssmi_f15':1}

    # Read WRF domain
    wrf_file = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = UD.read_wrf_domain( wrf_file )

    # Read data
    Hx_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Obs_Hx/MW/'+DAtime
    Tb_file = glob.glob( Hx_dir + '/mean_obs_res_d03*' + sensor + '*.txt' )
    d_all = read_Tb(Tb_file, sensor, dict_ss_len, d_wrf_d03)

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(Storm, DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # Prepare to calculate RMSE between the Hx and Yo
    metric = np.full(shape=(2,2), fill_value=None) # default type: none

    # ------------------ Plot -----------------------
    f, ax=plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(4,4.5), dpi=500)

    # Customize colormap
    max_T=20
    min_T=-20
    #min_RWB = 0
    #newRWB = Util_Vis.newRWB(max_T, min_T, min_RWB)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    # Loop over low frequency and high frequency if available
    # 1st row: lf; 2nd row: hf
    for input_it in range(2):
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
    
        # Filter Tb at this channel
        if len(ch_num) == 2:
            ch_idx = d_all['ch_obs'] == ch_num[input_it]
        else:
            ch_idx = d_all['ch_obs'] == ch_num[0]

        Lat_obs_ch = d_all['lat_obs'][ch_idx]
        Lon_obs_ch = d_all['lon_obs'][ch_idx]
        Yo_obs_ch = d_all['Yo_obs'][ch_idx]
        Yb_obspace = d_all['meanYb_obs'][ch_idx]
        Ya_obspace = d_all['meanYa_obs'][ch_idx]
        
        # Find the obs over land
        if d_all['ch_obs'][ch_idx][0] == d_lowf[sensor]:
            is_ocean = globe.is_ocean(Lat_obs_ch, Lon_obs_ch)
            mask_x = is_ocean
            #mask_x = np.full((np.size(Lon_obs_ch), ), True)
        else:
            mask_x = np.full((np.size(Lon_obs_ch), ), True)

        # HXb - Obs
        ax[i,0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i,0].coastlines(resolution='10m', color='black',linewidth=0.5)
        if plot_scatter:
            xb_s = ax[i,0].scatter(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],2.5,c=Yb_obspace[mask_x]-Yo_obs_ch[mask_x],\
                            edgecolors='none', cmap='bwr', vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
        else:
            masked_Tb = np.ma.masked_array(Yb_obspace, mask=mask_x)
            bounds = np.linspace(min_T,max_T,7)
            xb_s = ax[i,0].tricontourf(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],Yb_obspace[mask_x]-Yo_obs_ch[mask_x],\
                            cmap='bwr',vmin=min_T, vmax=max_T, levels=bounds, transform=ccrs.PlateCarree(), extend='both' )
        # add a reference point
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
            ax[i,0].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
        # calculate metric
        metric[i,0] = Bias(Yb_obspace[mask_x],Yo_obs_ch[mask_x])

        # HXa - Obs
        ax[i,1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax[i,1].coastlines(resolution='10m', color='black',linewidth=0.5)
        if plot_scatter:
            xa_s = ax[i,1].scatter(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],2.5,c=Ya_obspace[mask_x]-Yo_obs_ch[mask_x],\
                            edgecolors='none', cmap='bwr', vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
        else:
            bounds = np.linspace(min_T,max_T,7)
            xa_s = ax[i,1].tricontourf(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],Ya_obspace[mask_x]-Yo_obs_ch[mask_x],\
                            cmap='bwr',vmin=min_T, vmax=max_T, levels=bounds, transform=ccrs.PlateCarree(), extend='both' )
        # add a reference point
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
            ax[i,1].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
        # calculate metric
        metric[i,1] = Bias(Ya_obspace[mask_x],Yo_obs_ch[mask_x])


    # Colorbar
    cb_ticks = np.linspace(min_T, max_T, 9, endpoint=True)
    caxes = f.add_axes([0.2, 0.05, 0.6, 0.02])
    cbar = f.colorbar(xa_s, ticks=cb_ticks, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=8)
    #plt.text( 0.8, 0.7, 'Brightness Temperature (K)', fontsize=6, transform=transAxes)

    #subplot title
    font = {'size':8,}
    if metric[0,0] is None: 
        ax[0,0].set_title('LF:H(Xb)-Yo', font, fontweight='bold')
    else:
        metric_str = '%.2f' % metric[0,0]
        ax[0,0].set_title('LF:H(Xb)-Yo '+metric_str, font, fontweight='bold')
    
    if metric[0,1] is None:    
        ax[0,1].set_title('LF:H(Xa)-Yo', font, fontweight='bold')
    else:
        metric_str = '%.2f' % metric[0,1]
        ax[0,1].set_title('LF:H(Xa)-Yo '+metric_str, font, fontweight='bold')

    if metric[1,0] is None:    
        ax[1,0].set_title('HF:H(Xb)-Yo', font, fontweight='bold')
    else:
        metric_str = '%.2f' % metric[1,0]
        ax[1,0].set_title('HF:H(Xb)-Yo '+metric_str, font, fontweight='bold')

    if metric[1,1] is None:    
        ax[1,1].set_title('HF:H(Xa)-Yo', font, fontweight='bold')
    else:
        metric_str = '%.2f' % metric[1,1]
        ax[1,1].set_title('HF:H(Xa)-Yo '+metric_str, font, fontweight='bold')
    
    f.suptitle(Storm+': '+Exper_name+' Wrong AZ', fontsize=8, fontweight='bold')

    # Tick labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for i in range(2):
        for j in range(2):
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
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
            gl.xlabel_style = {'size': 5}
            gl.ylabel_style = {'size': 6}

    figure_des=small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/MW/Obspace/'+DAtime+'_'+sensor+'_Obspace_Diff_range20.png'
    plt.savefig(figure_des, dpi=400)
    print('Saving the figure: ', figure_des)


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'HARVEY'
    Exper_name = 'JerryRun/MW_THO/'
    Num_ens = 60

    start_time_str = '201708221200'
    end_time_str = '201708221200'
    Consecutive_times = False

    plot_scatter = True 
    # ------------------------------------------------------ 

    # Create MW DA times
    if not Consecutive_times:
        MW_times = ['201708221200',]
    else:
        MW_times = []
        exist_MW_times = os.listdir(big_dir+Storm+'/'+Exper_name+'/Obs_Hx/MW/20*')
        for it in exist_MW_times:
            if datetime.strptime(it,"%Y%m%d%H%M") < datetime.strptime(start_time_str,"%Y%m%d%H%M") \
                            and datetime.strptime(it,"%Y%m%d%H%M") > datetime.strptime(end_time_str,"%Y%m%d%H%M"):
                    continue
            else:
                MW_times.append( it )

    # Iterate thru each DAtime and plot Tb field
    for DAtime in MW_times:
        #obs_file_name = 'microwave_d03_' + DAtime + '_so'
        dict_ss_ch = ROMW.getSensor_Ch( small_dir+'Preprocess_Obs/toEnKFobs/MW/HARVEY/microwave_d03_201708221200_so_wrong')
        # Iterate thru each sensor and plot Tb comparison
        for sensor in dict_ss_ch:
            for imem in range(Num_ens):
                plot_Tb( Storm, Exper_name, DAtime, sensor,imem ) 
                #plot_Tb_diff( Storm, Exper_name, DAtime, sensor, dict_ss_len) 


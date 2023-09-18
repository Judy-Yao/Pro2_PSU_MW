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
import time

import Util_data as UD
import Read_Obspace_MW as ROMW

def read_memTb( Expers,wrf_dirs, DAtime, sensor, mem):

    # number of columns of each record
    ncol = 8
    d_MW_all = {}
    for iExper in Expers:
        idx_exper = Expers.index( iExper )
        mem_id = '{:03d}'.format(mem+1)
        input_file = wrf_dirs[idx_exper]+"/input_mem"+mem_id+'_d03_'+DAtime+'.tb.'+sensor+'.crtm.conv.txt'
        output_file = wrf_dirs[idx_exper]+"/output_mem"+mem_id+'_d03_'+DAtime+'.tb.'+sensor+'.crtm.conv.txt'
       
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

        d_Tb =  {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'ch_obs':ch_obs, 'Yo_obs':Yo_obs, 'Yb_obs':Yb_obs, 'Ya_obs':Ya_obs}
        d_MW_all[iExper] = d_Tb

    return d_MW_all

def read_meanTb( Expers,wrf_dirs, DAtime, sensor ):

    d_MW_all = {}
    # read out variables (not sorted)
    for iExper in Expers:
        idx_exper = Expers.index( iExper )
        Tb_file = wrf_dirs[idx_exper]+"/mean_obs_res_d03_" + DAtime + '.tb.' +  sensor + '.crtm.conv.txt'
        
        lat_obs = []
        lon_obs = []
        ch_obs = []
        Yo_obs = []
        meanYb_obs = []
        meanYa_obs = []

        # Read records
        print('Reading ', Tb_file)
        with open(Tb_file) as f:
            next(f)
            all_lines = f.readlines()

        for line in all_lines:
            split_line = line.split()
            lat_obs.append( float(split_line[0]) )
            lon_obs.append( float(split_line[1]) )
            ch_obs.append( int(split_line[2]) )
            Yo_obs.append( float(split_line[3]) )
            meanYb_obs.append( float(split_line[4]) )
            meanYa_obs.append( float(split_line[5]) )
        
        # Create a list tuples and sort
        mwso_tup = []
        for i in range(len(lat_obs)):
            mwso_tup.append( (lat_obs[i],lon_obs[i],ch_obs[i],Yo_obs[i],meanYb_obs[i],meanYa_obs[i]) )
        mwso_tup.sort(key=lambda a: (a[2], a[0], a[1], a[3]))
        mwso_list = [list( iso ) for iso in mwso_tup]

        lat = []
        lon = []
        ch = []
        Yo = []
        meanYb = []
        meanYa = [] 
        for i in range( np.shape(mwso_list)[0] ):
            lat.append( mwso_list[i][0] )
            lon.append( mwso_list[i][1] )
            ch.append( mwso_list[i][2] )
            Yo.append( mwso_list[i][3]  )
            meanYb.append( mwso_list[i][4] )
            meanYa.append( mwso_list[i][5] )

        lat = np.array( lat )
        lon = np.array( lon )
        ch = np.array( ch )
        Yo = np.array( Yo )
        meanYb = np.array( meanYb )
        meanYa = np.array( meanYa )
        d_Tb =  {'lat_obs':lat, 'lon_obs':lon, 'ch_obs':ch, 'Yo_obs':Yo, 'Yb_obs':meanYb, 'Ya_obs':meanYa}
        d_MW_all[iExper] = d_Tb
    
    # sanity check if different experiments have the same order of data
    assert d_MW_all[Expers[0]]['lat_obs'].all() == d_MW_all[Expers[1]]['lat_obs'].all()
    assert d_MW_all[Expers[0]]['Yo_obs'].all() == d_MW_all[Expers[1]]['Yo_obs'].all()
    return d_MW_all

def plot_Tb_Tbdiff( DAtime, sensor, mem):

    ch_num = [int(ich) for ich in dict_ss_ch[sensor]]
    # Define the low and high frequency for each sensor
    d_lowf = {'atms_npp':0, 'amsr2_gcom-w1':7, 'gmi_gpm':3, 'mhs_n19':0, 'mhs_n18':0, 'mhs_metop-a':0, 'mhs_metop-b':0, 'saphir_meghat':0, 'ssmis_f16': 13, 'ssmis_f17': 13, 'ssmis_f18': 13, 'ssmi_f15':1}

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = UD.read_wrf_domain( wrf_file )

    # Read Tbs of obs, Hxb, Hxa
    wrf_dirs = []
    for ie in Expers:
        wrf_dirs.append( big_dir+Storm+'/'+Exper_name+'/Obs_Hx/MW/'+ie )
    if mem == 60:
        d_all = read_meanTb( Expers,wrf_dirs, DAtime, sensor )
    else:
        d_all = read_memTb( Expers,wrf_dirs, DAtime, sensor, mem)

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(Storm, DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # ------------------ Plot -----------------------
    f, ax=plt.subplots(2, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(6,4), dpi=500)

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
            ch_idx = d_all[Expers[0]]['ch_obs'] == ch_num[input_it]
        else:
            ch_idx = d_all[Expers[0]]['ch_obs'] == ch_num[0]

        Lat_obs_ch = d_all[Expers[0]]['lat_obs'][ch_idx]
        Lon_obs_ch = d_all[Expers[0]]['lon_obs'][ch_idx]
        Yo_obs_ch = d_all[Expers[0]]['Yo_obs'][ch_idx]
        Yb_diff = d_all[Expers[0]]['Yb_obs'][ch_idx] - d_all[Expers[1]]['Yb_obs'][ch_idx]
        Ya_diff = d_all[Expers[0]]['Ya_obs'][ch_idx] - d_all[Expers[1]]['Ya_obs'][ch_idx]

        # Find the obs over land
        if d_all[Expers[0]]['ch_obs'][ch_idx][0] == d_lowf[sensor]:
            #is_ocean = globe.is_ocean(Lat_obs_ch, Lon_obs_ch)
            #mask_x = is_ocean
            mask_x = np.full((np.size(Lon_obs_ch), ), True)
        else:
            mask_x = np.full((np.size(Lon_obs_ch), ), True)
    
        #--- Obs
        max_T=300
        min_T=80
        min_Jet=150
        MWJet = Util_Vis.newJet(max_T=300, min_T=80, min_Jet=150)
        ax[i,0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i,0].coastlines(resolution='10m', color='black',linewidth=0.5)
        xo = ax[i,0].scatter(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],2.5,c=Yo_obs_ch[mask_x],\
            edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
        # add a reference point
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
            ax[i,0].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
        # calculate metric
        #metric[i,0] = Bias(Yb_obspace[mask_x],Yo_obs_ch[mask_x])
        # Colorbar
        caxes = f.add_axes([0.12, 0.05, 0.25, 0.02])
        xo_bar = f.colorbar(xo,ax=ax[1,0:1],orientation="horizontal", cax=caxes)
        xo_bar.ax.tick_params(labelsize=8)

        # Hxb: Exper 1 - Exper 2
        min_incre = -20
        max_incre = 20
        ax[i,1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax[i,1].coastlines(resolution='10m', color='black',linewidth=0.5)
        if plot_scatter:
            xb = ax[i,1].scatter(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],2.5,c=Yb_diff[mask_x],\
                edgecolors='none', cmap='bwr', vmin=min_incre, vmax=max_incre, transform=ccrs.PlateCarree())
        else:
            bounds = np.linspace(min_T,max_T,7)
            xb = ax[i,1].tricontourf(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],Yb_diff[mask_x],\
                            cmap='bwr',vmin=min_T, vmax=max_T, levels=bounds, transform=ccrs.PlateCarree(), extend='both' )
        # add a reference point
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
            ax[i,1].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
        # calculate metric
        #metric[i,1] = Bias(Ya_obspace[mask_x],Yo_obs_ch[mask_x])

        # Hxa: Exper 1 - Exper 2
        ax[i,2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax[i,2].coastlines(resolution='10m', color='black',linewidth=0.5)
        if plot_scatter:
            xa = ax[i,2].scatter(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],2.5,c=Ya_diff[mask_x],\
                edgecolors='none', cmap='bwr', vmin=min_incre, vmax=max_incre, transform=ccrs.PlateCarree())
        else:
            bounds = np.linspace(min_T,max_T,7)
            xa = ax[i,2].tricontourf(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],Ya_diff[mask_x],\
                            cmap='bwr',vmin=min_T, vmax=max_T, levels=bounds, transform=ccrs.PlateCarree(), extend='both' )
        # add a reference point
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
            ax[i,2].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
        # calculate metric
        #metric[i,1] = Bias(Ya_obspace[mask_x],Yo_obs_ch[mask_x])
        # Colorbar
        caxes = f.add_axes([0.4, 0.05, 0.5, 0.02])
        cb_diff_ticks = np.linspace(min_incre, max_incre, 9, endpoint=True)
        cbar = f.colorbar(xa, ax=ax[1,1:], ticks=cb_diff_ticks, orientation="horizontal", cax=caxes)
        cbar.ax.tick_params(labelsize=8)


    #subplot title
    font = {'size':8,}
    ax[0,0].set_title('Yo', font, fontweight='bold')
    ax[0,1].set_title('H(Xb) Diff', font, fontweight='bold')
    ax[0,2].set_title('H(Xa) Diff', font, fontweight='bold')

    title = Storm+': '+Exper_name+' '+'{:03d}'.format(mem+1)+'\ Serial - Random Obs'
    f.suptitle( title, fontsize=8, fontweight='bold')
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

    if mem == 60:
        figure_des=small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/MW/Obspace/mean_'+DAtime+'_'+sensor+'_Obspace_expers.png'
    else: 
        figure_des=small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/MW/Obspace/Ens/'+'{:03d}'.format(mem+1)+'_'+DAtime+'_'+sensor+'_Obspace_expers.png'
    plt.savefig(figure_des, dpi=400)
    print('Saving the figure: ', figure_des)




if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'HARVEY'
    Exper_name = 'JerryRun/MW_THO/'
    Expers = ['201708221200_correct_azimuth_DA_CRTM','201708221200_correct_azimuth_DA_CRTM_randomObs',]
    Num_ens = 60

    start_time_str = '201708221200'
    end_time_str = '201708221200'
    Consecutive_times = False

    plot_scatter = True
    plot_mean_diff = True
    plot_mem_diff = False
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

    for DAtime in MW_times:
        #obs_file_name = 'microwave_d03_' + DAtime + '_so'
        dict_ss_ch = ROMW.getSensor_Ch( small_dir+'Preprocess_Obs/toEnKFobs/MW/HARVEY/microwave_d03_201708221200_so_correct')

        if plot_mean_diff: # Plot difference of ensemble mean between two experiments
            dict_ss_len = {} # sensor: len_records_before (used to check if the length of readed records matches the length of pre-processed records)
            # Iterate thru each sensor and calculate mean of Yb and Ya         
            for sensor in dict_ss_ch:
                print('------------ Calculate mean of Hx --------------')
                #Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/MW/'+DAtime
                for ie in Expers:
                    Hx_dir = big_dir+Storm+'/'+Exper_name+'/'+'/Obs_Hx/MW/'+ie
                    if sensor == 'gmi_gpm':
                        len_pre_all = 0
                        len_pre_processed = ROMW.write_mean_eachSensor( Hx_dir, 'gmi_gpm_lf' )
                        len_pre_all = len_pre_all + len_pre_processed

                        len_pre_processed = ROMW.write_mean_eachSensor( Hx_dir, 'gmi_gpm_hf' )
                        len_pre_all = len_pre_all + len_pre_processed
                        dict_ss_len[sensor] = len_pre_all
                    else:
                        len_pre_processed = ROMW.write_mean_eachSensor( Hx_dir, sensor )
                        #time.sleep(10) # wait enough long time to get the file written
                        dict_ss_len[sensor] = len_pre_processed

            # Iterate thru each sensor and plot Tb comparison
            for sensor in dict_ss_ch:
                plot_Tb_Tbdiff( DAtime, sensor, 60)

        elif plot_mem_diff:

            for sensor in dict_ss_ch:
                for imem in range(Num_ens):
                    plot_Tb_Tbdiff( DAtime, sensor, imem) 






























 

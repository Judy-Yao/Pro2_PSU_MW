#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os,fnmatch
import glob
import numpy as np
import Util_Vis
import netCDF4 as nc
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from global_land_mask import globe
import math
import time
from datetime import datetime

import Util_data as UD
import Util_Vis
import Diagnostics as Diag

def RMSE(simu, obs):
    return np.sqrt( ((simu - obs) ** 2).mean() )

def Bias(simu, obs):
    return  np.sum((simu - obs),0)/np.size(obs,0)


def SensorCh_to_Freq( sensor, ch):

    if 'amsr2_' in sensor:
        if ch == '7' or ch == 7:
            return '18.7GHzV' 
        elif ch == '13' or ch == 13:
             return '89GHzV'
    elif 'atms_' in sensor:
        if ch == '18' or ch == 18:
            return '183.31+-7GHzH'
    elif 'gmi_' in sensor:
        if ch == '3' or ch == 3:
            return '18.7GHzV'
        elif ch == '13' or ch == 13:
            return '183.31+/-7GHzV'
    elif 'mhs_' in sensor:
        if ch == '5' or ch == 5:
            return '190.31GHzV' 
    elif 'saphir_' in sensor:
        if ch == '5' or ch == 5:
            return '183.31+-6.8GHz'
    elif 'ssmi_' in sensor:
        if ch == '1' or ch == 1:
            return 'fcdr_tb19v'
        elif ch == '6' or ch == 6:
            return 'fcdr_tb85v'
    elif 'ssmis_' in sensor:
        if ch == '13' or ch == 13:
            return '19.35GHzV'
        elif ch == '9' or ch == 9:
            return '183.31+-6.6GHzH'

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

    # Special treatment to gmi_gpm sensor
    gmi_ch = [key_sensor for key_sensor in dict_ss_ch if 'gmi_gpm' in key_sensor]
    if len( gmi_ch ) >= 1: 
        rem_key = ['gmi_gpm_lf','gmi_gpm_hf']
        add_ss = 'gmi_gpm'
        add_ch_num = [dict_ss_ch[ikey][0] for ikey in rem_key]
        [dict_ss_ch.pop(ikey) for ikey in rem_key]
        dict_ss_ch[add_ss] = add_ch_num
        return dict_ss_ch
    ##idx_gmi = list(dict_ss_ch.key()).index(key_sensor) for key_sensor in dict_ss_ch if 'gmi_gpm' in key_sensor
    else:
        return dict_ss_ch


def write_mean_eachSensor( Hx_dir, sensor ):

    # Number of ensemble members
    num_ens = 60
    
    # number of columns of each record
    ncol = 8

    # List the Yb and Ya files
    file_yb = sorted( glob.glob(Hx_dir + '/input_mem0*' + sensor + '*.txt') ) 
    file_ya = sorted( glob.glob(Hx_dir + '/output_mem0*' + sensor + '*.txt') )

    # Sanity check the number of calculated ens
    if np.size(file_yb) != num_ens:
        raise ValueError('Wait until all of the ens is calculated!')
    if np.size(file_ya) != num_ens:
        raise ValueError('Wait until all of the ens is calculated!')

    # Read obs from a member as the control variable
    tmp_control = np.fromfile( file_yb[0], sep=' ' ) # read member 1
    Lat_all = [ "{0:.3f}".format(item) for item in tmp_control[0::ncol] ] #latitude
    Lon_all = [ "{0:.3f}".format(item) for item in tmp_control[1::ncol] ] #longitude
    Chnum_all = [ int(item) for item in tmp_control[2::ncol] ]
    Control_yo = [ "{0:.3f}".format(item) for item in tmp_control[3::ncol] ] #observed Tb
    FOV_aScan = [ "{0:.3f}".format(item) for item in tmp_control[5::ncol] ] #FOV along scan
    FOV_cScan = [ "{0:.3f}".format(item) for item in tmp_control[6::ncol] ] #FOV cross scan 
    Azimuth_angle = [ "{0:.3f}".format(item) for item in tmp_control[7::ncol] ] #Azimuth angle

    # ---- Read prior Tb from the ens ----
    sum_yb = np.zeros( shape=np.shape(Control_yo) )
    # Iterate thru input ens
    for ifile in file_yb:
        print('Reading the file: ' + ifile)
        tmp = np.fromfile( ifile, sep=' ' )
        
        # Sanity check
        yo = [ "{0:.3f}".format(item) for item in tmp[3::ncol] ]  #tmp[3::ncol]
        if not ( yo == Control_yo ): 
            raise ValueError('Records of this member do not match the control member!')
        sum_yb = sum_yb + tmp[4::ncol]
    
    Yb_all_mean = np.round( (sum_yb / num_ens), 3 )

    # ---- Read posterior Tb from the ens -----
    sum_ya = np.zeros( shape=np.shape(Control_yo) )
    # Iterate thru output ens
    for ifile in file_ya:
        #print('Reading the file: ' + ifile)
        tmp = np.fromfile( ifile, sep=' ' )

        # Sanity check
        yo = [ "{0:.3f}".format(item) for item in tmp[3::ncol] ] 
        
        if not ( yo == Control_yo ): 
            raise ValueError('Records of this member do not match the control member!')  
        sum_ya = sum_ya + tmp[4::ncol] 
    
    Ya_all_mean = np.round( (sum_ya / num_ens), 3 )

    # Stack each list into an array
    all_attrs = np.column_stack( (Lat_all, Lon_all, Chnum_all, Control_yo, Yb_all_mean, Ya_all_mean, FOV_aScan, FOV_cScan, Azimuth_angle) )

    # ---- Write to file and save it to the disk ----
    header = ['Lat','Lon','Ch_num','Tb_obs','Tb_conv_Yb','Tb_conv_Ya','efov_aScan','efov_cScan','Azimuth_angle']
    file_name = ifile.replace("output_mem060_d03", "mean_obs_res_d03" ) 
    with open(file_name,'w') as f:
        # Add header 
        f.write('\t'.join( item.rjust(8) for item in header ) + '\n' ) 
        # Write the record to the file serially
        len_records = np.shape( all_attrs )[0]
        for irow in range( len_records ):
            irecord =  [str(item) for item in all_attrs[irow,:] ]
            f.write('\t'.join( item.rjust(8) for item in irecord ) + '\n') 

    return( len_records )  

def read_Tb(Tb_file, sensor, dict_ss_len, d_wrf_d03):

    lat_obs = []
    lon_obs = []
    ch_obs = []
    Yo_obs = []
    meanYb_obs = []
    meanYa_obs = []

    # Define the wrf domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    # number of columns of each record
    #ncol = 9

    # Read records
    for ifile in Tb_file:
        print('Reading ', ifile)
        with open(ifile) as f:
            next(f)
            all_lines = f.readlines() 
        
        for line in all_lines:
            split_line = line.split()
            
            read_lat = float(split_line[0])
            read_lon = float(split_line[1])
            if read_lat <= lat_max and read_lat >= lat_min and read_lon <= lon_max and read_lon >= lon_min:
                lat_obs.append( read_lat )
                lon_obs.append( read_lon )
                ch_obs.append( int(split_line[2]) )
                Yo_obs.append( float(split_line[3]) ) 
                meanYb_obs.append( float(split_line[4]) )
                meanYa_obs.append( float(split_line[5]) )
            else:
                continue

    #if np.size(lat_obs) != dict_ss_len[sensor]:
    #    raise ValueError('The length of post-processed file is not equal to the pre-processed file!')

    lat_obs = np.array( lat_obs )
    lon_obs = np.array( lon_obs )
    ch_obs = np.array( ch_obs )
    Yo_obs = np.array( Yo_obs )
    meanYb_obs = np.array( meanYb_obs )
    meanYa_obs = np.array( meanYa_obs )

    dict_Tb_all = {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'ch_obs':ch_obs, 'Yo_obs':Yo_obs, 'meanYb_obs':meanYb_obs, 'meanYa_obs':meanYa_obs}
    return dict_Tb_all


def plot_Tb(Storm, Exper_name, DAtime, sensor, dict_ss_len):
   
    ch_num = [int(ich) for ich in dict_ss_ch[sensor]]
    # Define the low and high frequency for each sensor
    d_lowf = {'atms_npp':0, 'amsr2_gcom-w1':7, 'gmi_gpm':3, 'mhs_n19':0, 'mhs_n18':0, 'mhs_metop-a':0, 'mhs_metop-b':0, 'saphir_meghat':0, 'ssmis_f16': 13, 'ssmis_f17': 13, 'ssmis_f18': 13, 'ssmi_f15':1}

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = UD.read_wrf_domain( wrf_file )

    # Read data
    Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/MW/'+DAtime
    Tb_file = glob.glob( Hx_dir + '/mean_obs_res_d03*' + sensor + '*.txt' )
    d_all = read_Tb(Tb_file, sensor, dict_ss_len, d_wrf_d03)

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
        Yb_obspace = d_all['meanYb_obs'][ch_idx]
        Ya_obspace = d_all['meanYa_obs'][ch_idx] 

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
  
        # Annotation
        if len(ch_num) == 2:
            if ch_num[input_it] == d_lowf[sensor]:
                f.text( 0.03,0.6,SensorCh_to_Freq(sensor,ch_num[input_it]),rotation='vertical')
            else:
                f.text( 0.03,0.15,SensorCh_to_Freq(sensor,ch_num[input_it]),rotation='vertical')
        else:
            if ch_num[0] == d_lowf[sensor]:
                f.text( 0.03,0.6,s=SensorCh_to_Freq(sensor,ch_num[0]),rotation='vertical',)
            else:
                f.text( 0.03,0.15,SensorCh_to_Freq(sensor,ch_num[0]),rotation='vertical',)

    # Colorbar
    caxes = f.add_axes([0.2, 0.05, 0.6, 0.02])
    cbar = f.colorbar(cs, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=8)
    #plt.text( 0.8, 0.7, 'Brightness Temperature (K)', fontsize=6,)

    #subplot title
    matplotlib.rcParams['mathtext.fontset'] = 'custom'
    matplotlib.rcParams['mathtext.bf'] = 'STIXGeneral:italic:bold'
    font = {'size':9,}
    ax[0,0].set_title('Yo', font, fontweight='bold')
    ax[0,1].set_title(r'$\mathbf{\overline{H(Xb)}}$', font, )
    ax[0,2].set_title(r'$\mathbf{\overline{H(Xa)}}$', font) 

    title = Storm+': '+Exper_name+'\n'+DAtime+'          ('+sensor+')'
    f.suptitle(title, fontsize=7, fontweight='bold')
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
            gl.xlabel_style = {'size': 6}
            gl.ylabel_style = {'size': 6} 
   
    figure_des=small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/MW_Obspace/'+DAtime+'_'+sensor+'_Obspace.png'
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
        metric[i,0] = RMSE(Yb_obspace[mask_x],Yo_obs_ch[mask_x])

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
        metric[i,1] = RMSE(Ya_obspace[mask_x],Yo_obs_ch[mask_x])


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
    
    f.suptitle(Storm+': '+Exper_name, fontsize=8, fontweight='bold')

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
    Storm = 'JOSE'
    DA = 'IR+MW'
    MP = 'WSM6'

    start_time_str = '201709050000'
    end_time_str = '201709070000'
    Consecutive_times = True

    plot_full = True
    plot_diff = False
    plot_scatter = True 
    # ------------------------------------------------------ 

    # Create experiment names
    Exper_name =  UD.generate_one_name( Storm,DA,MP )
    Exper_obs =  UD.generate_one_name( Storm,'IR+MW',MP )
    # Create MW DA times
    if not Consecutive_times:
        MW_times = ['201709062200',]
    else:
        MW_times = []
        exist_MW_times = sorted(fnmatch.filter(os.listdir( big_dir+Storm+'/'+Exper_name+'/Obs_Hx/MW/' ),'20*'))
        for it in exist_MW_times:
            if datetime.strptime(it,"%Y%m%d%H%M") < datetime.strptime(start_time_str,"%Y%m%d%H%M") \
                            or datetime.strptime(it,"%Y%m%d%H%M") > datetime.strptime(end_time_str,"%Y%m%d%H%M"):
                    continue
            else:
                MW_times.append( it )

    # Create plot dirs
    if plot_full:
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/MW_Obspace/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

    if plot_diff:
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/MW_Obspace_Diff/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

    # Iterate thru each DAtime and plot Tb field
    print( MW_times )
    for DAtime in MW_times:
        obs_file_name = 'microwave_' + DAtime + '_so'
        dict_ss_ch = getSensor_Ch( small_dir+Storm+'/Obs_y/MW/Processed_2nd_time/'+Exper_obs+'/'+obs_file_name )

        dict_ss_len = {} # sensor: len_records_before (used to check if the length of readed records matches the length of pre-processed records)
        # Iterate thru each sensor and calculate mean of Yb and Ya         
        for sensor in dict_ss_ch:
            print('------------ Calculate mean of Hx --------------')
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/MW/'+DAtime
            print(Hx_dir)
            if sensor == 'gmi_gpm':
                len_pre_all = 0
                len_pre_processed = write_mean_eachSensor( Hx_dir, 'gmi_gpm_lf' )
                len_pre_all = len_pre_all + len_pre_processed

                len_pre_processed = write_mean_eachSensor( Hx_dir, 'gmi_gpm_hf' )
                len_pre_all = len_pre_all + len_pre_processed
                dict_ss_len[sensor] = len_pre_all
            else:
                len_pre_processed = write_mean_eachSensor( Hx_dir, sensor )
                #time.sleep(10) # wait enough long time to get the file written
                dict_ss_len[sensor] = len_pre_processed 

        # Iterate thru each sensor and plot Tb comparison
        for sensor in dict_ss_ch:
            print(sensor)
            print(dict_ss_ch[sensor])
            print('------------ Plot ----------------------')
            if plot_full:
                plot_Tb( Storm, Exper_name, DAtime, sensor, dict_ss_len) 
            if plot_diff:
                plot_Tb_diff( Storm, Exper_name, DAtime, sensor, dict_ss_len) 


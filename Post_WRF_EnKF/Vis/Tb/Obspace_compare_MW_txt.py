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



def RMSE(simu, obs):
    return np.sqrt( ((simu - obs) ** 2).mean() )


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

def read_bestrack(btk_file, DAtime):

    with open(btk_file) as f:
        all_lines = f.readlines()

    # Process all of records to our format/unit 
    btk_lat = []
    btk_lon = []
    for line in all_lines:
        # split one record into different parts
        split_line = line.split()
        # Read time
        btk_time = split_line[2].replace(',','') + '00'
        if btk_time == DAtime:
            # Read latitude
            lat_line = split_line[6].replace(',','')
            if 'N' in lat_line:
                btk_lat.append(float(lat_line.replace('N',''))/10)
            else:
                btk_lat.append(0-float(lat_line.replace('S',''))/10)
            # Read longitute
            lon_line = split_line[7].replace(',','')
            if 'W' in lon_line:
                btk_lon.append(0-float(lon_line.replace('W',''))/10)
            else:
                btk_lon.append(float(lon_line.replace('E',''))/10)

            return btk_lon, btk_lat


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
        print('Reading the file: ' + ifile)
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
    file_name = ifile.replace("output_mem060", "mean" ) 
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
    wrf_file = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = read_wrf_domain( wrf_file )

    # Read data
    Hx_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Obs_Hx/MW/'+DAtime
    Tb_file = glob.glob( Hx_dir + '/mean_d03*' + sensor + '*.txt' )
    d_all = read_Tb(Tb_file, sensor, dict_ss_len, d_wrf_d03)

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = read_TCvitals('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # Read location from best-track
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        Best_track_file = os.listdir('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/' + Storm + '/Post_Storm_btk')
        btk_lon, btk_lat = read_bestrack('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/' + Storm + '/Post_Storm_btk/' + Best_track_file[0], DAtime)
        print( 'Location from Best-track: ', btk_lon, btk_lat )


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

        if len(ch_num) == 2:
            ch_idx = d_all['ch_obs'] == ch_num[input_it]
        else:
            ch_idx = d_all['ch_obs'] == ch_num[0]

        Lat_obs_ch = d_all['lat_obs'][ch_idx] 
        Lon_obs_ch = d_all['lon_obs'][ch_idx]
        Yo_obs_ch = d_all['Yo_obs'][ch_idx]
        Yb_obspace = d_all['meanYb_obs'][ch_idx]
        Ya_obspace = d_all['meanYa_obs'][ch_idx] 

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
    caxes = f.add_axes([0.2, 0.97, 0.6, 0.02])
    cbar = f.colorbar(cs, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)
    #plt.text( 0.8, 0.7, 'Brightness Temperature (K)', fontsize=6, transform=transAxes)
    
    #subplot title
    font = {'size':8,}
    ax[0,0].set_title('Yo', font, fontweight='bold')
    ax[0,1].set_title('H(Xb)', font, fontweight='bold')
    ax[0,2].set_title('H(Xa)', font, fontweight='bold')

    # Axis labels
    #lon_gridlines = list(range(math.floor(lon_min)-1, math.ceil(lon_max)+1,1))
    #lat_gridlines = list(range(math.floor(lat_min)-1, math.ceil(lat_max)+1,1))
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
    
    plt.savefig('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Vis_analyze/Tb/MW/Obspace/'+DAtime+'_'+sensor+'_Obspace.png', dpi=300)
  

def plot_Tb_diff(Storm, Exper_name, DAtime, sensor, dict_ss_len):

    ch_num = [int(ich) for ich in dict_ss_ch[sensor]]
    # Define the low and high frequency for each sensor
    d_lowf = {'atms_npp':0, 'amsr2_gcom-w1':7, 'gmi_gpm':3, 'mhs_n19':0, 'mhs_n18':0, 'mhs_metop-a':0, 'mhs_metop-b':0, 'saphir_meghat':0, 'ssmis_f16': 13, 'ssmis_f17': 13, 'ssmis_f18': 13, 'ssmi_f15':1}

    # Read WRF domain
    wrf_file = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = read_wrf_domain( wrf_file )

    # Read data
    Hx_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Obs_Hx/MW/'+DAtime
    Tb_file = glob.glob( Hx_dir + '/mean_d03*' + sensor + '*.txt' )
    d_all = read_Tb(Tb_file, sensor, dict_ss_len, d_wrf_d03)

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = read_TCvitals('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # Read location from best-track
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        Best_track_file = os.listdir('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/' + Storm + '/Post_Storm_btk')
        btk_lon, btk_lat = read_bestrack('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/' + Storm + '/Post_Storm_btk/' + Best_track_file[0], DAtime)
        print( 'Location from Best-track: ', btk_lon, btk_lat )

    # Prepare to calculate RMSE between the Hx and Yo
    rmse = np.full(shape=(2,2), fill_value=None) # default type: none

    # ------------------ Plot -----------------------
    f, ax=plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(4,4.5), dpi=500)

    # Customize colormap
    max_T=50
    min_T=-50
    #min_RWB = 0
    #newRWB = Util_Vis.newRWB(max_T, min_T, min_RWB)

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

        if len(ch_num) == 2:
            ch_idx = d_all['ch_obs'] == ch_num[input_it]
        else:
            ch_idx = d_all['ch_obs'] == ch_num[0]

        Lat_obs_ch = d_all['lat_obs'][ch_idx]
        Lon_obs_ch = d_all['lon_obs'][ch_idx]
        Yo_obs_ch = d_all['Yo_obs'][ch_idx]
        Yb_obspace = d_all['meanYb_obs'][ch_idx]
        Ya_obspace = d_all['meanYa_obs'][ch_idx]

        #if ch_num[input_it] == d_lowf[sensor_short]:
        if d_all['ch_obs'][ch_idx][0] == d_lowf[sensor]:
            is_ocean = globe.is_ocean(Lat_obs_ch, Lon_obs_ch)
            mask_x = is_ocean
        else:
            mask_x = np.full((np.size(Lon_obs_ch), ), True)


        # HXb - Obs
        ax[i,0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i,0].coastlines(resolution='10m', color='black',linewidth=0.5)
        ax[i,0].scatter(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],2.5,c=Yb_obspace[mask_x]-Yo_obs_ch[mask_x],edgecolors='none', cmap='RdBu_r', vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
            ax[i,0].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
        
        rmse[i,0] = RMSE(Yb_obspace, Yo_obs_ch)

        # HXa - Obs
        ax[i,1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax[i,1].coastlines(resolution='10m', color='black',linewidth=0.5)
        cs = ax[i,1].scatter(Lon_obs_ch[mask_x], Lat_obs_ch[mask_x],2.5,c=Ya_obspace[mask_x]-Yo_obs_ch[mask_x],\
                 edgecolors='none', cmap='RdBu_r', vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
        if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
            ax[i,1].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())

        rmse[i,1] = RMSE(Ya_obspace, Yo_obs_ch)


    # Colorbar
    cb_ticks = np.linspace(min_T, max_T, 5, endpoint=True)
    caxes = f.add_axes([0.2, 0.97, 0.6, 0.02])
    cbar = f.colorbar(cs, ticks=cb_ticks, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)
    #plt.text( 0.8, 0.7, 'Brightness Temperature (K)', fontsize=6, transform=transAxes)

    #subplot title
    font = {'size':8,}
    if rmse[0,0] is None: 
        ax[0,0].set_title('LF:H(Xb)-Yo', font, fontweight='bold')
    else:
        rmse_str = '%.2f' % rmse[0,0]
        ax[0,0].set_title('LF:H(Xb)-Yo '+rmse_str, font, fontweight='bold')
    
    if rmse[0,1] is None:    
        ax[0,1].set_title('LF:H(Xa)-Yo', font, fontweight='bold')
    else:
        rmse_str = '%.2f' % rmse[0,1]
        ax[0,1].set_title('LF:H(Xa)-Yo '+rmse_str, font, fontweight='bold')

    if rmse[1,0] is None:    
        ax[1,0].set_title('HF:H(Xb)-Yo', font, fontweight='bold')
    else:
        rmse_str = '%.2f' % rmse[1,0]
        ax[1,0].set_title('HF:H(Xb)-Yo '+rmse_str, font, fontweight='bold')

    if rmse[1,1] is None:    
        ax[1,1].set_title('HF:H(Xa)-Yo', font, fontweight='bold')
    else:
        rmse_str = '%.2f' % rmse[1,1]
        ax[1,1].set_title('HF:H(Xa)-Yo '+rmse_str, font, fontweight='bold')
    
    f.suptitle(Storm+': '+Exper_name, fontsize=8, fontweight='bold')

    # Tick labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for i in range(2):
        for j in range(2):
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


    plt.savefig('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Vis_analyze/Tb/MW/Obspace/'+DAtime+'_'+sensor+'_Diff_Obspace_range50.png', dpi=300)
    print('Saving the figure: ', '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Vis_analyze/Tb/MW/Obspace/'+DAtime+'_'+sensor+'_Diff_Obspace.png')



if __name__ == '__main__':

    Storm = 'HARVEY'
    Exper_name =  'J_DA+J_WRF+J_init'#'JerryRun/MW_THO/'#'J_DA+J_WRF+J_init'
    Exper_obs = 'J_DA+J_WRF+J_init'

    #MW_times = ['201708230800','201708231200']
    MW_times = ['201708221200',]#'201708221300','201708221900','201708222000','201708222100','201708222300','201708230000','201708230800','201708231200']
    #MW_times = ['201709161800','201709161900','201709162100','201709162200','201709162300','201709170100','201709170400','201709170500','201709170700','201709170900','201709171000','201709171100','201709171300','201709171700']

#    Plot_scatter = True # Plot Tbs on scattered grid points OR plt.contourf

    # Iterate thru each DAtime
    for DAtime in MW_times:
        #obs_file_name = 'microwave_d03_' + DAtime + '_so'
        obs_file_name = 'microwave_' + DAtime + '_so'
        #dict_ss_ch = getSensor_Ch( '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/Obs_y/MW/Processed_1st_time/'+obs_file_name )
        dict_ss_ch = getSensor_Ch( '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/Obs_y/MW/Processed_2nd_time/'+Exper_obs+'/'+obs_file_name )

        dict_ss_len = {} # sensor: len_records_before (used to check if the length of readed records matches the length of pre-processed records)
        # Iterate thru each sensor and calculate mean of Yb and Ya         
        for sensor in dict_ss_ch:
            print('------------ Calculate mean of Hx --------------')
            Hx_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Obs_Hx/MW/'+DAtime
            if sensor == 'gmi_gpm':
                len_pre_all = 0
                len_pre_processed = write_mean_eachSensor( Hx_dir, 'gmi_gpm_lf' )
                time.sleep(10) # wait enough long time to get the file written
                len_pre_all = len_pre_all + len_pre_processed

                len_pre_processed = write_mean_eachSensor( Hx_dir, 'gmi_gpm_hf' )
                time.sleep(10)
                len_pre_all = len_pre_all + len_pre_processed
                dict_ss_len[sensor] = len_pre_all
            else:
                len_pre_processed = write_mean_eachSensor( Hx_dir, sensor )
                time.sleep(10) # wait enough long time to get the file written
                dict_ss_len[sensor] = len_pre_processed 

        # Iterate thru each sensor and plot Tb comparison
        for sensor in dict_ss_ch:
            print(sensor)
            print(dict_ss_ch[sensor])
            print('------------ Plot ----------------------')
            plot_Tb( Storm, Exper_name, DAtime, sensor, dict_ss_len) 
            plot_Tb_diff( Storm, Exper_name, DAtime, sensor, dict_ss_len) 


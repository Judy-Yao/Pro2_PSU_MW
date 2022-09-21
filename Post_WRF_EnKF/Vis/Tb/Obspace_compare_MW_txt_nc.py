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
import matlab.engine


def read_obs(obs_file, sensor):

    with open(obs_file) as tmp:
        obs_all = tmp.readlines()

    Ch_obs = []
    Lat_obs = []
    Lon_obs = []
    Yo_obs = []
    for line in obs_all:
        line_split = line.split()
        if sensor in line_split[1]:
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


def read_simu_Tb(Hxb_files, Hxa_files):
    
    idx_file = 0
    if len(Hxb_files) >1:
        ncdir = nc.Dataset(Hxb_files[0],'r')
        dimension = ncdir.variables['XLAT'].shape
        dimension_y = dimension[1]
        dimension_x = dimension[2]

        Lat_x = np.zeros([len(Hxb_files), dimension_y, dimension_x])
        Lon_x = np.zeros([len(Hxb_files), dimension_y, dimension_x])
        Ch_x = np.zeros(len(Hxb_files))
        Yb_x = np.zeros([len(Hxb_files), dimension_y, dimension_x])
        Ya_x = np.zeros([len(Hxb_files), dimension_y, dimension_x])

        for Hxb_file in Hxb_files:
            print('Hxb_file:' + Hxb_file)
            ncdir = nc.Dataset(Hxb_file, 'r')

            Lat_x[idx_file,:,:] =  ncdir.variables['XLAT'][0,:,:]
            Lon_x[idx_file,:,:] =  ncdir.variables['XLONG'][0,:,:]
            Ch_x[idx_file] = ncdir.variables['ch'][:] 
            Yb_x[idx_file,:,:] =  ncdir.variables['Brightness_Temperature'][0,:,:,:]
            
            ncdir = nc.Dataset(Hxa_files[idx_file])
            Ya_x[idx_file,:,:] =  ncdir.variables['Brightness_Temperature'][0,:,:,:]

            idx_file = idx_file + 1
    else:
        Hxb_file = Hxb_files[0]
        print('Hxb_file:' + Hxb_file)
        ncdir = nc.Dataset(Hxb_file, 'r')

        Lat_x = ncdir.variables['XLAT'][0,:,:] #latitude: XLAT(time, y, x)
        Lon_x = ncdir.variables['XLONG'][0,:,:] #longitude: XLONG(time, y, x)
        Ch_x = ncdir.variables['ch'][:] # int ch(ch)
    
        Yb_x  = ncdir.variables['Brightness_Temperature'][0,:,:,:] #BackgroundTB:float Brightness_Temperature(time, ch, y, x)    

        Hxa_file = Hxa_files[0]
        ncdir = nc.Dataset(Hxa_file, 'r') 
        Ya_x = ncdir.variables['Brightness_Temperature'][0,:,:,:]
        # Comment: in matlab, data-read order is reversed on the data-storage order.
   
    dict_simu_Tb = {'Ch_x': Ch_x, 'Lat_x': Lat_x, 'Lon_x': Lon_x, 'Yb_x': Yb_x, 'Ya_x': Ya_x}
    return dict_simu_Tb


    
def plot_Tb(Storm, Exper_name, DAtime, sensor):
   
    ch_num = [int(ich) for ich in dict_ss_ch[sensor]]
    # Define the low and high frequency for each sensor
    d_lowf = {'atms_npp':0, 'amsr2_gcom-w1':7, 'gmi_gpm':3, 'mhs_n19':0, 'mhs_n18':0, 'mhs_metop-a':0, 'mhs_metop-b':0, 'saphir_meghat':0, 'ssmis_f16': 13, 'ssmis_f17': 13, 'ssmis_f18': 13, 'ssmi_f15':1}

    # Read obs data
    obs_file_name = 'microwave_' + DAtime + '_so'
    d_obs = read_obs('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/Obs_y/MW/'+obs_file_name, sensor)
    # Read simulated Tbs 
    Hx_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Obs_Hx/MW/'+DAtime
    if sensor == 'gmi_gpm':
        Hxb_file = [Hx_dir+'/wrf_enkf_input_d03_mean_'+DAtime+'_tb.gmi_gpm_lf.crtm.nc', Hx_dir+'/wrf_enkf_input_d03_mean_'+DAtime+'_tb.gmi_gpm_hf.crtm.nc']
        Hxa_file = [Hx_dir+'/wrf_enkf_output_d03_mean_'+DAtime+'_tb.gmi_gpm_lf.crtm.nc', Hx_dir+'/wrf_enkf_output_d03_mean_'+DAtime+'_tb.gmi_gpm_hf.crtm.nc']
    else:
        Hxb_file = [Hx_dir+'/wrf_enkf_input_d03_mean_'+DAtime+'_tb.'+sensor+'.crtm.nc',]
        Hxa_file = [Hx_dir+'/wrf_enkf_output_d03_mean_'+DAtime+'_tb.'+sensor+'.crtm.nc',]  
    
    d_simu = read_simu_Tb( Hxb_file, Hxa_file )
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
    lat_min = np.amin(d_simu['Lat_x'].flatten())
    lat_max = np.amax(d_simu['Lat_x'].flatten())
    lon_min = np.amin(d_simu['Lon_x'].flatten())
    lon_max = np.amax(d_simu['Lon_x'].flatten())


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

        # ---------- obs ---------------------
        if len(ch_num) == 2:
            ch_idx = d_obs['Ch_obs'] == ch_num[input_it]
        else:
            ch_idx = d_obs['Ch_obs'] == ch_num[0]

        Lat_obs_ch = d_obs['Lat_obs'][ch_idx] 
        Lon_obs_ch = d_obs['Lon_obs'][ch_idx]
        Yo_obs_ch = d_obs['Yo_obs'][ch_idx]
   
        #if ch_num[input_it] == d_lowf[sensor_short]:
        #if d_obs['Ch_obs'][ch_idx][0] == d_lowf[sensor]:
        #    is_ocean = globe.is_ocean(Lat_obs_ch, Lon_obs_ch)
        #    mask_x = is_ocean
        #else:
        #    mask_x = np.full((np.size(Lon_obs_ch), ), True)
   
        ax[i,0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax[i,0].coastlines(resolution='10m', color='black',linewidth=0.5)
        #ax[i,0].scatter(Lon_obs_ch[mask_x],Lat_obs_ch[mask_x],2.5,c=Yo_obs_ch[mask_x],edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
        ax[i,0].scatter(Lon_obs_ch,Lat_obs_ch,2.5,c=Yo_obs_ch,edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())

        # ------------ simulated Tb ----------------------
        if len(ch_num) == 2:
            for ich in d_simu['Ch_x'].tolist():
                if ich == ch_num[input_it]:
                    Ch_idx = d_simu['Ch_x'].tolist().index(ich)
        else:
            for ich in d_simu['Ch_x'].tolist():
                if ich == ch_num[0]:
                    Ch_idx = d_simu['Ch_x'].tolist().index(ich)      
        
        if sensor == 'gmi_gpm':
            Lat_x_ch = d_simu['Lat_x'][Ch_idx,:,:].flatten()
            Lon_x_ch = d_simu['Lon_x'][Ch_idx,:,:].flatten()
            Yb_x_ch = d_simu['Yb_x'][Ch_idx,:,:].flatten()
            Ya_x_ch = d_simu['Ya_x'][Ch_idx,:,:].flatten()
           
            # interpolate simulated Tbs to obs physical space
            #Yb_obspace, Ya_obspace = matlab_Interpolant(Lat_obs_ch, Lon_obs_ch, Yo_obs_ch, Lat_x_ch, Lon_x_ch, Yb_x_ch, Ya_x_ch)
            mYb_obspace = eng.griddata(matlab.double(Lon_x_ch.tolist()), matlab.double(Lat_x_ch.tolist()), matlab.double(Yb_x_ch.tolist()), matlab.double(Lon_obs_ch.tolist()), matlab.double(Lat_obs_ch.tolist()))
            Yb_obspace = np.array(mYb_obspace._data)
            mYa_obspace = eng.griddata(matlab.double(Lon_x_ch.tolist()), matlab.double(Lat_x_ch.tolist()), matlab.double(Ya_x_ch.tolist()), matlab.double(Lon_obs_ch.tolist()), matlab.double(Lat_obs_ch.tolist()))
            Ya_obspace = np.array(mYa_obspace._data)
            #if d_simu['Ch_x'][Ch_idx] == d_lowf[sensor]:
            #    is_ocean = globe.is_ocean(Lat_x_ch.flatten(), Lon_x_ch.flatten())  
            #    mask_x = is_ocean
            #else:
            #    mask_x = np.full((np.size(Lat_x_ch.flatten()), ), True)
      
            ax[i,1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            ax[i,1].coastlines(resolution='10m', color='black',linewidth=0.5)
            #ax[i,1].scatter(Lon_obs_ch[mask_x], Lat_obs_ch.flatten()[mask_x],2.5,c=Yb_obspace[mask_x],\
            #    edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree()) 
            ax[i,1].scatter(Lon_obs_ch, Lat_obs_ch,2.5,c=Yb_obspace,\
                 edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree()) 
            if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
                ax[i,1].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())        
                #ax[i,1].scatter(btk_lon, btk_lat, s=3, marker='o', facecolors='none', edgecolors='black', transform=ccrs.PlateCarree())

            ax[i,2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            ax[i,2].coastlines(resolution='10m', color='black',linewidth=0.5)
            #cs = ax[i,2].scatter(Lon_obs_ch[mask_x], Lat_obs_ch[mask_x],2.5,c=Ya_obspace[mask_x],\
            #    edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
            cs = ax[i,2].scatter(Lon_obs_ch, Lat_obs_ch,2.5,c=Ya_obspace,\
                 edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
            if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
                ax[i,2].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
                #ax[i,2].scatter(btk_lon, btk_lat, s=3, marker='o', facecolors='none', edgecolors='black', transform=ccrs.PlateCarree())
        else:
            Lat_x_ch = d_simu['Lat_x'][:,:].flatten()
            Lon_x_ch = d_simu['Lon_x'][:,:].flatten()
            Yb_x_ch = d_simu['Yb_x'][Ch_idx,:,:].flatten()
            Ya_x_ch = d_simu['Ya_x'][Ch_idx,:,:].flatten()

            # interpolate simulated Tbs to obs physical space
            #Yb_obspace, Ya_obspace = matlab_Interpolant(Lat_obs_ch, Lon_obs_ch, Yo_obs_ch, Lat_x_ch, Lon_x_ch, Yb_x_ch, Ya_x_ch)
            mYb_obspace = eng.griddata(matlab.double(Lon_x_ch.tolist()), matlab.double(Lat_x_ch.tolist()), matlab.double(Yb_x_ch.tolist()), matlab.double(Lon_obs_ch.tolist()), matlab.double(Lat_obs_ch.tolist()))
            Yb_obspace = np.array(mYb_obspace._data)
            mYa_obspace = eng.griddata(matlab.double(Lon_x_ch.tolist()), matlab.double(Lat_x_ch.tolist()), matlab.double(Ya_x_ch.tolist()), matlab.double(Lon_obs_ch.tolist()), matlab.double(Lat_obs_ch.tolist())) 
            Ya_obspace = np.array(mYa_obspace._data)
            #if d_simu['Ch_x'][Ch_idx] == d_lowf[sensor]:
            #    is_ocean = globe.is_ocean(d_simu['Lat_x'].flatten(), d_simu['Lon_x'].flatten())
            #    mask_x = is_ocean
            #else:
            #    mask_x = np.full((np.size(d_simu['Lon_x'].flatten()), ), True)
      
            ax[i,1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            ax[i,1].coastlines(resolution='10m', color='black',linewidth=0.5)
            #ax[i,1].scatter(Lon_obs_ch[mask_x], Lat_obs_ch[mask_x],2.5,c=Yb_obspace[mask_x],\
            #    edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
            ax[i,1].scatter(Lon_obs_ch, Lat_obs_ch,2.5,c=Yb_obspace,\
                 edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
            if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
                ax[i,1].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
                #ax[i,1].scatter(btk_lon, btk_lat, s=3, marker='o', facecolors='none', edgecolors='black', transform=ccrs.PlateCarree())

            ax[i,2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            ax[i,2].coastlines(resolution='10m', color='black',linewidth=0.5)
            #cs = ax[i,2].scatter(Lon_obs_ch[mask_x], Lat_obs_ch[mask_x],2.5,c=Ya_obspace[mask_x],\
            #    edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
            cs = ax[i,2].scatter(Lon_obs_ch, Lat_obs_ch,2.5,c=Ya_obspace,\
                 edgecolors='none', cmap=MWJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
            if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
                ax[i,2].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
                #ax[i,2].scatter(btk_lon, btk_lat, s=3, marker='o', facecolors='none', edgecolors='black', transform=ccrs.PlateCarree())
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
    lon_ticks = list(range(math.ceil(lon_min), math.ceil(lon_max),2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),2)) 
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
    

if __name__ == '__main__':

    Storm = 'MARIA'
    Exper_name = 'newWRF_MW_THO'

    MW_times = ['201709161800','201709161900','201709162100','201709162200','201709162300','201709170100','201709170400','201709170500','201709170700','201709170900','201709171000','201709171100','201709171300','201709171700']

    eng = matlab.engine.start_matlab() # start a new matlab process
    for DAtime in MW_times:
        obs_file_name = 'microwave_' + DAtime + '_so'
        dict_ss_ch = getSensor_Ch( '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/Obs_y/MW/'+obs_file_name )
        
        for sensor in dict_ss_ch:
            print(sensor)
            print(dict_ss_ch[sensor])
            plot_Tb( Storm, Exper_name, DAtime, sensor ) 

    eng.quit()


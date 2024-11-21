import os
import glob
import numpy as np
import netCDF4 as nc
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as patches
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from global_land_mask import globe
import math
from datetime import datetime, timedelta
import time
from scipy.interpolate import interp2d

import Util_data as UD
import Util_Vis
import matlab.engine
import Diagnostics as Diag


def RMSE(simu, obs):
    return np.sqrt( ((simu - obs) ** 2).mean() )

def Bias(simu, obs):
    return  np.sum((simu - obs),0)/np.size(obs,0)

def mean_Yo_Hx(simu, obs):
    return  np.sum((obs - simu),0)/np.size(obs,0)

# ------------------------------------------------------------------------------------------------------
#           Object: Tbs and their attributes; Operation: Read and Process 
# ------------------------------------------------------------------------------------------------------

# Read thinned IR data from SO files
def read_obs(obs_file, d_wrf_d03):

    # Define the wrf domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    with open(obs_file) as tmp:
        obs_all = tmp.readlines()

    Ch_obs = []
    Lat_obs = []
    Lon_obs = []
    Yo_obs = []
    for line in obs_all:
        line_split = line.split()

        read_lat = float(line_split[3])
        read_lon = float(line_split[4])
        if read_lat <= lat_max and read_lat >= lat_min and read_lon <= lon_max and read_lon >= lon_min:      
            Ch_obs.append( line_split[2])
            Lat_obs.append( read_lat )
            Lon_obs.append( read_lon )
            Yo_obs.append(float(line_split[5]))

    Ch_obs = np.array(Ch_obs)
    Lat_obs = np.array(Lat_obs)
    Lon_obs = np.array(Lon_obs)
    Yo_obs = np.array(Yo_obs)
    # Note: If you want to build up your matrix one column at a time,
    # you might be best off to keep it in a list until it is finished, and only then convert it into an array.    

    dict_obs_all = {'Ch_obs': Ch_obs, 'Lat_obs': Lat_obs, 'Lon_obs': Lon_obs, 'Yo_obs': Yo_obs}
    return dict_obs_all


# Average the ensemble of H(x)
def write_mean_bin ( DAtime, sensor, Hx_dir, ch_list):

    print("Initiate the function to average the ensemble...")
    start_time=time.process_time()

    # Number of ensemble members
    num_ens = 60

    # Dimension of the domain
    xmax = 297
    ymax = 297

    # List the Yb and Ya files
    file_yb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem0*.bin') )
    file_ya = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_output_mem0*.bin') )

    # Sanity check the number of calculated ens
    if np.size(file_yb) != num_ens:
        raise ValueError('Wait until all of the ens is calculated!')
    if np.size(file_ya) != num_ens:
        raise ValueError('Wait until all of the ens is calculated!')

    # Read attributes from a member
    tmp_control = np.fromfile( file_yb[0],dtype='<f4') # <: little endian; f: float; 4: 4 bytes
    n_ch = len(tmp_control)/(xmax*ymax) - 2
    n_ch = int(n_ch)
    if n_ch != len(ch_list):
        print('Error!! # of channels in data is '+str(n_ch))
    tmp_data = tmp_control[:].reshape(n_ch+2,ymax,xmax) 
    Lon_all =  [ "{0:.3f}".format(item) for item in tmp_data[0,:,:].flatten() ] #longitude
    Lat_all = [ "{0:.3f}".format(item) for item in tmp_data[1,:,:].flatten() ] #latitude
    Chnum_all = [ ch_list[0] for i in range( len(Lon_all))  ] 
 
    # ---- Read prior Tb from the ens ----
    sum_yb = np.zeros( shape=np.shape(Lon_all) )
    # Iterate thru input ens
    for ifile in file_yb:
        print('Reading the file: ' + ifile)
        tmp = np.fromfile( ifile,dtype='<f4')
        
        # Sanity check
        n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
        tmp_data = tmp[:].reshape(n_ch+2,ymax,xmax)
        sum_yb = sum_yb + tmp_data[2,:,:].flatten()

    Yb_all_mean = np.round( (sum_yb / num_ens), 3 )

    # ---- Read posterior Tb from the ens -----
    sum_ya = np.zeros( shape=np.shape(Lon_all) )
    # Iterate thru input ens
    for ifile in file_ya:
        print('Reading the file: ' + ifile)
        tmp = np.fromfile( ifile,dtype='<f4')

        # Sanity check
        n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
        tmp_data = tmp[:].reshape(n_ch+2,ymax,xmax)
        sum_ya = sum_ya + tmp_data[2,:,:].flatten()

    Ya_all_mean = np.round( (sum_ya / num_ens), 3 )

    # Stack each list into an array
    all_attrs = np.column_stack( (Lat_all, Lon_all, Chnum_all, Yb_all_mean, Ya_all_mean ) )

    # ---- Write to file and save it to the disk ----
    header = ['Lat','Lon','Ch_num','Tb_Yb_ms','Tb_Ya_ms']
    file_name = Hx_dir + "/mean_model_res_d03_" + DAtime + '_' +  sensor + '.txt'
    with open(file_name,'w') as f:
        # Add header 
        f.write('\t'.join( item.rjust(5) for item in header ) + '\n' )
        # Write the record to the file serially
        len_records = np.shape( all_attrs )[0]
        for irow in range( len_records ):
            irecord =  [str(item) for item in all_attrs[irow,:] ]
            f.write('\t'.join( item.rjust(5) for item in irecord ) + '\n')

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')


# Read the ensemble mean of Hx at model resolution 
def read_Hx_mean( Hx_file ):

    lat_x = []
    lon_x = []
    ch_x = []
    meanYb = []
    meanYa = []

    with open(Hx_file) as f:
        next(f) 
        all_lines = f.readlines()

    for line in all_lines:
        split_line = line.split()
        lat_x.append( float(split_line[0]) )
        lon_x.append( float(split_line[1]) )
        ch_x.append( split_line[2] )
        meanYb.append( float(split_line[3]) )
        meanYa.append( float(split_line[4]) )

    lat_x = np.array( lat_x )
    lon_x = np.array( lon_x )
    ch_x = np.array( ch_x )
    meanYb = np.array( meanYb )
    meanYa = np.array( meanYa )

    dict_Tb_all = {'Lat_x':lat_x, 'Lon_x':lon_x, 'Ch_x': ch_x, 'Yb_x':meanYb, 'Ya_x':meanYa}
    return dict_Tb_all


# Interpolate IR Tbs in model resolution to obs locations and Write it to a txt file
def interp_simu_to_obs_matlab( Hx_dir, sensor, DAtime, d_obs):

    print("Initiate the function to interpolate simulated Tbs to obs locations...")
    start_time=time.process_time()

    # Read simulated Tbs
    Hx_file = Hx_dir + "/mean_model_res_d03_" + DAtime + '_' +  sensor + '.txt' 
    print('Reading ensemble mean of Hx: ', Hx_file,'...')
    d_simu = read_Hx_mean( Hx_file )

    # Initiate the container for interpolated obs
    Ch_obspace_all = []
    Lat_obspace_all = []
    Lon_obspace_all =[]
    Yb_obspace_all = []
    Ya_obspace_all = []
    Yo_all = []

    # Start a matlab process
    eng = matlab.engine.start_matlab()

    for ich in ch_list:

        print('Channel number: ', ich)
        #Ch_idx_obs = d_obs['Ch'] == ich
        Lat_obs_ch = d_obs['lat']
        Lon_obs_ch = d_obs['lon']
        Yo_obs_ch = d_obs['obs']

        Ch_idx_x = d_simu['Ch_x'] == ich
        Lat_x_ch = d_simu['Lat_x'][Ch_idx_x]
        Lon_x_ch = d_simu['Lon_x'][Ch_idx_x]
        Yb_x_ch = d_simu['Yb_x'][Ch_idx_x]
        Ya_x_ch = d_simu['Ya_x'][Ch_idx_x]

        print('Number of NaN in Yb_x_ch', sum(np.isnan(Yb_x_ch)))
        print('Number of NaN in Ya_x_ch', sum(np.isnan(Ya_x_ch)))
        # interpolate simulated Tbs to obs location
        mYb_obspace = eng.griddata(matlab.double(Lon_x_ch.tolist()), matlab.double(Lat_x_ch.tolist()), matlab.double(Yb_x_ch.tolist()), matlab.double(Lon_obs_ch.tolist()), matlab.double(Lat_obs_ch.tolist()) )
        Yb_obspace = mYb_obspace._data
        print('Number of NaN in Yb_obspace', sum(np.isnan(Yb_obspace)))
        mYa_obspace = eng.griddata(matlab.double(Lon_x_ch.tolist()), matlab.double(Lat_x_ch.tolist()), matlab.double(Ya_x_ch.tolist()), matlab.double(Lon_obs_ch.tolist()), matlab.double(Lat_obs_ch.tolist()) )
        Ya_obspace = mYa_obspace._data
        print('Number of NaN in Ya_obspace', sum(np.isnan(Ya_obspace)))

        # Add values to the container
        Ch_obspace_ch = np.full( np.shape(Yo_obs_ch), ich )
        Ch_obspace_all = Ch_obspace_all + [ item  for item in Ch_obspace_ch]
        Lat_obspace_all =  Lat_obspace_all + ["{0:.3f}".format(item) for item in Lat_obs_ch]
        Lon_obspace_all = Lon_obspace_all + ["{0:.3f}".format(item) for item in Lon_obs_ch]
        Yb_obspace_all =  Yb_obspace_all + ["{0:.3f}".format(item) for item in Yb_obspace]
        Ya_obspace_all =  Ya_obspace_all + ["{0:.3f}".format(item) for item in Ya_obspace]
        Yo_all = Yo_all + ["{0:.3f}".format(item) for item in Yo_obs_ch]


    # End the matlab process
    eng.quit()

    # Stack each list into an array
    all_attrs = np.column_stack( (Lat_obspace_all, Lon_obspace_all, Ch_obspace_all, Yo_all, Yb_obspace_all, Ya_obspace_all) )

    # ---- Write to file and save it to the disk ----
    header = ['Lat','Lon','Ch_num','Tb_obs','Tb_Yb_obs','Tb_Ya_obs']
    file_name = Hx_dir + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    with open(file_name,'w') as f:
        # Add header 
        f.write('\t'.join( item.rjust(6) for item in header ) + '\n' )
        # Write the record to the file serially
        len_records = np.shape( all_attrs )[0]
        for irow in range( len_records ):
            irecord =  [str(item) for item in all_attrs[irow,:] ]
            f.write('\t'.join( item.rjust(6) for item in irecord ) + '\n')

    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')

    return( len_records )

# Read variables at obs resolution/location
def read_Tb_obsRes(Tb_file, sensor ):

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

    #if np.size(lat_obs) != dict_ss_len[sensor]:
    #    raise ValueError('The length of post-processed file is not equal to the pre-processed file!')

    lat_obs = np.array( lat_obs )
    lon_obs = np.array( lon_obs )
    ch_obs = np.array( ch_obs )
    Yo_obs = np.array( Yo_obs )
    meanYb_obs = np.array( meanYb_obs )
    meanYa_obs = np.array( meanYa_obs )
    print('Number of NaN in meanYa_obs', sum(np.isnan(meanYa_obs)))

    dict_Tb_all = {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'ch_obs':ch_obs, 'Yo_obs':Yo_obs, 'meanYb_obs':meanYb_obs, 'meanYa_obs':meanYa_obs}
    return dict_Tb_all

# Read variables at model resolution/location
def read_Tb_modelRes(Tb_file, sensor ):

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
        meanYb_obs.append( float(split_line[3]) )
        meanYa_obs.append( float(split_line[4]) )

    #if np.size(lat_obs) != dict_ss_len[sensor]:
    #    raise ValueError('The length of post-processed file is not equal to the pre-processed file!')

    lat_obs = np.array( lat_obs )
    lon_obs = np.array( lon_obs )
    ch_obs = np.array( ch_obs )
    meanYb_obs = np.array( meanYb_obs )
    meanYa_obs = np.array( meanYa_obs )
    print('Number of NaN in meanYa_model', sum(np.isnan(meanYa_obs)))

    dict_Tb_all = {'lat_model':lat_obs, 'lon_model':lon_obs, 'ch_model':ch_obs,'meanYb_model':meanYb_obs, 'meanYa_model':meanYa_obs}
    return dict_Tb_all

# ------------------------------------------------------------------------------------------------------
#           Object: things that might be of interest to the plotting  
# ------------------------------------------------------------------------------------------------------

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


def plot_Tb(Storm, Exper_name, Hx_dir, DAtime, sensor, ch_list ):

    # Read WRF domain
    wrf_file = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = read_wrf_domain( wrf_file )

    # Read Tbs of obs, Hxb, Hxa
    Tb_file = Hx_dir + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    d_all = read_Tb_obsRes(Tb_file, sensor )  

    # impose other conditions
    if limit:
        condi1 = d_all['meanYb_obs'] - d_all['Yo_obs'] < -3 #d_all['Yo_obs'] <= 210
        idx_a = np.where( condi1 )[0]
        condi2 = d_all['Yo_obs'] <= 220 #d_all['meanYb_obs'] <= 210
        idx_b = np.where( condi2 )[0]
        idx_x = list(set(idx_a)&set(idx_b))
    else:
        idx_x = range(len(d_all['Yo_obs']))

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir, Storm, DAtime )

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

    ax[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    ax[0].scatter(d_all['lon_obs'],d_all['lat_obs'],1.5,c=d_all['Yo_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())

    ax[1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    ax[1].scatter(d_all['lon_obs'][idx_x], d_all['lat_obs'][idx_x],1.5,c=d_all['meanYb_obs'][idx_x],\
                edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    #if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
    #    ax[1].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())

    ax[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax[2].scatter(d_all['lon_obs'][idx_x], d_all['lat_obs'][idx_x],1.5,c=d_all['meanYa_obs'][idx_x],\
                edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    #if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
    #    ax[2].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())

    # Colorbar
    caxes = f.add_axes([0.2, 0.1, 0.6, 0.02])
    cbar = f.colorbar(cs, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)

    #subplot title
    matplotlib.rcParams['mathtext.fontset'] = 'custom'
    matplotlib.rcParams['mathtext.bf'] = 'STIXGeneral:italic:bold'
    font = {'size':9,}
    ax[0].set_title('Yo', font, fontweight='bold')
    ax[1].set_title(r'$\mathbf{\overline{H(Xb)}}$', font, )
    ax[2].set_title(r'$\mathbf{\overline{H(Xa)}}$', font)

    #title for all
    if not limit:
        f.suptitle(Storm+': '+Exper_name+'\n'+DAtime, fontsize=6, fontweight='bold')
    else:
        title_name = Storm+': '+Exper_name+' '+DAtime+'\n H(Xb)-obs < -3K (obs<=220k)'
        f.suptitle(title_name, fontsize=6, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(3):
        gl = ax[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
       
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
        gl.xlabel_style = {'size': 4}
        gl.ylabel_style = {'size': 6}

    if not limit:
        des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IRch8_Obspace/'+DAtime+'_'+sensor+'_Obspace.png'
    else:
        des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IRch8_Obspace/'+DAtime+'_'+sensor+'_Obspace_limit.png'
    plt.savefig(des_name,dpi=300)
    print('Saving the figure: ',des_name)


def plot_Tb_diff(Storm, Exper_name, Hx_dir, DAtime, sensor, ch_list ):

    # Read WRF domain
    wrf_file = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = read_wrf_domain( wrf_file )

    # Read Tbs of obs, Hxb, Hxa
    Tb_file = Hx_dir + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
    d_all = read_Tb_obsRes(Tb_file, sensor )

    # impose other conditions
    if limit:
        condi1 = abs(d_all['meanYb_obs'] - d_all['Yo_obs']) <= 3 #d_all['Yo_obs'] <= 210
        idx_a = np.where( condi1 )[0]
        condi2 = d_all['Yo_obs'] <= 210 #d_all['meanYb_obs'] <= 210
        idx_b = np.where( condi2 )[0]
        idx_x = list(set(idx_a)&set(idx_b))
    else:
        idx_x = range(len(d_all['Yo_obs']))

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir, Storm, DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # Prepare to calculate Bias and RMSE between the Hx and Yo
    metric = np.full((2,2), fill_value=None) # default type: none

    # ------------------ Plot -----------------------
    fig, axs=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,2.5), dpi=400)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']
    
    # Set the map
    for i in range(3):
        axs.flat[i].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        axs.flat[i].coastlines(resolution='10m', color='black',linewidth=0.5)
    
    # Obs
    min_obs = 185
    max_obs = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )
    obs_s = axs.flat[0].scatter(d_all['lon_obs'],d_all['lat_obs'],1.5,c=d_all['Yo_obs'],edgecolors='none', cmap=IRcmap, vmin=min_obs, vmax=max_obs,transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        axs.flat[0].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='darkviolet', transform=ccrs.PlateCarree())
    #axs.flat[0].add_patch(patches.Polygon(path,facecolor='none',edgecolor='white',linewidth=0.5 ))
    # Colorbar
    caxes = fig.add_axes([0.12, 0.1, 0.25, 0.02])
    obs_bar = fig.colorbar(obs_s,ax=axs[0],orientation="horizontal", cax=caxes)
    obs_bar.ax.tick_params(labelsize=6)

    max_T=15
    min_T=-15
    # HXb - Obs & HXa - Obs
    if plot_scatter:
        xb_s = axs.flat[1].scatter(d_all['lon_obs'], d_all['lat_obs'],1.5,c=d_all['meanYb_obs']-d_all['Yo_obs'],\
                        edgecolors='none', cmap='bwr', vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree() ) #RdBu_r
        xa_s = axs.flat[2].scatter(d_all['lon_obs'], d_all['lat_obs'],1.5,c=d_all['meanYa_obs']-d_all['Yo_obs'],\
                        edgecolors='none', cmap='bwr', vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    else:
        bounds = np.linspace(min_T,max_T,7)
        ncdir = nc.Dataset( wrf_file )
        xb_s = axs.flat[1].tricontourf(d_all['lon_obs'], d_all['lat_obs'], d_all['meanYb_obs']-d_all['Yo_obs'], cmap='bwr', \
                        vmin=min_T, vmax=max_T, levels=bounds, transform=ccrs.PlateCarree(), extend='both' )
        xa_s = axs.flat[2].tricontourf(d_all['lon_obs'], d_all['lat_obs'], d_all['meanYa_obs']-d_all['Yo_obs'], cmap='bwr', \
                        vmin=min_T, vmax=max_T, levels=bounds, transform=ccrs.PlateCarree(), extend='both' )
    # add a reference point
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        axs.flat[1].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
        axs.flat[2].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
    # add a box 
    #axs.flat[1].add_patch(patches.Polygon(path,facecolor='none',edgecolor='black',linewidth=0.5 )) 
    #axs.flat[2].add_patch(patches.Polygon(path,facecolor='none',edgecolor='black',linewidth=0.5 ))

    # Colorbar
    caxes = fig.add_axes([0.4, 0.1, 0.5, 0.02])
    cb_diff_ticks = np.linspace(min_T, max_T, 7, endpoint=True)
    cbar = fig.colorbar(xa_s, ax=axs[1:], ticks=cb_diff_ticks, orientation="horizontal", cax=caxes, extend='both')
    cbar.ax.tick_params(labelsize=6)
    
    # Calculate metrics and annotate them
    metric[0,0] = Bias(d_all['meanYb_obs'], d_all['Yo_obs'] )
    metric[0,1] = Bias(d_all['meanYa_obs'], d_all['Yo_obs'] )
    metric[1,0] = RMSE(d_all['meanYb_obs'], d_all['Yo_obs'] )
    metric[1,1] = RMSE(d_all['meanYa_obs'], d_all['Yo_obs'] ) 

    fig.text( 0.42,0.15,'Bias:'+'%.2f' % metric[0,0],rotation='horizontal',fontsize=6,fontweight='bold')
    fig.text( 0.52,0.15,'; RMSE:'+'%.2f' % metric[1,0],rotation='horizontal',fontsize=6,fontweight='bold')
    fig.text( 0.66,0.15,'Bias:'+'%.2f' % metric[0,1],rotation='horizontal',fontsize=6,fontweight='bold')
    fig.text( 0.76,0.15,'; RMSE:'+'%.2f' % metric[1,1],rotation='horizontal',fontsize=6,fontweight='bold')

    #subplot title
    font = {'size':8,}
    axs.flat[0].set_title('CH'+ch_list[0]+': Yo', font, fontweight='bold')
    axs.flat[1].set_title(r'$\mathbf{\overline{H(Xb)}}$'+'-Yo', font, fontweight='bold')
    axs.flat[2].set_title(r'$\mathbf{\overline{H(Xa)}}$'+'-Yo', font, fontweight='bold')

    #title for all
    fig.suptitle(Storm+': '+Exper_name+'\n'+DAtime, fontsize=6, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(3):
        gl = axs.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5,alpha=0.7,color='gray',linestyle='--')

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
        gl.xlabel_style = {'size': 4}
        gl.ylabel_style = {'size': 6}

    if plot_scatter:
        figure_des=small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IRch8_Obspace_Diff/'+DAtime+'_'+sensor+'_Obspace_Diff_scatter.png'
    else:
        figure_des=small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IRch8_Obspace_Diff/'+DAtime+'_'+sensor+'_Obspace_Diff_contourf.png'
    plt.savefig(figure_des, dpi=400)
    print('Saving the figure: ', figure_des)
    plt.close()

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'MARIA'
    DA = 'CONV'
    MP = 'WSM6'
 
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709160000'
    end_time_str = '201709160000'
    Consecutive_times = True
    
    # limitations
    limit = False

    Interp_to_obs = True
    plot_full = False
    plot_diff = False
    plot_scatter = False
    # ------------------------------------------------------   


    # Create experiment names

    Exper_name = UD.generate_one_name( Storm,DA,MP )
    Exper_obs =  UD.generate_one_name( Storm,'IR',MP )
    
    if not Consecutive_times:
        IR_times = ['201708230400','201708231200']
        #['201709051200','201709051800','201709060000','201709060600','201709061200','201709061800','201709070000']
        #['201709030000','201709030600','201709031200','201709031800','201709040000','201709040600','201709041200','201709041800','201709050000']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Process Tbs to obs location
    if Interp_to_obs:
        for DAtime in IR_times:
            # Read assimilated obs 
            file_Diag = big_dir+Storm+'/'+Exper_obs+'/run/'+DAtime+'/enkf/d03/fort.10000'
            d_obs = Diag.Find_IR( file_Diag, fort_v )

            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            # Calculate mean of H(Xb) and H(Xa)
            print('------------ Calculate mean of Hx in model resolution --------------')
            write_mean_bin( DAtime, sensor, Hx_dir, ch_list )

            # Interpolate HX in model resolution to obs location AND write it to a txt file
            print('------------ Interpolate Hx in model resolution to obs location --------------')
            interp_simu_to_obs_matlab( Hx_dir, sensor, DAtime, d_obs)
            #time.sleep(60)

    # Plot Tb
    if plot_full:
        # Create plot dir
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IRch8_Obspace/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

        start_time=time.process_time()
        for DAtime in IR_times:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            print('------------ Plot ----------------------')         
            print('DAtime: '+ DAtime)
            plot_Tb( Storm, Exper_name, Hx_dir, DAtime, sensor, ch_list) 
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')

    # Plot Tb diff
    if plot_diff:
        # Create plot dir
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IRch8_Obspace_Diff/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

        start_time=time.process_time()
        for DAtime in IR_times:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            print('------------ Plot ----------------------')
            print('DAtime: '+ DAtime)
            plot_Tb_diff( Storm, Exper_name, Hx_dir, DAtime, sensor, ch_list)
        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')







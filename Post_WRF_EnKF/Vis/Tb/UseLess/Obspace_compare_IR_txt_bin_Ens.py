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

def RMSE(simu, obs):
    return np.sqrt( ((simu - obs) ** 2).mean() )

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


# Read crtm calcuated IR data from binary file
def read_simu_Tb(Hxb_file, Hxa_file, ch_list):
    
    xmax = 297
    ymax = 297
    print('Hxb_file:' + Hxb_file)
    print('Hxa_file:' + Hxa_file)
    print('xmax, ymax: '+str(xmax)+' '+str(ymax))

    Hxb_data = np.fromfile(Hxb_file,dtype='<f4') # <: little endian; f: float; 4: 4 bytes
    n_ch = len(Hxb_data)/(xmax*ymax) - 2
    n_ch = int(n_ch)
    if n_ch != len(ch_list):
        print('Error!! # of channels in data is '+str(n_ch))
    Hxb_sim = Hxb_data[:].reshape(n_ch+2,ymax,xmax)

    dict_simu_Tb = {}
    dict_simu_Tb['Lon_x'] = Hxb_sim[0,:,:].flatten()
    dict_simu_Tb['Lat_x'] = Hxb_sim[1,:,:].flatten()
    dict_simu_Tb['Ch_x'] = ch_list[0]
    dict_simu_Tb['Yb_x'] = Hxb_sim[2,:,:].flatten()

    Hxa_data = np.fromfile(Hxa_file,dtype='<f4')
    n_ch = len(Hxa_data)/(xmax*ymax) - 2
    n_ch = int(n_ch)
    if n_ch != len(ch_list):
        print('Error!! # of channels in data is '+str(n_ch))
    Hxa_sim = Hxa_data[:].reshape(n_ch+2,ymax,xmax)
    dict_simu_Tb['Ya_x'] = Hxa_sim[2,:,:].flatten() 
    
    return dict_simu_Tb


# Average the ensemble of H(x)
def write_mean_bin ( DAtime, sensor, Hx_dir, ch_list ):

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
        tmp = np.fromfile( file_yb[0],dtype='<f4')
        
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
        tmp = np.fromfile( file_ya[0],dtype='<f4')

        # Sanity check
        n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
        tmp_data = tmp[:].reshape(n_ch+2,ymax,xmax)
        sum_ya = sum_ya + tmp_data[2,:,:].flatten()

    Ya_all_mean = np.round( (sum_ya / num_ens), 3 )

    # Stack each list into an array
    all_attrs = np.column_stack( (Lat_all, Lon_all, Chnum_all, Yb_all_mean, Ya_all_mean ) )

    # ---- Write to file and save it to the disk ----
    header = ['Lat','Lon','Ch_num','Tb_Yb','Tb_Ya']
    file_name = Hx_dir + "/mean_model_res_d03" + DAtime + '_' +  sensor + '.txt'
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
def interp_simu_to_obs_matlab( Hxb, sensor, DAtime, obs_file_dir, d_wrf_d03 ):

    print("Initiate the function to interpolate simulated Tbs to obs locations...")
    start_time=time.process_time()

    # Read obs data
    d_obs = read_obs( obs_file_dir,d_wrf_d03 )
    
    # Load simulated Tbs
    Hxa = Hxb.replace( 'input','output') 
    print('Reading the Hxb: ' + Hxb)
    d_simu = read_simu_Tb(Hxb, Hxa, ch_list)

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
        Ch_idx_obs = d_obs['Ch_obs'] == ich
        Lat_obs_ch = d_obs['Lat_obs'][Ch_idx_obs]
        Lon_obs_ch = d_obs['Lon_obs'][Ch_idx_obs]
        Yo_obs_ch = d_obs['Yo_obs'][Ch_idx_obs]

        Ch_idx_x = d_simu['Ch_x'] == ich
        Lat_x_ch = d_simu['Lat_x'][Ch_idx_x]
        Lon_x_ch = d_simu['Lon_x'][Ch_idx_x]
        Yb_x_ch = d_simu['Yb_x'][Ch_idx_x]
        Ya_x_ch = d_simu['Ya_x'][Ch_idx_x]

        print('Number of NaN in Yb_x_ch', sum(np.isnan(Yb_x_ch)))
        print('Number of NaN in Ya_x_ch', sum(np.isnan(Ya_x_ch)))
        # interpolate simulated Tbs to obs location
        mYb_obspace = eng.griddata(matlab.double(Lon_x_ch.tolist()), matlab.double(Lat_x_ch.tolist()), matlab.double(Yb_x_ch.tolist()), matlab.double(Lon_obs_ch.tolist()), matlab.double(Lat_obs_ch.tolist()) )
        Yb_obspace = np.array(mYb_obspace._data)
        print('Number of NaN in Yb_obspace', sum(np.isnan(Yb_obspace)))
        mYa_obspace = eng.griddata(matlab.double(Lon_x_ch.tolist()), matlab.double(Lat_x_ch.tolist()), matlab.double(Ya_x_ch.tolist()), matlab.double(Lon_obs_ch.tolist()), matlab.double(Lat_obs_ch.tolist()) )
        Ya_obspace = np.array(mYa_obspace._data)
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
    file_name = Hxb.replace('TB_GOES_CRTM_input','Interp_Tb')
    file_name = file_name.replace('bin','txt')
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
def read_allTb(Tb_file, sensor ):

    lat_obs = []
    lon_obs = []
    ch_obs = []
    Yo_obs = []
    meanYb_obs = []
    meanYa_obs = []

    # number of columns of each record
    #ncol = 9

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


# ------------------------------------------------------------------------------------------------------
#           Object: things that might be of interest to the plotting  
# ------------------------------------------------------------------------------------------------------

# Storm center produced by ATCF of the NWP
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

# Plotting domain
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


def plot_Tb(Storm, Exper_name, ifile, DAtime, sensor, ch_list ):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = read_wrf_domain( wrf_file )

    # Read Tbs of obs, Hxb, Hxa
    d_all = read_allTb(ifile, sensor )  

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = read_TCvitals('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

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
    ax[1].scatter(d_all['lon_obs'], d_all['lat_obs'],1.5,c=d_all['meanYb_obs'],\
                edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax[1].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())

    ax[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax[2].scatter(d_all['lon_obs'], d_all['lat_obs'],1.5,c=d_all['meanYa_obs'],\
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

    head_tail = os.path.split( ifile )
    mem = head_tail[1].replace('_d03_2017-08-22_12:00.txt','')
    plt.savefig('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IR/obspace_60mem/'+DAtime+'_'+sensor+'_'+mem+'.png', dpi=300)
    #print('Saving the figure: ', '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IR/Obspace/'+DAtime+'_'+sensor+'_Obspace.png') 


def plot_Tb_diff(Storm, Exper_name, Hx_dir, DAtime, sensor, ch_list ):

    # Read WRF domain
    wrf_file = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = read_wrf_domain( wrf_file )

    # Read Tbs of obs, Hxb, Hxa
    Tb_file = Hx_dir + "/mean_obs_res_d03" + DAtime + '_' +  sensor + '.txt'
    d_all = read_allTb(Tb_file, sensor )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat = read_TCvitals('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/TCvitals/'+Storm+'_tcvitals', DAtime)
        print( 'Location from TCvital: ', tc_lon, tc_lat )

    # Prepare to calculate RMSE between the Hx and Yo
    rmse = np.full(2, fill_value=None) # default type: none

    # ------------------ Plot -----------------------
    f, ax=plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,3), dpi=400)

    # Customize colormap
    max_T=50
    min_T=-50
    #min_RWB = 0
    #newRWB = Util_Vis.newRWB(max_T, min_T, min_RWB)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    # HXb - Obs
    ax[0].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    ax[0].scatter(d_all['lon_obs'], d_all['lat_obs'],1.5,c=d_all['meanYb_obs']-d_all['Yo_obs'],\
                edgecolors='none', cmap='RdBu_r', vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax[0].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())
    
    rmse[0] = RMSE(d_all['meanYb_obs'], d_all['Yo_obs'] )

    # HXa - Obs
    ax[1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax[1].scatter(d_all['lon_obs'], d_all['lat_obs'],1.5,c=d_all['meanYa_obs']-d_all['Yo_obs'],\
                edgecolors='none', cmap='RdBu_r', vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax[1].scatter(tc_lon, tc_lat, s=3, marker='*', edgecolors='white', transform=ccrs.PlateCarree())
    
    rmse[1] = RMSE(d_all['meanYa_obs'], d_all['Yo_obs'] )
    
    # Colorbar
    cb_ticks = np.linspace(min_T, max_T, 5, endpoint=True)
    caxes = f.add_axes([0.2, 0.97, 0.6, 0.02])
    cbar = f.colorbar(cs, ticks=cb_ticks, orientation="horizontal", cax=caxes)
    cbar.ax.tick_params(labelsize=6)

    #subplot title
    font = {'size':8,}
    if rmse[0] is None:
        ax[0].set_title('Ch'+ch_list[0]+':H(Xb)-Yo ', font, fontweight='bold')
    else:
        rmse_str = '%.2f' % rmse[0]
        ax[0].set_title('Ch'+ch_list[0]+':H(Xb)-Yo '+rmse_str, font, fontweight='bold')

    if rmse[1] is None:
        ax[1].set_title('Ch'+ch_list[0]+':H(Xa)-Yo ', font, fontweight='bold')
    else:
        rmse_str = '%.2f' % rmse[1]
        ax[1].set_title('Ch'+ch_list[0]+':H(Xa)-Yo '+rmse_str, font, fontweight='bold')
 

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(2):
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

    plt.savefig('/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IR/Obspace/'+DAtime+'_'+sensor+'_Obspace_Diff_range50.png', dpi=300)
    print('Saving the figure: ', '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IR/Obspace/'+DAtime+'_'+sensor+'_Obspace_Diff.png')





if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    Storm = 'HARVEY'
    Exper_name = 'JerryRun/IR_THO'
    sensor = 'abi_gr'
    ch_list = ['8',]
    start_time_str = '201709160000' 
    end_time_str = '201709160000'
    Interp_to_obs = False
    Consecutive_times = False
    If_plot = True
    num_ens = 60

    if not Consecutive_times:
        IR_times = ['201708221200',]#'201708221800','201708230000','201708230600','201708231200']
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        IR_times = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]


    # Process Tbs to obs location
    if Interp_to_obs:
        for DAtime in IR_times:
            # Get obs's sensor/channel info
            obs_file_name = 'radiance_d03_' + DAtime + '_so'
            obs_file_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/Obs_y/IR/'+obs_file_name

            Hx_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            # List the Yb and Ya files
            file_yb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem0*.bin') ) 
            for Hxb in file_yb:
                print('Reading the Hxb: ' + Hxb)
                Hxa = Hxb.replace( 'input','output')
                # Interpolate HX in model resolution to obs location AND write it to a txt file
                print('------------ Interpolate Hx in model resolution to obs location --------------')
                # Read WRF domain
                wrf_file = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
                d_wrf_d03 = read_wrf_domain( wrf_file )
                interp_simu_to_obs_matlab( Hxb, sensor, DAtime, obs_file_dir, d_wrf_d03 )
                #time.sleep(60)

    # Plot
    if If_plot:
        for DAtime in IR_times:
            Hx_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            print('------------ Plot ----------------------')         
            print('DAtime: '+ DAtime)
            # List the interpolated files
            interp_files = sorted( glob.glob(Hx_dir + '/Interp_Tb_mem0*.txt') ) 
            for ifile in interp_files:
                print('Plotting ' + ifile)
                plot_Tb( Storm, Exper_name, ifile, DAtime, sensor, ch_list) 
            #plot_Tb_diff( Storm, Exper_name, Hx_dir, DAtime, sensor, ch_list) 
    















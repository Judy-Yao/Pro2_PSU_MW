
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


# Interpolate IR Tbs in model resolution to obs locations and Write it to a txt file
def interp_simu_to_obs_matlab( Hxb, sensor, DAtime, d_obs, ):

    print("Initiate the function to interpolate simulated Tbs to obs locations...")
    start_time=time.process_time()

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
        
        #Ch_idx_obs = d_obs['Ch_obs'] == ich
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

    dict_Tb_all = {'lat_obs':lat_obs, 'lon_obs':lon_obs, 'ch_obs':ch_obs, 'Yo_obs':Yo_obs, 'meanYb_obs':meanYb_obs, 'meanYa_obs':meanYa_obs}
    return dict_Tb_all

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

    des_path = big_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IR/obspace_60mem/'+DAtime+'_'+sensor+'_'+mem+'.png'
    plt.savefig( des_path, dpi=300 )
    print( 'Saving the figure: ', des_path )
    plt.close()

def plot_Tb_diff( Storm, Exper_name, ifile, DAtime, sensor, ch_list ):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = UD.read_wrf_domain( wrf_file )

    # Read Tbs of obs, Hxb, Hxa
    d_all = read_allTb(ifile, sensor )

    # Read location from TCvitals
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18']):
        tc_lon, tc_lat, tc_slp = UD.read_TCvitals(small_dir, Storm, DAtime )

    # Prepare to calculate Bias and RMSE between the Hx and Yo
    metric = np.full((2,2), fill_value=None) # default type: none

    # ------------------ Plot -----------------------
    f, ax=plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,2.5), dpi=300)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    # Obs
    min_T = 185
    max_T = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )
    ax[0].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax[0].coastlines(resolution='10m', color='black',linewidth=0.5)
    obs_s = ax[0].scatter(d_all['lon_obs'],d_all['lat_obs'],1.5,c=d_all['Yo_obs'],edgecolors='none', cmap=IRcmap, vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax.flat[0].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='darkviolet', transform=ccrs.PlateCarree())
    # Colorbar
    caxes = f.add_axes([0.12, 0.1, 0.25, 0.02])
    obs_bar = f.colorbar(obs_s,ax=ax[0],orientation="horizontal", cax=caxes)
    obs_bar.ax.tick_params(labelsize=6)

    max_T=15
    min_T=-15
    # HXb - Obs & HXa - Obs
    ax[1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[1].coastlines(resolution='10m', color='black',linewidth=0.5)
    xb_s = ax[1].scatter(d_all['lon_obs'],d_all['lat_obs'],1.5,c=d_all['meanYb_obs']-d_all['Yo_obs'],\
                edgecolors='none', cmap='bwr', vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())

    ax[2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[2].coastlines(resolution='10m', color='black',linewidth=0.5)
    xa_s = ax[2].scatter(d_all['lon_obs'],d_all['lat_obs'],1.5,c=d_all['meanYa_obs']-d_all['Yo_obs'],\
                edgecolors='none', cmap='bwr', vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    # add a reference point
    if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
        ax[1].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
        ax[2].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='black', transform=ccrs.PlateCarree())
    # Colorbar
    caxes = f.add_axes([0.4, 0.1, 0.5, 0.02])
    cb_diff_ticks = np.linspace(min_T, max_T, 7, endpoint=True)
    cbar = f.colorbar(xa_s, ax=ax[1:], ticks=cb_diff_ticks, orientation="horizontal", cax=caxes, extend='both')
    cbar.ax.tick_params(labelsize=6)

    # Calculate metrics and annotate them
    metric[0,0] = Bias(d_all['meanYb_obs'], d_all['Yo_obs'] )
    metric[0,1] = Bias(d_all['meanYa_obs'], d_all['Yo_obs'] )
    metric[1,0] = RMSE(d_all['meanYb_obs'], d_all['Yo_obs'] )
    metric[1,1] = RMSE(d_all['meanYa_obs'], d_all['Yo_obs'] )

    f.text( 0.42,0.15,'Bias:'+'%.2f' % metric[0,0],rotation='horizontal',fontsize=6,fontweight='bold')
    f.text( 0.52,0.15,'; RMSE:'+'%.2f' % metric[1,0],rotation='horizontal',fontsize=6,fontweight='bold')
    f.text( 0.66,0.15,'Bias:'+'%.2f' % metric[0,1],rotation='horizontal',fontsize=6,fontweight='bold')
    f.text( 0.76,0.15,'; RMSE:'+'%.2f' % metric[1,1],rotation='horizontal',fontsize=6,fontweight='bold')

    #subplot title
    font = {'size':8,}
    ax.flat[0].set_title('CH'+ch_list[0]+': Yo', font, fontweight='bold')
    ax.flat[1].set_title(r'$\mathbf{\overline{H(Xb)}}$'+'-Yo', font, fontweight='bold')
    ax.flat[2].set_title(r'$\mathbf{\overline{H(Xa)}}$'+'-Yo', font, fontweight='bold')

    #title for all
    f.suptitle(Storm+': '+Exper_name+'\n'+DAtime, fontsize=6, fontweight='bold')

    #subplot title
    head_tail = os.path.split( ifile )
    mem = head_tail[1].replace('Interp_Tb_mem','')
    mem = mem.replace('_d03_2017-09-03_00:00.txt','')

    font = {'size':8,}
    ax[0].set_title('Yo', font, fontweight='bold')
    ax[1].set_title(mem+': '+r'$H_{Tb}(Xb)$', font, fontweight='bold')
    ax[2].set_title(mem+': '+r'$H_{Tb}(Xa)$', font, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
    for j in range(3):
        gl = ax[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')

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

    des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IR_60mem_Diff/'+DAtime+'_'+mem+'_Obspace_Diff_mspace.png'
    plt.savefig( des_name, dpi=300)
    plt.close()
    print('Saving the figure: ', des_name)



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # -------- Configuration -----------------
    Storm = 'IRMA'
    DA = 'CONV'
    MP = 'THO'

    # Time range set up
    start_time_str = '201709030000'
    end_time_str = '201709030000'
    Consecutive_times = True

    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']
    num_ens = 60

    start_time_str = '201709030000' 
    end_time_str = '201709030000'
    Consecutive_times = True

    Interp_to_obs = False
    
    If_plot = False
    If_plot_diff = True
    # -----------------------------------------
    
    # Experiment name
    Exper_name = UD.generate_one_name( Storm,DA,MP )
    Exper_obs = UD.generate_one_name( Storm,'IR',MP )

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
            file_Diag = big_dir+Storm+'/'+Exper_obs+'/run/'+DAtime+'/enkf/d03/fort.10000'
            d_obs = Diag.Find_IR( file_Diag, fort_v )

            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            # List the Yb and Ya files
            Hxb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_*') ) 
            print('------------ Interpolate Hx in model resolution to obs location --------------')
            for ifile in Hxb:
                print('Reading the Hxb: ' + ifile)
                # Interpolate HX in model resolution to obs location AND write it to a txt file
                interp_simu_to_obs_matlab( ifile, sensor, DAtime, d_obs, )
                #time.sleep(60)

    # Plot
    if If_plot:
        for DAtime in IR_times:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            print('------------ Plot ----------------------')         
            print('DAtime: '+ DAtime)
            # List the interpolated files
            interp_files = sorted( glob.glob(Hx_dir + '/Interp_Tb_mem0*.txt') ) 
            for ifile in interp_files:
                print('Plotting ' + ifile)
                plot_Tb( Storm, Exper_name, ifile, DAtime, sensor, ch_list)

    if If_plot_diff:
        for DAtime in IR_times:
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            print('------------ Plot ----------------------')
            print('DAtime: '+ DAtime)
            # List the interpolated files
            interp_files = sorted( glob.glob(Hx_dir + '/Interp_Tb_mem0*.txt') )
            for ifile in interp_files:
                print('Plotting ' + ifile)
                plot_Tb_diff( Storm, Exper_name, ifile, DAtime, sensor, ch_list)

    















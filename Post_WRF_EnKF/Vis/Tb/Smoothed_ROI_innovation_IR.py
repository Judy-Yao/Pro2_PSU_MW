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
from scipy.spatial import cKDTree

import Util_data as UD
import Util_Vis
import Obspace_compare_IR_txt_bin as ROIR
import Diagnostics as Diag

def convolution_IR_obsRes( lon_obs,lat_obs,diff_obs,roi ):

    #print('diff_obs:',np.shape( diff_obs ))

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    ncdir = nc.Dataset(wrf_file, 'r')
    Lat_x = ncdir.variables['XLAT'][0,:,:] #latitude: XLAT(time, y, x)
    Lon_x = ncdir.variables['XLONG'][0,:,:] #longitude: XLONG(time, y, x)

    # Define grid for analysis
    min_lon = np.amin(Lon_x) #- 0.5
    max_lon = np.amax(Lon_x) #+ 0.5
    min_lat = np.amin(Lat_x) #- 0.5
    max_lat = np.amax(Lat_x) #+ 0.5

    # Approximation of ROI in grid_resolution 
    mask_size = np.floor( roi/ (model_reso) )

    # Flatten the grid for easier indexing
    grid_points = np.c_[Lon_x.ravel(), Lat_x.ravel()]
    #print( 'grid_points:',np.shape( grid_points ) )

    # Use KDTree for efficient nearest-neighbor search
    tree = cKDTree(np.c_[lon_obs, lat_obs])

    # Find points within each grid cell
    distances, indices = tree.query(grid_points, k=1)

    # Map diff_obserature values to the grid
    diff_obs_grid = np.full(Lon_x.shape, np.nan)
    diff_obs_grid.ravel()[np.arange(len(indices))] = diff_obs[indices]

    # Create a 2D Gaussian kernel
    sigma = mask_size / 10  # Standard deviation (adjustable)
    ax = np.arange(-(mask_size // 2), mask_size // 2 + 1)
    x, y = np.meshgrid(ax, ax)
    gaussian_kernel = np.exp(-(x**2 + y**2) / (2 * sigma**2))
    gaussian_kernel /= gaussian_kernel.sum()  # Normalize to sum to 1

    # Initialize output grid for mean diff_obseratures
    mean_diff_obs = np.full_like(diff_obs_grid, np.nan)
    mean_lon = np.full_like(diff_obs_grid, np.nan)
    mean_lat = np.full_like(diff_obs_grid, np.nan)

    # Iterate over the grid and apply a sliding mask
    for i in range(int(mask_size // 2), int(diff_obs_grid.shape[0] - mask_size // 2)):
        for j in range(int(mask_size // 2), int(diff_obs_grid.shape[1] - mask_size // 2)):

            # Extract the mask region
            mask_region = diff_obs_grid[
                i - int(mask_size // 2) : int(i + mask_size // 2 + 1),
                j - int(mask_size // 2) : int(j + mask_size // 2 + 1),
            ]

            # Apply Gaussian weights, ignoring NaN values
            valid_mask = ~np.isnan(mask_region)
            if valid_mask.any():
                #print('mask_region',np.shape(mask_region[valid_mask]))
                #print('gaussian kernel',np.shape(gaussian_kernel[valid_mask]))
                mean_diff_obs[i, j] = np.sum(
                    mask_region[valid_mask] * gaussian_kernel[valid_mask]
                ) / np.sum(gaussian_kernel[valid_mask])

            ## Compute the mean, ignoring NaN values
            #mean_diff_obs[i, j] = np.nanmean(mask_region)

            # Assign the location (center of the mask)
            mean_lon[i, j] = Lon_x[i, j]
            mean_lat[i, j] = Lat_x[i, j]

    d_smooth = {'smt_lon':mean_lon,'smt_lat':mean_lat,'smt_diff':mean_diff_obs}
    return d_smooth


def plot_test( d_smooth ):

    # Read WRF domain
    wrf_file = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
    d_wrf_d03 = ROIR.read_wrf_domain( wrf_file )

    f, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6.5,6), dpi=300)

    # Define the domain
    lat_min = d_wrf_d03['lat_min']
    lat_max = d_wrf_d03['lat_max']
    lon_min = d_wrf_d03['lon_min']
    lon_max = d_wrf_d03['lon_max']

    #Define Tb threshold
    min_T = -15
    max_T = 15
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)
    # plot differnet ROI
    for iroi in ROI:
        roi = ROI[iroi]['nonhydro']
        if int(roi) == 200:
            cs = ax.scatter(d_smooth[roi]['smt_lon'],d_smooth[roi]['smt_lat'],5,marker='s',c=d_smooth[roi]['smt_diff'],edgecolors='none', cmap='bwr', vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())
        elif int(roi) == 30:
            cs = ax.scatter(d_smooth[roi]['smt_lon'],d_smooth[roi]['smt_lat'],5,c=d_smooth[roi]['smt_diff'],edgecolors='none', cmap='bwr', vmin=min_T, vmax=max_T,transform=ccrs.PlateCarree())

    # Colorbar
    cbaxes = f.add_axes([0.9, 0.1, 0.02, 0.8])
    cbar = f.colorbar(cs, cax=cbaxes,fraction=0.046, pad=0.04, )
    cbar.ax.tick_params(labelsize=15)
    #caxes = f.add_axes([0.2, 0.1, 0.6, 0.02])
    #cbar = f.colorbar(cs, orientation="horizontal", cax=caxes)
    #cbar.ax.tick_params(labelsize=6)

    #subplot title
    matplotlib.rcParams['mathtext.fontset'] = 'custom'
    matplotlib.rcParams['mathtext.bf'] = 'STIXGeneral:italic:bold'
    font = {'size':9,}
    #ax.set_title('Yo', font, fontweight='bold')

    #title for all
    #f.suptitle(title_name, fontsize=6, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(1):
        gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.7, linestyle='--')

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
        gl.xlabel_style = {'size': 15}
        gl.ylabel_style = {'size': 15}

    des_name = 'test.png'
    #des_name = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Tb/IRch8_Obspace/'+DAtime+'_'+sensor+'_Obspace_limit.png'
    plt.savefig(des_name,dpi=300)
    print('Saving the figure: ',des_name)







if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR'
    MP = 'THO'

    model_reso = 3
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs','Hroi']

    start_time_str = '201709030600'
    end_time_str = '201709030600'
    Consecutive_times = True

    # limitations
    #ROI = {'small':{'hydro':'30','nonhydro':'30'}} 
    ROI = {'large':{'hydro':'0','nonhydro':'200'},'small':{'hydro':'30','nonhydro':'30'}} # KM: for hydrometer variables, for non-hydrometeor variables
    plot_scatter = True
    # ------------------------------------------------------   


    # Create experiment names
    Exper_name = 'CONV+IR-THO_onlyTCvitals' #UD.generate_one_name( Storm,DA,MP )
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


    # Average the Tbs in the area impacted by ROI
    for DAtime in IR_times:

        print('DAtime: '+ DAtime)
        # Read assimilated IR Tbs and their ROI
        file_Diag = big_dir+Storm+'/'+Exper_obs+'/run/'+DAtime+'/enkf/d03/fort.10000'
        d_so = Diag.Find_IR( file_Diag, fort_v )
        lon_so = d_so['lon']
        lat_so = d_so['lat']

        # Assemble the so data with a specific ROI
        #so_loc = {}
        #for iroi in ROI.keys(): # fort.10000 only records hydrometeor variables
        #    idx_roi = d_so['Hroi'] == int( ROI[iroi]['hydro'] )/model_reso 
            #so_loc[ROI[iroi]['hydro']] = [(a,b) for a, b in zip(lon_so[idx_roi], lat_so[idx_roi])] 

        # Read obs, Hxb, and Hxa of Tbs
        Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
        Tb_file = Hx_dir + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
         
        d_obsRes = ROIR.read_Tb_obsRes(Tb_file, sensor )
        lon_obs = d_obsRes['lon_obs']
        lat_obs = d_obsRes['lat_obs']
        #obs_minus_hxb = d_obsRes['Yo_obs'] - d_obsRes['meanYb_obs']
        obs_minus_hxa = d_obsRes['Yo_obs'] - d_obsRes['meanYa_obs']
        #obs_list = [(a, b) for a, b in zip(lon_obs, lat_obs)]

        d_smt_diff_hxb = {}
        # Apply convolution filters
        for iroi in ROI:
            # Find the index of each element in obs_list appearing in so_list
            #idx_exist = [
            #    next((i for i, obs in enumerate(obs_list) if np.array_equal(obs, element)), np.nan)
            #    for element in so_loc[ROI[iroi]['hydro']]
            #]
            #idx_exist = np.array( idx_exist )

            #!!!!!!!!!!! Assume no random order is applied !!!!!!!!!!!!!!!!1
            idx_roi = d_so['Hroi'] == int( ROI[iroi]['hydro'] )/model_reso
            idx_exist = idx_roi
            d_smt_diff_hxb[ROI[iroi]['nonhydro']] = convolution_IR_obsRes( lon_obs[idx_exist],lat_obs[idx_exist],obs_minus_hxa[idx_exist],int(ROI[iroi]['nonhydro']) )

        # Plot
        plot_test( d_smt_diff_hxb )










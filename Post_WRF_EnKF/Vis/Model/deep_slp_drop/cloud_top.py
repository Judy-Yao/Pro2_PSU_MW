#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

from numba import njit, prange
import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
from wrf import getvar, ll_to_xy
import math
import scipy as sp
import scipy.ndimage
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
from scipy import interpolate
import statistics

from Track_xbxa import read_HPI_model
import Util_data as UD
import Util_Vis
import Obspace_IR as ROIR


# ------------------------------------------------------------------------------------------------------
#           Operation: Perform Operations
# ------------------------------------------------------------------------------------------------------
@njit(parallel=True)
def cal_ctp(get_ctp,ver_coor,qhydro,layer):

    ctp = get_ctp
    ctp[:] = np.nan
    ip_valid = np.zeros( xmax*ymax )
    ip_valid[:] = np.nan # collect idx of point that has non-Nan value
    ip_seq = np.arange( xmax*ymax )

    for ip in prange(ctp.shape[0]):
        for il in range( ver_coor.shape[0]-1,0,-1 ):
            if (il-layer-1) < 0:
                break
            if qhydro[il,ip] <= 1e-6 and (qhydro[il-layer-1:il,ip] > 1e-6).all(): # see deep cloud qhydro[il-21:il,ip] > 1e-6
                slope = (ver_coor[il, ip]-ver_coor[il-1, ip])/(qhydro[il, ip]-qhydro[il-1, ip])
                ctp[ip] = ver_coor[il-1, ip]+(1e-6 - qhydro[il-1, ip])*slope
                ip_valid[ip] = ip_seq[ip]
                break
    return ctp,ip_valid

def total_qh( wrf_file, hydros):

    # Vertical dimension
    nLevel = 42
    # Domain dimension
    xmax = 297
    ymax = 297

    qhydro = np.zeros( [nLevel, xmax*ymax] )
    for var_name in hydros:
        ncdir = nc.Dataset( wrf_file, 'r')
        var = ncdir.variables[var_name][0,:,:,:]
        var = var.reshape( var.shape[0],-1 )
        qhydro = qhydro + var
    
    # Filter out zero elements
    #layer_all = qhydro[40,:]
    #non_zeros = [x for x in layer_all if x != 0]
    
    return qhydro

def plot_obs_ctp( wrf_files, d_ctp, d_obs ):

    # Read WRF domain
    d_wrf_d03 = UD.read_wrf_domain( wrf_files[1] )

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
    obs_s = axs.flat[0].scatter(d_obs['lon_obs'],d_obs['lat_obs'],1.5,c=d_obs['Yo_obs'],edgecolors='none', cmap=IRcmap, vmin=min_obs, vmax=max_obs,transform=ccrs.PlateCarree())
    #if any( hh in DAtime[8:10] for hh in ['00','06','12','18'] ):
    #    axs.flat[0].scatter(tc_lon, tc_lat, s=1, marker='*', edgecolors='darkviolet', transform=ccrs.PlateCarree())
    #axs.flat[0].add_patch(patches.Polygon(path,facecolor='none',edgecolor='white',linewidth=0.5 ))
    # Colorbar
    caxes = fig.add_axes([0.12, 0.1, 0.25, 0.02])
    obs_bar = fig.colorbar(obs_s,ax=axs[0],orientation="horizontal", cax=caxes)
    obs_bar.ax.tick_params(labelsize=6)

    # Model Cloud top

    # vertical coordinate
    if use_pressure:
        min_vc = 50
        max_vc = 850
        bounds = [100,200,300,400,500,600,700,800,900]
    else:
        min_vc = 0
        max_vc = 20
        bounds = np.linspace(min_vc,max_vc,9) 

    # plot
    if plot_scatter:
        if use_pressure:
            axs[1].scatter(d_ctp['lon'], d_ctp['lat'], s=1.5, c=d_ctp['input_ctp'], edgecolors='none', cmap='magma_r', vmin=min_vc, vmax=max_vc, transform=ccrs.PlateCarree())
            cs = axs[2].scatter(d_ctp['lon'], d_ctp['lat'], s=1.5, c=d_ctp['output_ctp'], edgecolors='none', cmap='magma_r', vmin=min_vc, vmax=max_vc, transform=ccrs.PlateCarree())
            caxes = fig.add_axes([0.4, 0.1, 0.5, 0.02])
            cbar = fig.colorbar(cs,ax=axs[1:],ticks=range(min_vc,max_vc+50,100),orientation="horizontal", cax=caxes)
            cbar.ax.tick_params(labelsize=6)
            cbar.ax.invert_xaxis()
        else:
            axs[1].scatter(d_ctp['lon'], d_ctp['lat'], s=1, c=d_ctp['input_ctp'], edgecolors='none', cmap='magma_r', vmin=min_vc, vmax=max_vc, transform=ccrs.PlateCarree())
            cs = axs[2].scatter(d_ctp['lon'], d_ctp['lat'], s=1, c=d_ctp['output_ctp'], edgecolors='none', cmap='magma_r', vmin=min_vc, vmax=max_vc, transform=ccrs.PlateCarree())
            caxes = fig.add_axes([0.4, 0.1, 0.5, 0.02])
            cbar = fig.colorbar(cs,ax=axs[1:],ticks=bounds,orientation="horizontal", cax=caxes)
            cbar.ax.tick_params(labelsize=6)
    else:
        cs = axs[1].contourf(d_ctp['lon'].reshape(xmax,ymax), d_ctp['lat'].reshape(xmax,ymax), d_ctp['input_ctp'].reshape(xmax,ymax),cmap='magma_r',levels=bounds,extend='both',transform=ccrs.PlateCarree())
        cs = axs[2].contourf(d_ctp['lon'].reshape(xmax,ymax), d_ctp['lat'].reshape(xmax,ymax), d_ctp['output_ctp'].reshape(xmax,ymax),cmap='magma_r',levels=bounds,extend='both',transform=ccrs.PlateCarree())
        caxes = fig.add_axes([0.4, 0.1, 0.5, 0.02])
        color_ticks = bounds
        cbar = fig.colorbar(cs,ax=axs[1:],ticks=color_ticks,orientation="horizontal", cax=caxes)
        bounds_str =  [ str(item) for item in color_ticks ]
        cbar.ax.set_xticklabels( bounds_str)
        cbar.ax.tick_params(labelsize=6)
        #cbar.ax.invert_xaxis()

    #Subplot title
    font = {'size':6,}
    axs.flat[0].set_title('Ch'+ch_list[0]+':Yo', font, fontweight='bold')
    ## format the statistics
    axs.flat[1].set_title('CTP(Xb)', font, fontweight='bold')
    axs.flat[2].set_title('CTP(Xa)', font, fontweight='bold') 

    # Set title
    if DeepCloud and ShallowCloud:
        title_name = Storm+': '+Exper_name+' @ '+DAtime+'  \nCTP threshold: 1e-6 kg/kg'
    elif DeepCloud and not ShallowCloud:
        title_name = Storm+': '+Exper_name+' @ '+DAtime+'\n Deep Clouds\n(condition: at the boundary,20 consecutive bottom layers with Qhydro > 1e-6)'
    elif not DeepCloud and ShallowCloud:
        title_name = Storm+': '+Exper_name+' @ '+DAtime+'\n Shallow Clouds\n(condition: at the boundary,NOT 20 consecutive bottom layers with Qhydro > 1e-6)'
    else:
        raise ValueError('Such condition does not exist!')
    fig.suptitle(title_name, fontsize=6, fontweight='bold')


    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(3):
        gl = axs.flat[j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5,alpha=0.7,color='gray',linestyle='--')

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

    if plot_scatter and use_pressure:
        figure_des=plot_dir+DAtime+'_'+MP+'_cloud_top_scatter_pres.png'
    elif plot_scatter and not use_pressure:
        figure_des=plot_dir+DAtime+'_'+MP+'_cloud_top_scatter_height.png'
    elif not plot_scatter and use_pressure:
        figure_des=plot_dir+DAtime+'_'+MP+'_cloud_top_contour_pres.png'
    else:
        figure_des=plot_dir+DAtime+'_'+MP+'_cloud_top_contour_height.png'
    plt.savefig(figure_des, dpi=400)
    print('Saving the figure: ', figure_des)


def find_ct( wrf_file ):

    # Domain dimension
    xmax = 297
    ymax = 297

    # Read vertical coordinate
    ncdir = nc.Dataset( wrf_file, 'r')

    if use_pressure:
        p = ncdir.variables['P'][0,:,:,:] # perturbation
        pb = ncdir.variables['PB'][0,:,:,:]
        tmp_p = (p + pb)/100
        pres = tmp_p.reshape( tmp_p.shape[0],-1 )
        ver_coor = pres
    else:
        PHB = ncdir.variables['PHB'][0,:,:,:]
        PH = ncdir.variables['PH'][0,:,:,:]
        geoHkm = (PHB+PH)/9.8/1000 # in km
        geoHkm = geoHkm.reshape( geoHkm.shape[0],-1)
        geoHkm_half_eta = (geoHkm[:-1]+geoHkm[1:])/2
        geoHkm_half_eta = np.ma.getdata(geoHkm_half_eta)
        ver_coor = geoHkm_half_eta


    # Find points of interest: hydrometeor value as threshold
    qhydro = total_qh( wrf_file,hydros )
    #idx_p = range(xmax*ymax)

    # Calculate CTP
    get_ctp = np.zeros( xmax*ymax )
    if DeepCloud and ShallowCloud:
        ctp,idx_p = cal_ctp(get_ctp,ver_coor,qhydro,0)
        idx_p = idx_p.astype(int)
    elif DeepCloud and not ShallowCloud:
        ctp,idx_p = cal_ctp(get_ctp,ver_coor,qhydro,layer_depth)
        idx_p = idx_p.astype(int)
    elif not DeepCloud and ShallowCloud:
        # deep cloud
        ctp_deep,ip_deep = cal_ctp(get_ctp,ver_coor,qhydro,layer_depth)
        ip_deep = ip_deep.astype(int)
        # all sky
        ctp_all,ip_all = cal_ctp(get_ctp,ver_coor,qhydro,0)
        ip_all = ip_all.astype(int)
        # for shallow cloud
        mask_all = ~np.isnan(ip_all)
        ip_all = ip_all[mask_all]
        mask_deep = ~np.isnan(ip_deep)
        ip_deep = ip_deep[mask_deep]
        ip_shallow = list(set(ip_all).difference(set(ip_deep)))
        ctp_sw = np.zeros( xmax*ymax )
        ctp_sw[:] = np.nan
        for i in ip_shallow:
            ctp_sw[i] = ctp_all[i]
        ctp = ctp_sw
        idx_p = ip_shallow
    else:
        raise ValueError('Such condition does not exist!')
    return ctp,idx_p


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    model_resolution = 3000 #m
    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR'
    MP = 'THO'
    hydros =  ['QCLOUD','QRAIN','QICE','QSNOW','QGRAUP']

    # obs
    sensor = 'abi_gr'
    ch_list = ['8',]

    # Time range set up
    start_time_str = '201709030000'
    end_time_str = '201709030000'
    Consecutive_times = True

    # If study: deep_slp_incre
    deep_slp_incre = False
    incre_slp_th = 0 # threshold of increment, unit:hpa
    radius_th = 200 #km

    # how to define cloud top
    DeepCloud = True
    ShallowCloud = True
    layer_depth = 20
    use_pressure = False
    plot_scatter = False 
 
    # Dimension  
    xmax = 297
    ymax = 297 
    # --------------------------------------------------

    # Create experiment names
    Exper_name =  UD.generate_one_name( Storm,DA,MP )

    if not Consecutive_times:
        DAtimes = ['201709180000',]
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]

    # Identify the min slp pressure in both background and analysis
    if deep_slp_incre:
        # ----- Read min slp from model-----------------
        HPI_models = {}
        DAtimes_dir = [big_dir+Storm+'/'+Exper_name+'/fc/'+it for it in DAtimes]
        file_kinds = ['wrf_enkf_input_d03_mean','wrf_enkf_output_d03_mean']
        for ifk in file_kinds:
            idx = file_kinds.index( ifk )
            HPI_models[ifk] = read_HPI_model( Storm, Exper_name, ifk, DAtimes_dir )
        incre_slp = np.array(HPI_models['wrf_enkf_output_d03_mean']['min_slp']) - np.array(HPI_models['wrf_enkf_input_d03_mean']['min_slp'])

    # Cloud top
    for DAtime in DAtimes:
        start_time=time.process_time()
        if deep_slp_incre:
    
            pass

        else:
            idx_t = DAtimes.index( DAtime )
            print('At '+DAtime)
            # Read IR obs
            Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
            Tb_file = Hx_dir + "/mean_obs_res_d03_" + DAtime + '_' +  sensor + '.txt'
            d_obs = ROIR.read_Tb_obsRes(Tb_file, sensor )

            # Read model-related: cloud top pressure for xb and xa
            d_ctp = {}
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            wrf_files = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']    
            # Calculate cloud top and attach the correpsonding location 
            ctp,idx_p = find_ct( wrf_files[0] )
            d_ctp['input_ctp'] = ctp
            d_ctp['input_idx'] = idx_p
            ctp,idx_p = find_ct( wrf_files[1] )
            d_ctp['output_ctp'] = ctp
            d_ctp['output_idx'] = idx_p
            #print(np.amax( d_ctp['input_ctp'] ))
            #print(np.amax( d_ctp['output_ctp'] ))
            # Read WRF domain
            ncdir = nc.Dataset( wrf_files[0], 'r')
            xlat = ncdir.variables['XLAT'][0,:,:].flatten()
            xlon = ncdir.variables['XLONG'][0,:,:].flatten()
            d_ctp['lat'] = xlat
            d_ctp['lon'] = xlon

            # ------ Plot -------------------
            plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/Cloud_top/'
            plotdir_exists = os.path.exists( plot_dir )
            if plotdir_exists == False:
                os.mkdir(plot_dir) 
            plot_obs_ctp( wrf_files, d_ctp, d_obs )

        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')





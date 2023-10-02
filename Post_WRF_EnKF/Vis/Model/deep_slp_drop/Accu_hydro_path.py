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

from Track_xbxa import read_HPI_model
import Util_data as UD
import Util_Vis

# ------------------------------------------------------------------------------------------------------
#           Operation: Perform Operations
# ------------------------------------------------------------------------------------------------------
# Calculate hydro mass at any grid point
@njit(parallel=True)
def hydro_mass( full_p, tv, geoHm, q ):
    R = 287.06
    res = np.zeros( (q.shape[0],q.shape[1]), )
    res[:] = np.nan
    for im in prange( q.shape[1] ):
        for il in range( q.shape[0] ):
            zdiff = geoHm[il+1,im] - geoHm[il,im]
            res[il,im] = (full_p[il,im]/(R * tv[il,im])) * q[il,im] * zdiff #* model_resolution**2
    # make sure all values are reasonable
    assert res.any() != np.nan
    return res

# Calculate the vertical coordinate where threshold value is (The furthest depth of cloud that IR can see)
@njit(parallel=True)
def cal_cdir(get_cdir,ver_coor,accu_qhydro,accu_th):

    cd_ir = get_cdir
    cd_ir[:] = np.nan
    ip_valid = np.zeros( xmax*ymax )
    ip_valid[:] = np.nan # collect idx of point that has non-Nan value
    ip_seq = np.arange( xmax*ymax )

    for ip in prange(cd_ir.shape[0]):
        for il in range( ver_coor.shape[0]-1,0,-1 ):
            if (il-1) < 0:
                break
            if accu_qhydro[il,ip] <= accu_th and accu_qhydro[il-1,ip] > accu_th:
                slope = (ver_coor[il, ip]-ver_coor[il-1, ip])/(accu_qhydro[il, ip]-accu_qhydro[il-1, ip])
                cd_ir[ip] = ver_coor[il-1, ip]+(accu_th - accu_qhydro[il-1, ip])*slope
                ip_valid[ip] = ip_seq[ip]
                break
    return cd_ir,ip_valid


# Read variables at model resolution/location
def read_Tb_modelRes(Tb_file, sensor ):

    lat_m = []
    lon_m = []
    Yo_obs = []
    meanYb_m = []
    meanYa_m = []

    # Read records
    print('Reading ', Tb_file)
    with open(Tb_file) as f:
        next(f)
        all_lines = f.readlines()

    for line in all_lines:
        split_line = line.split()
        lat_m.append( float(split_line[0]) )
        lon_m.append( float(split_line[1]) )
        meanYb_m.append( float(split_line[3]) )
        meanYa_m.append( float(split_line[4]) )

    lat_m = np.array( lat_m )
    lon_m = np.array( lon_m )
    meanYb_m = np.array( meanYb_m )
    meanYa_m = np.array( meanYa_m )

    dict_Tb_all = {'lat_m':lat_m, 'lon_m':lon_m, 'meanYb_m':meanYb_m, 'meanYa_m':meanYa_m}
    return dict_Tb_all





# ------------------------------------------------------------------------------------------------------
#           Operation: Calculate accumulated hydro mass
# ------------------------------------------------------------------------------------------------------
def compute_accu_hydromass( wrf_files, hydros, idx_t, HPI_models=None):

    d_hydro = {}

    # Find points of interest
    if deep_slp_incre:
        # Find the min slp in enkf input as the anchor point
        # Find the circled area with the min slp as the center
        output_file = wrf_files[1]
        ncdir = nc.Dataset( output_file, 'r')
        anchor_ij = ll_to_xy(ncdir, HPI_models['wrf_enkf_output_d03_mean']['lat'][idx_t], HPI_models['wrf_enkf_output_d03_mean']['lon'][idx_t])
        idx_area = UD.find_circle_area_model_ij( output_file, anchor_ij.values[0], anchor_ij.values[1], radius_th, model_resolution/1000) # in km
    else:
        idx_area = range(xmax*ymax)
        
    # Calculate fields needed for mass of hydrometeors
    full_p = np.zeros( [len(wrf_files), nLevel, len(idx_area)] )
    tv_k = np.zeros( [len(wrf_files), nLevel, len(idx_area)] )
    geoHm = np.zeros( [len(wrf_files), nLevel+1, len(idx_area)] )
    for wrf_file in wrf_files:
        ifile = wrf_files.index(wrf_file)
        ncdir = nc.Dataset( wrf_file, 'r')
        # full pressure
        p = ncdir.variables['P'][0,:,:,:] # perturbation
        pb = ncdir.variables['PB'][0,:,:,:]
        tmp_p = p + pb
        tmp_p = tmp_p.reshape( tmp_p.shape[0],-1 )
        full_p[ifile,:,:] = tmp_p[:,idx_area]
        # geopotential height
        ph = ncdir.variables['PH'][0,:,:,:] # perturbation
        phb = ncdir.variables['PHB'][0,:,:,:]
        tmp_geoHm = (ph+phb)/9.81 # in meter
        tmp_geoHm = tmp_geoHm.reshape( tmp_geoHm.shape[0],-1 )
        geoHm[ifile,:,:] = tmp_geoHm[:,idx_area]
        # virtual temperature
        tv = getvar(ncdir,'tv',units='K')
        tmp_tv = tv.values
        tmp_tv = tmp_tv.reshape( tmp_tv.shape[0],-1 )
        tv_k[ifile,:,:] = tmp_tv[:,idx_area]
    # Accumulate the mass of hydrometeor 
    for var_name in hydros:
        accu_hydro = np.zeros( [len(wrf_files), nLevel, len(idx_area)] )
        accu_hydro[:] = np.nan
        for wrf_file in wrf_files:
            ifile = wrf_files.index(wrf_file)
            ncdir = nc.Dataset( wrf_file, 'r')
            ivar = hydros.index(var_name)
            var = ncdir.variables[var_name][0,:,:,:]
            var = var.reshape( var.shape[0],-1 )
            # calculate the mass of hydrometeor at any given point for a file
            tmp = hydro_mass( full_p[ifile,:,:], tv_k[ifile,:,:], geoHm[ifile,:,:], var[:,idx_area] )
            # accumulate
            if path_from_top:
                accu_tmp = np.zeros( [len(idx_area)] )
                for il in range(nLevel-1,-1,-1): # 41,...,0
                    accu_tmp[:] = accu_tmp[:] + tmp[il,:]
                    accu_hydro[ifile,il,:] = accu_tmp[:] 
            else:
                pass
        d_hydro[var_name] = accu_hydro

    # Calculate the vertical coordinate
    if ver_use_press:
        return d_hydro, full_p
    else:
        return d_hydro, geoHm

# ------------------------------------------------------------------------------------------------------
#           Operation: Plot hydro mass
# ------------------------------------------------------------------------------------------------------
def Plot_DM_all_oneTime( DAtime,v_interest,d_hydro,ver_coor ):

    if ver_use_press:
        ver_coor = np.mean( ver_coor, axis=0 )
        ver_coor = np.mean( ver_coor, axis=1)/100 # convert to hPa
    else:
        ver_coor = np.mean( ver_coor, axis=0)
        ver_coor = np.mean( ver_coor, axis=1 )/1000 # convert to km
        ver_coor = (ver_coor[:-1]+ver_coor[1:])/2

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)

    # Manually set discrete values on x and y axis and interpolate data to these values
    ## x axis: range of accumulated water path value
    if Storm == 'IRMA' and MP == 'THO':
        x_range = np.arange(0,500.5,50)
    else:
        x_range = np.arange(0,290.5,50) #2100.5
    x_axis_rg = range(len(x_range))
    f_xinterp = interpolate.interp1d( x_range, x_axis_rg)
    ## y axis: model vertical coordinate
    if ver_use_press:
        pass
    else:
        y_range = np.arange(0,31,1)
        y_axis_rg = range(len(y_range))
        f_yinterp = interpolate.interp1d( y_range, y_axis_rg)
        loc_iny = f_yinterp( ver_coor )

    labels = ['Xb_','Xa_']
    lstyle = ['--','-']
    Color = ['#f032e6','#911eb4','#4363d8','#f58231','#469990']

    for ifile in range( len(wrf_files) ):
        #PF_all = np.zeros( (nLevel,) )
        for var in hydros:
            idx = hydros.index( var )
            PF_x = np.mean(d_hydro[var][ifile,:,:],axis=1)*1000 # convert from kg/m2 to g/m2
            loc_inx = f_xinterp( PF_x )
            ax.plot( loc_inx,loc_iny,Color[idx],linewidth=3,label=labels[ifile]+var,linestyle=lstyle[ifile] )
            #PF_all = PF_all + PF_x
        #loc_inx = f_xinterp( PF_all )
        #if ifile == 0:
        #    ax.plot( loc_inx,loc_iny,'black',linewidth=2,label=labels[ifile]+'All',linestyle=lstyle[ifile] )
        #else:
        #    ax.plot( loc_inx,loc_iny,'black',linewidth=2,label=labels[ifile]+'All',linestyle=lstyle[ifile] )

    # Plot a line indicating the threshold value
    th_loc = f_xinterp(accu_th)
    ax.axvline(x=th_loc,color='black',linestyle='-',linewidth=2)

    # set lables
    ax.legend(loc='upper right',fontsize='10')
    # set X label
    if Storm == 'IRMA' and MP == 'THO':
        xticks_loc = list(x_axis_rg[::2])
        xlabel_like = list(x_range[::2])
        ax.set_xlim(xmin=0,xmax=f_xinterp(500.5))
    else:
        xticks_loc = list(x_axis_rg[::1])
        xlabel_like = list(x_range[::1])
        ax.set_xlim(xmin=0,xmax=f_xinterp(250))
    xticks_loc.insert(0,th_loc)
    xlabel_like.insert(0,accu_th)
    ax.set_xticks( xticks_loc )
    ax.set_xticklabels( [str(it) for it in xlabel_like],fontsize=15 )
    ax.set_xlabel('Accumulated Mass (gram m-2)',fontsize=15)
    ax.set_xlim(xmin=0)
    # set Y label
    if ver_use_press:
        pass
    else:
        ylabel_like = [0.0,5.0,10.0,15.0,20.0]
        yticks = []
        list_y_range = list(y_range)
        for it in ylabel_like:
            yticks.append( list_y_range.index(it) )
        ax.set_yticks( yticks )
        ax.set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
        ax.set_ylabel('Height (km)',fontsize=15)
        ax.set_ylim(ymin=0,ymax=20.5) # cut off data above 25km

    # Set title
    title_name = 'Domain-mean Profile: Top-to-down Accumulated Hydro Mass'
    ax.set_title( title_name,fontweight="bold",fontsize='12' )
    fig.suptitle(Storm+': '+Exper_name+' '+DAtime , fontsize=10, fontweight='bold')

    # Save the figure
    save_des = plot_dir+DAtime+'_DM_VP_accumulated_each_hydro.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


    return None


def Plot_depth_IRsee( DAtime, wrf_files, d_hydro, ver_coor ):

    # Read WRF domain
    d_wrf_d03 = UD.read_wrf_domain( wrf_files[1] )
    ncdir = nc.Dataset(wrf_files[1], 'r')
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()

    # Read simulated Tb of Xb and Xa
    Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+DAtime+'/'
    Tb_file = Hx_dir + "/mean_model_res_d03_" + DAtime + '_' +  sensor + '.txt'
    d_Tb = read_Tb_modelRes(Tb_file, sensor )

    # Calculate the vertical coordinate where IR can see the furthest depth of cloud 
    tmp = d_hydro['QCLOUD']+d_hydro['QRAIN']+d_hydro['QICE']+d_hydro['QSNOW']+d_hydro['QGRAUP']    
    cdirs = np.zeros( (len(wrf_files),xmax*ymax) ) 
    for ifile in wrf_files:
        get_cdir = np.zeros( (xmax*ymax,) )
        idx = wrf_files.index( ifile )
        cdir,idx_p = cal_cdir(get_cdir,ver_coor[idx,:,:],tmp[idx,:,:]*1000,accu_th)
        if ver_use_press:
            cdirs[idx,:] = cdir/100 # convert to hPa
        else:
            cdirs[idx,:] = cdir/1000 # conver to km

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(8.5,8), dpi=400)

    # domain area
    lat_min = d_wrf_d03['lat_min']
    lon_min = d_wrf_d03['lon_min']
    lat_max = d_wrf_d03['lat_max']
    lon_max = d_wrf_d03['lon_max']

    # map projection
    for i in range(4):
        ax.flat[i].set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax.flat[i].coastlines (resolution='10m', color='black', linewidth=1)

    # Simulated Obs    
    min_obs = 185
    max_obs = 325
    IRcmap = Util_Vis.IRcmap( 0.5 )
    ax[0,0].scatter(d_Tb['lon_m'],d_Tb['lat_m'],1.5,c=d_Tb['meanYb_m'],edgecolors='none', cmap=IRcmap, vmin=min_obs, vmax=max_obs,transform=ccrs.PlateCarree())
    ir_f = ax[1,0].scatter(d_Tb['lon_m'],d_Tb['lat_m'],1.5,c=d_Tb['meanYa_m'],edgecolors='none', cmap=IRcmap, vmin=min_obs, vmax=max_obs,transform=ccrs.PlateCarree())
    # Colorbar
    cbaxes = fig.add_axes([0.00, 0.1, 0.03, 0.8])
    tb_range = np.linspace(min_obs,max_obs,8)
    IR_bar = fig.colorbar(ir_f,cax=cbaxes,ticks=tb_range,fraction=0.046, pad=0.04)
    #IR_bar.ax.set_ylabel('Brightness Temperature (K)')
    IR_bar.ax.tick_params(labelsize=12)
 
    # cloud depth that IR can see
    if ver_use_press:
        min_vc = 50
        max_vc = 850
        bounds = [100,200,300,400,500,600,700,800,900]
    else:
        min_vc = 0
        max_vc = 21
        bounds = np.linspace(min_vc,max_vc,8)
    if plot_scatter:
        ax[0,1].scatter(xlon,xlat,s=2,c=cdirs[0,:], edgecolors='none', cmap='magma_r', vmin=min_vc, vmax=max_vc, transform=ccrs.PlateCarree())
        cdir_f = ax[1,1].scatter(xlon,xlat,s=2,c=cdirs[1,:], edgecolors='none', cmap='magma_r', vmin=min_vc, vmax=max_vc, transform=ccrs.PlateCarree())
        cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
        cdir_bar = fig.colorbar(ticks=range(min_vc,max_vc+50,100),cax=cbaxes,fraction=0.046, pad=0.04)
        cdir_bar.ax.tick_params(labelsize=12)
    else:
        ax[0,1].contourf(xlon.reshape(xmax,ymax),xlat.reshape(xmax,ymax),cdirs[0,:].reshape(xmax,ymax),cmap='magma_r',levels=bounds,extend='both',transform=ccrs.PlateCarree())
        cdir_f = ax[1,1].contourf(xlon.reshape(xmax,ymax),xlat.reshape(xmax,ymax),cdirs[1,:].reshape(xmax,ymax),cmap='magma_r',levels=bounds,extend='both',transform=ccrs.PlateCarree()) 
        cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
        color_ticks = bounds
        cdir_bar = fig.colorbar(cdir_f,cax=cbaxes,fraction=0.046, pad=0.04)
        bounds_str =  [ str(item) for item in color_ticks ]
        cdir_bar.ax.set_yticklabels( bounds_str)
        cdir_bar.ax.tick_params(labelsize=12)


    # Set title
    title_name = Storm+': '+Exper_name+' '+DAtime+' \n Cloud Depth IR can See, TH='+str(accu_th)+' gram m-2'
    fig.suptitle(title_name, fontsize=12, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(2):
        for i in range(2):
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

            if j==0:
                gl.ylabels_left = True
                gl.ylabels_right = False
            else:
                gl.ylabels_left = False
                gl.ylabels_right = False

            if i == 1:
                gl.xlabels_top = False
                gl.xlabels_bottom = True
            else:
                gl.xlabels_top = False
                gl.xlabels_bottom = False

            gl.ylocator = mticker.FixedLocator(lat_ticks)
            gl.xlocator = mticker.FixedLocator(lon_ticks)
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 10}
            gl.ylabel_style = {'size': 10}


    # Save the figure
    if plot_scatter:
        save_des = plot_dir+DAtime+'_cloud_depth_IRcanSee_scatter.png'
    else:
        save_des = plot_dir+DAtime+'_cloud_depth_IRcanSee_contourf.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


    return None



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    model_resolution = 3000 #m

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR'
    MP = 'WSM6'
    hydros =  ['QCLOUD','QRAIN','QICE','QSNOW','QGRAUP']
    
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709040700'
    end_time_str = '201709040700'
    Consecutive_times = True

    deep_slp_incre = False
    incre_slp_th = 0 # threshold of increment, unit:hpa
    plot_circle = True
    radius_th = 200 # km

    path_from_top = True # accumulate water from top to down
    accu_th = 30.0 
    ver_use_press = False
    each_water = True
    domain_mean = True
    depth_IR_see = False
    plot_scatter = False
    # ------------------------------------------------------  
    # Dimension of the domain
    nLevel = 42
    xmax = 297
    ymax = 297

    # Create experiment names
    Exper_name =  UD.generate_one_name( Storm,DA,MP )

    # Times
    if not Consecutive_times:
        DAtimes = ['201709180000',]
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 1))]
        DAtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]


    if not deep_slp_incre:
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/Accu_HydroPath/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

        start_time=time.process_time()
        for DAtime in DAtimes:
            idx_t = DAtimes.index( DAtime )
            print('At '+DAtime)
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+DAtime+'/'
            wrf_files = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
            d_hydro,ver_coor = compute_accu_hydromass( wrf_files, hydros, idx_t) # hydro mass at any grid point

            # --- Condition
            if each_water:
                v_interest = hydros
            else:
                v_interest = ['liquid','ice','all_hydro']
            
            if domain_mean:
                Plot_DM_all_oneTime( DAtime,v_interest,d_hydro,ver_coor ) 

            if depth_IR_see:
                Plot_depth_IRsee( DAtime, wrf_files, d_hydro, ver_coor )

        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')




 

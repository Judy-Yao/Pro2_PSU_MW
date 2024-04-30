
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
from matplotlib import ticker
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
from scipy import interpolate
from fast_histogram import histogram2d as hist2d

#from Track_xbxa import read_HPI_model
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


# ------------------------------------------------------------------------
#                    Object: Tbs of Forecasts
# ------------------------------------------------------------------------

# Read crtm calcuated IR data from one binary file
def read_simu_IR_one(Hxb_file, ch_list):

    xmax = 297
    ymax = 297
    print('Hxb_file:' + Hxb_file)

    Hxb_data = np.fromfile(Hxb_file,dtype='<f4') # <: little endian; f: float; 4: 4 bytes
    n_ch = len(Hxb_data)/(xmax*ymax) - 2
    n_ch = int(n_ch)
    if n_ch != len(ch_list):
        print('Error!! # of channels in data is '+str(n_ch))
    Hxb_sim = Hxb_data[:].reshape(n_ch+2,ymax,xmax)

    dict_simu_Tb = {}
    dict_simu_Tb['Lon_x'] = Hxb_sim[0,:,:]
    dict_simu_Tb['Lat_x'] = Hxb_sim[1,:,:]
    dict_simu_Tb['Ch_x'] = ch_list[0]
    dict_simu_Tb['Tb_x'] = Hxb_sim[2,:,:]

    return dict_simu_Tb


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
    geoHm = np.zeros( [len(wrf_files), nLevel, len(idx_area)] )
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
        geoHkm_half_eta = (tmp_geoHm[:-1,:]+tmp_geoHm[1:,:])/2
        geoHm[ifile,:,:] = geoHkm_half_eta[:,idx_area]
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
            tmp = hydro_mass( full_p[ifile,:,:], tv_k[ifile,:,:], geoHm[ifile,:,:], var.data[:,idx_area] )
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
def Plot_DM_all_oneTime( FCtime,v_interest,d_hydro,ver_coor ):

    if ver_use_press:
        ver_coor = np.mean( ver_coor, axis=0 )
        ver_coor = np.mean( ver_coor, axis=1)/100 # convert to hPa
    else:
        ver_coor = np.mean( ver_coor, axis=0)
        ver_coor = np.mean( ver_coor, axis=1 )/1000 # convert to km

    # Set up figure
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)

    # Manually set discrete values on x and y axis and interpolate data to these values
    ## x axis: range of accumulated water path value
    #if Storm == 'IRMA' and MP == 'THO':
    #    x_range = np.arange(0,500.5,50)
    #else:
    #    x_range = np.arange(0,500.5,50) #2100.5
   
    x_range = list(np.logspace(0,4,num=50)) # start: 10**0, end:10**3
    x_range.insert(0,0)
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
        # Read simulated Tb of Xb and Xa
        Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+FCtime+'/'
        Tb_file = Hx_dir + "/mean_model_res_d03_" + FCtime + '_' +  sensor + '.txt'
        d_Tb = read_Tb_modelRes(Tb_file, sensor )
        print(np.shape(d_Tb['meanYb_m']))
        if ifile == 0: 
            condi = d_Tb['meanYb_m'] < 210
        else:
            condi = d_Tb['meanYa_m'] < 210
        idx_area = np.where( condi )[0]
       
        PF_all = np.zeros( (nLevel,) )
        for var in hydros:
            idx = hydros.index( var )
            print(np.shape(d_hydro[var][ifile,:,idx_area]))
            PF_x = np.mean(d_hydro[var][ifile,:,idx_area],axis=0)*1000 # convert from kg/m2 to g/m2
            loc_inx = f_xinterp( PF_x )
            ax.plot( loc_inx,loc_iny,Color[idx],linewidth=3,label=labels[ifile]+var,linestyle=lstyle[ifile] )
            PF_all = PF_all + PF_x
        print(np.max(PF_all))
        loc_inx = f_xinterp( PF_all )
        if ifile == 0:
            ax.plot( loc_inx,loc_iny,'black',linewidth=2,label=labels[ifile]+'All',linestyle=lstyle[ifile] )
        else:
            ax.plot( loc_inx,loc_iny,'black',linewidth=2,label=labels[ifile]+'All',linestyle=lstyle[ifile] )

    # Plot a line indicating the threshold value
    th_loc = f_xinterp(accu_th)
    ax.axvline(x=th_loc,color='black',linestyle='-',linewidth=2)

    # set lables
    ax.legend(loc='upper right',fontsize='10')

    # set X label
    xlabel_like = [0,1,10,20,100,1000,3162]
    xticks_loc = []
    for it in xlabel_like:
        xticks_loc.append( f_xinterp(it) )
    ax.set_xlim(xmin=f_xinterp(0),xmax=f_xinterp(3162))
    #if Storm == 'IRMA' and MP == 'THO':
    #    xticks_loc = list(x_axis_rg[::2])
    #    xlabel_like = list(x_range[::2])
    #    ax.set_xlim(xmin=0,xmax=f_xinterp(500.5))
    #else:
    #    xticks_loc = list(x_axis_rg[::1])
    #    xlabel_like = list(x_range[::1])
    #    ax.set_xlim(xmin=0,xmax=f_xinterp(250))
    #xticks_loc.insert(0,th_loc)
    #xlabel_like.insert(0,accu_th)
    ax.set_xticks( xticks_loc )
    ax.set_xticklabels( ['0',r'$10^{0}$',r'$10^{1}$','20',r'$10^{2}$',r'$10^{3}$',r'$10^{3.5}$'],fontsize=15 )
    #ax.set_xticklabels( [str(it) for it in xlabel_like],fontsize=15 )
    ax.set_xlabel('Accumulated Mass (gram m-2)',fontsize=15)
    ax.set_xlim(xmin=0)
    # set Y label
    if ver_use_press:
        pass
    else:
        ylabel_like = [10.0,12.5,15.0,17.5,20.0] #[0.0,5.0,10.0,15.0,20.0]
        yticks = []
        for it in ylabel_like:
            yticks.append( f_yinterp( it ) )
        ax.set_yticks( yticks )
        ax.set_yticklabels( [str(it) for it in ylabel_like],fontsize=15 )
        ax.set_ylabel('Height (km)',fontsize=15)
        ax.set_ylim(ymin=10,ymax=20.5) # cut off data above 25km

    # Set title
    title_name = 'Domain-mean Profile: Top-to-down Accumulated Hydro Mass'
    ax.set_title( title_name,fontweight="bold",fontsize='12' )
    fig.suptitle(Storm+': '+Exper_name+' '+FCtime , fontsize=10, fontweight='bold')

    # Save the figure
    save_des = plot_dir+FCtime+'_DM_VP_accumulated_each_hydro.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


    return None


def Plot_depth_IRsee( FCtime, wrf_files, d_hydro, ver_coor ):

    # Read WRF domain
    d_wrf_d03 = UD.read_wrf_domain( wrf_files[1] )
    ncdir = nc.Dataset(wrf_files[1], 'r')
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()

    # Read simulated Tb of Xb and Xa
    Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+FCtime+'/'
    Tb_file = Hx_dir + "/mean_model_res_d03_" + FCtime + '_' +  sensor + '.txt'
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
    title_name = Storm+': '+Exper_name+' '+FCtime+' \n Cloud Depth IR can See, TH='+str(accu_th)+' gram m-2'
    fig.suptitle(title_name, fontsize=12, fontweight='bold')

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))

    for j in range(2):
        for i in range(2):
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

            if j==0:
                gl.left_labels = True
                gl.right_labels = False
            else:
                gl.left_labels = False
                gl.right_labels = False

            if i == 1:
                gl.top_labels = False
                gl.bottom_labels = True
            else:
                gl.top_labels = False
                gl.bottom_labels = False

            gl.ylocator = mticker.FixedLocator(lat_ticks)
            gl.xlocator = mticker.FixedLocator(lon_ticks)
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 10}
            gl.ylabel_style = {'size': 10}


    # Save the figure
    if plot_scatter:
        save_des = plot_dir+FCtime+'_cloud_depth_IRcanSee_scatter.png'
    else:
        save_des = plot_dir+FCtime+'_cloud_depth_IRcanSee_contourf.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()

    return None




##################################################################
# One Scene
##################################################################

def cloudHeight_Tb_2D( cdirs ):

    # temporary parameters
    min_ch = 7.5
    max_ch = 19
    min_tb = 185
    max_tb = 250
    number_bins = 100

    # Read simulated Tb of Xb and Xa
    #Hx_dir = big_dir+Storm+'/'+Exper_name+'/Obs_Hx/IR/'+FCtime+'/'
    Hx_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/JerryRun/IR_THO/wrf_df/201708241200/'
    Tb_file = Hx_dir+'/CRTM_'+ wrfname + '.bin'
    d_Tb = read_simu_IR_one(Tb_file, ch_list)

    # bin
    d_count = hist2d(cdirs[0,:],d_Tb['Tb_x'].flatten(),range=[[min_ch,max_ch],[min_tb,max_tb]],bins=number_bins)
    # Compute log10 of non-zero elements, return NAN where elements are ZERO
    with np.errstate(divide='ignore', invalid='ignore'):
        dcount = np.where(d_count != 0, np.log10(d_count), np.nan)

    # Set up figure
    fig,ax = plt.subplots(1, 1, figsize=(10,6), dpi=300 )   

    # verification
    #ax.scatter(cdirs[0,:], d_Tb['Tb_x'], s=2)

    # new map
    # Customize the colormap
    color_intervals = [0,0.5,1,1.5,2,2.5,3,3.5]
    exist_cmap = plt.cm.rainbow #.reversed()
    colors = exist_cmap(np.linspace(0,1,len(color_intervals)-1))
    new_map = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals)-1)

    # Plot 2d histogram
    show = ax.imshow(np.transpose(dcount),cmap=new_map,aspect=0.5,vmin=0,vmax=3.5) 
    # Add color bar below the plot
    caxes = fig.add_axes([0.12, 0.88, 0.78, 0.03])
    color_bar = fig.colorbar(show,cax=caxes,orientation='horizontal')#ticks=bounds)
    color_bar.ax.xaxis.set_major_locator(ticker.FixedLocator([0,1,2,3]))
    color_bar.ax.xaxis.set_major_formatter(ticker.FixedFormatter(['$10^0$', '$10^1$', '$10^2$', '$10^3$']))
    color_bar.ax.tick_params(labelsize=15)

    # Ticks
    tick = list(range(0,number_bins+1))
    xticks = np.linspace(min_ch,max_ch,number_bins+1)#list(range(0,max_ch-min_ch+1,1)) 
    yticks = np.linspace(min_tb,max_tb,number_bins+1)#list(range(0,max_tb-min_tb+1,10))
    f_xip = interpolate.interp1d( xticks, tick)
    f_yip = interpolate.interp1d( yticks, tick)
    # plot line where y=y0 is
    ax.axhline( y = f_yip(210), color='grey', linestyle='dashed', linewidth=2)
    # plot a patch
    #cloud_th = 210
    #xc_left = f_xinterp( 180 )
    #xc_right = f_xinterp( cloud_th )
    #yc_bottom = f_yinterp( 15 )
    #yc_up = f_yinterp( 45 )
    #pp = plt.Rectangle( (xc_left,yc_bottom),xc_right-xc_left,yc_up-yc_bottom,alpha=0.2,facecolor='grey')
    #ax[j].add_patch(pp) 
    # labels
    xtick_labels = list(np.arange(min_ch,max_ch+1,2.5))
    ytick_labels = list(range(min_tb,max_tb+1,10))
    ax.set_xticks( [f_xip(it) for it in xtick_labels] )
    ax.set_yticks( [f_yip(it) for it in ytick_labels] )
    ax.set_xticklabels( [str(it) for it in xtick_labels],fontsize=15 )
    ax.set_yticklabels( [str(it) for it in ytick_labels],fontsize=15 )
    ax.set_xlim(xmin=0,xmax=number_bins)
    ax.set_ylim(ymin=0,ymax=number_bins)
    ax.set_xlabel('IR-sensitive cloud top height (km)',fontsize=16)
    ax.set_ylabel('Brightness Temperatures (K)',fontsize=15)


    des_name = '/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/JerryRun/IR_THO/wrf_df/201708241200/twoD_cloud_top_height_'+FCtime+'.png' #big_dir+Storm+'/'+Exper_name+'/cloud_top_height_'+MP+'_'+FCtime+'.png'
    plt.savefig( des_name, dpi=200)
    print('Saving the figure: ', des_name)  

#CRTM_wrfout_d03_2017-08-25_18:00:00.bin 

def spatialD_oneScene( wrf_files, cdirs ):

    # Read WRF domain
    domain = UD.read_wrf_domain( wrf_files[0] )
    ncdir = nc.Dataset(wrf_files[0], 'r')
    xlat = ncdir.variables['XLAT'][0,:,:].flatten()
    xlon = ncdir.variables['XLONG'][0,:,:].flatten()

    # ------------------ Plot -----------------------
    f, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5, sharex='all', sharey='all',  figsize=(5,5), dpi=400)

    # Define the domain
    lat_min = domain['lat_min']
    lat_max = domain['lat_max']
    lon_min = domain['lon_min']
    lon_max = domain['lon_max']

    # cloud depth that IR can see
    if ver_use_press:
        min_vc = 50
        max_vc = 850
        bounds = [100,200,300,400,500,600,700,800,900]
    else:
        bounds = [0,5,10,12,14,15,17,18,20]
        min_vc = np.amin(bounds)
        max_vc = np.amax(bounds)
        exist_cmap = plt.cm.magma.reversed()
        colors = exist_cmap(np.linspace(0,1,len(bounds)))
        new_map = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(bounds))
        #min_vc = 0
        #max_vc = 21
        #bounds = np.linspace(min_vc,max_vc,8)
    if plot_scatter:
        cdir_f = ax.scatter(xlon,xlat,s=2,c=cdirs[0,:], edgecolors='none', cmap=new_map, vmin=min_vc,vmax=max_vc,levels=bounds, transform=ccrs.PlateCarree())
        cbaxes = f.add_axes([0.125, 0.05, 0.775, 0.02])
        #cdir_bar = fig.colorbar(ticks=range(min_vc,max_vc+50,100),cax=cbaxes,fraction=0.046, pad=0.04)
        #cdir_bar.ax.tick_params(labelsize=12)
    else:
        cdir_f = ax.contourf(xlon.reshape(xmax,ymax),xlat.reshape(xmax,ymax),cdirs[0,:].reshape(xmax,ymax),cmap='magma_r',levels=bounds,extend='max',transform=ccrs.PlateCarree())
        caxes = f.add_axes([0.125, 0.05, 0.775, 0.02])
        color_ticks = bounds
        cdir_bar = f.colorbar(cdir_f,orientation="horizontal", cax=caxes,ticks=color_ticks)
        bounds_str =  [ str(item) for item in color_ticks ]
        cdir_bar.ax.set_xticklabels( bounds_str)
        cdir_bar.ax.tick_params(labelsize=10)

    # Mark the slp 
    #ax.scatter(lon_storm,lat_storm,20,'blue',marker='*',transform=ccrs.PlateCarree())

    #title for all
    f.suptitle(Storm+':'+Exper_name+'\n@'+FCtime,fontsize=11, fontweight='bold')
    #f.suptitle(Storm+' IC:Ref;FC:Ref \n@'+fctime,fontsize=11, fontweight='bold')

    #subplot title
    font = {'size':8,}
    #ax.set_title('Forecast--min slp: '+str("{0:.3f}".format(min_slp.values))+' hPa',font,fontweight='bold',fontsize=13)

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
        gl.xlabel_style = {'size': 10}
        gl.ylabel_style = {'size': 10}

    des_name = '/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/JerryRun/IR_THO/wrf_df/201708241200/cloud_top_height_'+FCtime+'.png' #big_dir+Storm+'/'+Exper_name+'/cloud_top_height_'+MP+'_'+FCtime+'.png'
    plt.savefig( des_name, dpi=200)
    print('Saving the figure: ', des_name)



def Plot_depth_IRsee_OneScene( FCtime, wrf_files, d_hydro, ver_coor ):

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

    # plot spatial distribution 
    if plot_spatialD:
        spatialD_oneScene( wrf_files, cdirs )
    if plot_2D_hist:
        cloudHeight_Tb_2D( cdirs )


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    model_resolution = 3000 #m

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'IR'
    MP = 'TuneWSM6'
    hydros =  ['QICE','QSNOW','QGRAUP']#['QCLOUD','QRAIN','QICE','QSNOW','QGRAUP']
    
    sensor = 'abi_gr'
    ch_list = ['8',]
    fort_v = ['obs_type','lat','lon','obs']

    start_time_str = '201709030100'
    end_time_str = '201709030100'
    time_step = 3600 # seconds
    Consecutive_times = True

    deep_slp_incre = False
    incre_slp_th = 0 # threshold of increment, unit:hpa
    plot_circle = True
    radius_th = 200 # km

    path_from_top = True # accumulate water from top to down
    accu_th = 20.0 # gram per meter squared 
    ver_use_press = False
    each_water = True
    domain_mean = True

    depth_IR_see = False
    if depth_IR_see:
        oneScene = True
    plot_scatter = False
    plot_spatialD = False
    plot_2D_hist = True

    # ------------------------------------------------------  
    # Dimension of the domain
    nLevel = 42
    xmax = 297
    ymax = 297

    # Create experiment names
    Exper_name = UD.generate_one_name( Storm,DA,MP ) #'JerryRun/IR_THO/wrf_df/201708241200'

    # Times
    if not Consecutive_times:
        FCtimes = ['201709180000',]
    else:
        time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / time_step
        time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+1, 3))]
        FCtimes = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]


    if not deep_slp_incre:
        # ------ Plot -------------------
        plot_dir = small_dir+Storm+'/'+Exper_name+'/Vis_analyze/Model/Accu_HydroPath/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)

        start_time=time.process_time()
        for FCtime in FCtimes:
            idx_t = FCtimes.index( FCtime )
            print('At '+FCtime)
            ## deterministic forecast
            #wrfname = 'wrfout_d03_'+str(FCtime[0:4])+'-'+str(FCtime[4:6])+'-'+str(FCtime[6:8])+'_'+str(FCtime[8:10])+':00:00'
            #wrf_files = ['/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/JerryRun/IR_THO/wrf_df/201708241200/'+wrfname,]
            wrf_dir = big_dir+Storm+'/'+Exper_name+'/fc/'+FCtime+'/'
            wrf_files = [wrf_dir+'/wrf_enkf_input_d03_mean',wrf_dir+'/wrf_enkf_output_d03_mean']
            d_hydro,ver_coor = compute_accu_hydromass( wrf_files, hydros, idx_t) # hydro mass at any grid point

            # --- Condition
            if each_water:
                v_interest = hydros
            else:
                v_interest = ['liquid','ice','all_hydro']
            
            if domain_mean:
                Plot_DM_all_oneTime( FCtime,v_interest,d_hydro,ver_coor ) 

            if depth_IR_see:
                if not oneScene:
                    Plot_depth_IRsee( FCtime, wrf_files, d_hydro, ver_coor )
                else:
                    Plot_depth_IRsee_OneScene( FCtime, wrf_files, d_hydro, ver_coor )

        end_time = time.process_time()
        print ('time needed: ', end_time-start_time, ' seconds')




 

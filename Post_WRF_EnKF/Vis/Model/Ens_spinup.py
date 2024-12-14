import os,sys,stat # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import math
import netCDF4 as nc
import matplotlib
import scipy as sp
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import numpy.ma as ma

import Util_data as UD

# ---------------------------------------------------------------------------------------------------------------
#    Operation: identify min PSFC for each ensemble member
# ------------------------------------------------------------------------------------------------------------------
# Read and idenfity the minimum PSFC in each ensemble member
def find_mpsfc( wrf_file ):

    ncdir = nc.Dataset( wrf_file )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    PSFC = ncdir.variables['PSFC'][0,:,:]

    mpsfc = np.nanmin( PSFC )/100 
    idx = np.nanargmin( PSFC )
    lat_mpsfc =  lat.flatten()[idx] 
    lon_mpsfc =  lon.flatten()[idx] 

    return {'mpsfc':mpsfc,'lat_mpsfc':lat_mpsfc,'lon_mpsfc':lon_mpsfc}

# Loop thru all ensemble members &
# Identify the min PSFC for that member &
# Save mPSFCs to a txt file
def identify_mpsfc_ens( wrf_dir, key=None):

    mpsfc_info = np.zeros( shape=(num_ens,4) )
    #mpsfc_info = np.nan

    for ie in range(num_ens):
        id_num = ie+1
        if key == 'start':
            wrf_file = wrf_dir+'wrfinput_d03_'+f"{ie+1:03}"
        else:
            dt = datetime.strptime(end_time_str[Storm], '%Y%m%d%H%M')
            wrf_time = dt.strftime('%Y-%m-%d_%H:%M:%S')
            wrf_file = wrf_dir+'wrfinput_d03_'+wrf_time+'_'+f"{ie+1:03}"
        d_mpsfc = find_mpsfc( wrf_file)
        mpsfc_info[ie,1] = "{0:.3f}".format( d_mpsfc['lon_mpsfc'] )
        mpsfc_info[ie,2] = "{0:.3f}".format(d_mpsfc['lat_mpsfc'] )
        mpsfc_info[ie,3] = "{0:.3f}".format(d_mpsfc['mpsfc'] )

    # May save the mpsfc info
    if If_save:
        header = ['ID','Lon','Lat','mpsfc']
        if key == 'start':
            des_path = wrf_dir+ start_time_str[Storm]+'_ens_mpsfc.txt'
        else:
            des_path = wrf_dir+ end_time_str[Storm]+'_ens_mpsfc.txt'
        with open(des_path,'w') as f:
            # Add header 
            f.write('\t'.join( item.rjust(6) for item in header ) + '\n' )
            # Write the record to the file serially
            len_records = np.shape( mpsfc_info )[0]
            for irow in range( len_records ):
                irecord = mpsfc_info[irow,:].astype(str)
                irecord[0] = f"{irow+1:03}"
                #print(irecord)
                f.write('\t'.join( item.rjust(6) for item in irecord ) + '\n')
    print('Save '+des_path)

    return None

# Read Ens mpsfc info from pre-calculated files
# e.g., DAtime/DAtime_enkf_output_mpsfc.txt 
def read_mpsfc_ens( txt_file  ):

    d_mpsfc = {}

    with open( txt_file ) as f:
        next(f)
        all_lines = f.readlines()

    for line in all_lines:
        split_line = line.split()
        id_ = split_line[0]
        lon = float(split_line[1])
        lat = float(split_line[2])
        mpsfc = float(split_line[3])
        d_mpsfc[id_] = (lon,lat,mpsfc)

    return d_mpsfc

# ---------------------------------------------------------------------------------------------------------------
#    Operation: identify min SLP for each ensemble member
# ------------------------------------------------------------------------------------------------------------------
def find_mslp( Storm,wrf_file,time_str):

    ncdir = nc.Dataset( wrf_file )
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # sea level pressure
    #slp = getvar(ncdir, 'slp')
    slp = UD.compute_slp( ncdir )
    #print( np.amax( abs(slp - slp_getvar.values )) )
    # original SLP
    slp_values = slp #slp.values
    slp_values[slp_values > 1030] = np.nan
    slp_original = slp_values
    # smoothed SLP
    slp_smt_values = sp.ndimage.gaussian_filter(slp, [1,1]) #[11,11]
    slp_smt_values[slp_smt_values > 1030] = np.nan
    slp_smooth = slp_smt_values
    # simulated storm center
    if Storm == 'HARVEY': # Harvey is special with its location near land!

        # eyeball where the storm is
        if time_str == '201708220000':
            lat_mask = lat <= 16.6
            lon_mask = lon <= -88.
            mask = lat_mask | lon_mask
        elif (time_str >= '201708221200') & (time_str <= '201708221200'):
            lat_mask = lat <= 18
            lon_mask = lon <= -91.5
            mask = lat_mask | lon_mask
        elif (time_str >= '201708221700') & (time_str <= '201708222300'):
            lat_mask = lat <= 18
            lon_mask =  (lon >= -88) | (lon <= -91.5)
            mask = lat_mask | lon_mask
        else:
            mask = lat <= 18
        slp_masked = ma.masked_array(slp_values, mask=mask)
        minslp = np.nanmin( slp_masked )

        slp_smooth_masked = ma.masked_array(slp_smt_values, mask=mask)
        idx = np.nanargmin( slp_smooth_masked )
        lat_minslp = lat.flatten()[idx]
        lon_minslp = lon.flatten()[idx]
        #print(lon_minslp,lat_minslp)
    else:
        minslp = np.nanmin( slp_values )
        idx = np.nanargmin( slp_smt_values )
        lat_minslp = lat.flatten()[idx]
        lon_minslp = lon.flatten()[idx]
    return {'mslp':minslp,'lat_minslp':lat_minslp,'lon_minslp':lon_minslp}

# Loop thru all ensemble members &
# Identify the mslp for that member &
# Save mslps to a txt file
def identify_mslp_ens( wrf_dir, key=None):

    mslp_info = np.zeros( shape=(num_ens,4) )
    #mslp_info = np.nan

    for ie in range(num_ens):
        id_num = ie+1
        if key == 'start':
            wrf_file = wrf_dir+'wrfinput_d03_'+f"{ie+1:03}"
            d_mslp = find_mslp( Storm,wrf_file,start_time_str[Storm])
        else:
            dt = datetime.strptime(end_time_str[Storm], '%Y%m%d%H%M')
            wrf_time = dt.strftime('%Y-%m-%d_%H:%M:%S')
            wrf_file = wrf_dir+'wrfinput_d03_'+wrf_time+'_'+f"{ie+1:03}"
            d_mslp = find_mslp( Storm,wrf_file,end_time_str[Storm])
        
        mslp_info[ie,1] = "{0:.3f}".format( d_mslp['lon_minslp'] )
        mslp_info[ie,2] = "{0:.3f}".format(d_mslp['lat_minslp'] )
        mslp_info[ie,3] = "{0:.3f}".format(d_mslp['mslp'] )

    # May save the mslp info
    if If_save:
        header = ['ID','Lon','Lat','mslp']
        if key == 'start':
            des_path = wrf_dir+ start_time_str[Storm]+'_ens_mslp.txt'
        else:
            des_path = wrf_dir+ end_time_str[Storm]+'_ens_mslp.txt'
        with open(des_path,'w') as f:
            # Add header 
            f.write('\t'.join( item.rjust(6) for item in header ) + '\n' )
            # Write the record to the file serially
            len_records = np.shape( mslp_info )[0]
            for irow in range( len_records ):
                irecord = mslp_info[irow,:].astype(str)
                irecord[0] = f"{irow+1:03}"
                #print(irecord)
                f.write('\t'.join( item.rjust(6) for item in irecord ) + '\n')
    print('Save '+des_path)

    return None

# Read Ens mslp info from pre-calculated files
# e.g., DAtime/DAtime_enkf_output_mslp.txt 
def read_mslp_ens( txt_file  ):

    d_mslp = {}

    with open( txt_file ) as f:
        next(f)
        all_lines = f.readlines()

    for line in all_lines:
        split_line = line.split()
        id_ = split_line[0]
        lon = float(split_line[1])
        lat = float(split_line[2])
        mslp =  float(split_line[3])
        d_mslp[id_] = (lon,lat,mslp)

    return d_mslp

# ---------------------------------------------------------------------------------------------------------------
#    Operation: hard-to-categorize functions
# ------------------------------------------------------------------------------------------------------------------
# Read var for the ens for experiments
def Read_vars_Ens( var_name ):
    d_st = {}
    d_ed = {}

    for imp in MP:
        d_st[imp] = {}
        d_ed[imp] = {}
        wrf_dir = big_dir+Storm+'/'+Exper_name[imp]+'/fc/'+start_time_str[Storm]+'/'
        if var_name == 'PSFC':
            print('Reading the minimum PSFC for each ens member...')
            txt_file = wrf_dir+ start_time_str[Storm]+'_ens_mpsfc.txt'
            d_st[imp] = read_mpsfc_ens( txt_file )
            txt_file = wrf_dir+ end_time_str[Storm]+'_ens_mpsfc.txt'
            d_ed[imp] = read_mpsfc_ens( txt_file )
        elif var_name == 'slp':
            print('Reading the mslp for each ens member...')
            txt_file = wrf_dir+ start_time_str[Storm]+'_ens_mslp.txt'
            d_st[imp] = read_mslp_ens( txt_file )
            txt_file = wrf_dir+ end_time_str[Storm]+'_ens_mslp.txt'
            d_ed[imp] = read_mslp_ens( txt_file )

    return d_st, d_ed

def read_wrf_domain( wrf_file ):

    ncdir = nc.Dataset(wrf_file, 'r')

    Lat_x = ncdir.variables['XLAT'][0,:,0] #latitude: XLAT(time, y, x)
    Lon_x = ncdir.variables['XLONG'][0,0,:] #longitude: XLONG(time, y, x)

    lat_min = np.min( Lat_x.flatten() )
    lat_max = np.max( Lat_x.flatten() )
    lon_min = np.min( Lon_x.flatten() )
    lon_max = np.max( Lon_x.flatten() )

    d_wrf_d = {'lat_min':lat_min,'lat_max':lat_max, 'lon_min':lon_min, 'lon_max':lon_max}
    return d_wrf_d

def identify_one_domain():

    lats_min = []
    lats_max = []
    lons_min = []
    lons_max = []
    for imp in MP:
        wrf_dir = big_dir+Storm+'/'+Exper_name[imp]+'/fc/'+start_time_str[Storm]+'/'
        dt = datetime.strptime(end_time_str[Storm], '%Y%m%d%H%M')
        end_time = dt.strftime('%Y-%m-%d_%H:%M:%S')
        wrf_files = [wrf_dir+'wrfinput_d03_001',wrf_dir+'wrfinput_d03_'+end_time+'_001']
        for ifile in wrf_files:
            dd = read_wrf_domain( ifile )
            lats_min.append(dd['lat_min'] )
            lats_max.append(dd['lat_max'] )
            lons_min.append(dd['lon_min'] )
            lons_max.append(dd['lon_max'] )
    min_lat = min(lats_min)
    min_lon = min(lons_min)
    max_lat = max(lats_max)
    max_lon = max(lons_max)

    return {'min_lat':min_lat,'min_lon':min_lon,'max_lat':max_lat,'max_lon':max_lon}

def identify_domain( key ):

    lats_min = []
    lats_max = []
    lons_min = []
    lons_max = []
    for imp in MP:
        wrf_dir = big_dir+Storm+'/'+Exper_name[imp]+'/fc/'+start_time_str[Storm]+'/'
        dt = datetime.strptime(end_time_str[Storm], '%Y%m%d%H%M')
        end_time = dt.strftime('%Y-%m-%d_%H:%M:%S')
        if key == 'input':
            wrf_file =  wrf_dir+'wrfinput_d03_001'
        else: 
            wrf_file =  wrf_dir+'wrfinput_d03_'+end_time+'_001' 
        dd = read_wrf_domain( wrf_file )
        lats_min.append(dd['lat_min'] )
        lats_max.append(dd['lat_max'] )
        lons_min.append(dd['lon_min'] )
        lons_max.append(dd['lon_max'] )
    min_lat = min(lats_min)
    min_lon = min(lons_min)
    max_lat = max(lats_max)
    max_lon = max(lons_max)

    return {'min_lat':min_lat,'min_lon':min_lon,'max_lat':max_lat,'max_lon':max_lon}

def blend_colors(color1, color2):
    # Convert hex to RGB
    rgb1 = mcolors.to_rgb(color1)
    rgb2 = mcolors.to_rgb(color2)

    # Take the mean of each RGB channel
    blended_rgb = [(c1 + c2) / 2 for c1, c2 in zip(rgb1, rgb2)]

    # Convert back to hex
    blended_hex = mcolors.to_hex(blended_rgb)
    return blended_hex


# ---------------------------------------------------------------------------------------------------------------
#    Operation: Plot movement of storm center in geomap
# ------------------------------------------------------------------------------------------------------------------
def plot_geo_move( d_st, d_ed ):

    # Find the domain area
    d_domain = identify_one_domain()
    if Storm == 'IRMA':
        lat_min = d_domain['min_lat']+3
        lon_min = d_domain['min_lon']+3
        lat_max = d_domain['max_lat']-4
        lon_max = d_domain['max_lon']-5
    elif Storm == 'JOSE':
        lat_min = d_domain['min_lat']
        lon_min = d_domain['min_lon']-1
        lat_max = d_domain['max_lat']
        lon_max = d_domain['max_lon']+0.5
        

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5,dpi=300)
    #fig, ax=plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0, 'hspace':0}, linewidth=0.5,figsize=(6.5,6), dpi=300) #(6.5,6)
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black',linewidth=0.5)

    # Normalize the variable for colormap
    if Storm == 'IRMA':
        norm = plt.Normalize(950, 1020)
    elif Storm == 'JOSE':
        norm = plt.Normalize(1008, 1015)
    cmap= {'WSM6':cm.Reds_r,'THO':cm.Blues_r}

    for imp in MP:
        # loop thru ids
        for im in d_st[imp].keys():
            # plot the start and end points
            lon_st = d_st[imp][im][0]
            lon_ed = d_ed[imp][im][0]
            lat_st = d_st[imp][im][1]
            lat_ed = d_ed[imp][im][1]
            color_st = cmap[imp](norm(d_st[imp][im][2]))
            color_ed = cmap[imp](norm(d_ed[imp][im][2]))
            ax.plot(lon_st,lat_st, 's', color=color_st, markersize=2, transform=ccrs.PlateCarree())
            ax.plot(lon_ed,lat_ed, 'o', color=color_ed, markersize=2, transform=ccrs.PlateCarree())
            # draw an arrow connecting the points
            #mean_color = blend_colors(color_st, color_ed)
            ax.annotate('', xy=(lon_ed, lat_ed), xytext=(lon_st, lat_st),
                   arrowprops=dict(arrowstyle="->", color=color_st, lw=0.5), transform=ccrs.PlateCarree())

    # Add a colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap['THO'], norm=norm)
    sm.set_array([])  # Dummy array for colorbar
    cax_above = fig.add_axes([0.15, 0.1, 0.7, 0.02])  # [left, bottom, width, height]
    cbar1 = fig.colorbar(sm, cax=cax_above, orientation='horizontal',label='THO', extend='both')

    sm = plt.cm.ScalarMappable(cmap=cmap['WSM6'], norm=norm)
    sm.set_array([])  # Dummy array for colorbar
    cax_below = fig.add_axes([0.15, 0.95, 0.7, 0.02])  # [left, bottom, width, height]
    cbar2 = fig.colorbar(sm, cax=cax_below, orientation='horizontal',label='WSM6', extend='both')
    cbar2.set_ticks([])

    #subplot title
    font = {'size':12,}
    if var_name == 'slp':
        title_name = 'Minimum SLP (hPa)'
    elif var_name == 'PSFC':
        title_name = 'Minimum PSFC (hPa)'
    ax.set_title( title_name, font, fontweight='bold')

    #title for all
    tt_name = Storm+': movement over 12-hr spinup'
    fig.suptitle(tt_name, fontsize=10)

    # Axis labels
    lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
    lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),1))
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    gl.top_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.right_labels = False

    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 12}

    if var_name == 'slp':
        save_des = small_dir+Storm+'/'+Exper_name[MP[0]]+'/Vis_analyze/Model/Ens_MSLP_time/Ens_spinup_MSLP_geomove.png'
    elif var_name == 'PSFC':
        save_des = small_dir+Storm+'/'+Exper_name[MP[0]]+'/Vis_analyze/Model/Ens_MSLP_time/Ens_spinup_MPSFC_geomove.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()


# ---------------------------------------------------------------------------------------------------------------
#    Operation: Plot storm center for input and output
# ------------------------------------------------------------------------------------------------------------------
def plot_geo_compare( d_st, d_ed ):

    fig, ax=plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, gridspec_kw = {'wspace':0.1, 'hspace':0.5}, linewidth=0.5,figsize=(6.5,4.5), dpi=300)
    for i in range(2):
        ax[i].coastlines(resolution='10m', color='black',linewidth=0.5)
    # Set domain
    if Storm == 'IRMA':
        lat_min_st = 18
        lat_max_st = 20
        lon_min_st = -44
        lon_max_st = -42
    else:
        d_domain = identify_domain( 'input')
        lat_min_st = d_domain['min_lat']-0.5
        lon_min_st = d_domain['min_lon']-0.5
        lat_max_st = d_domain['max_lat']+0.5
        lon_max_st = d_domain['max_lon']+0.5
    ax[0].set_extent([lon_min_st,lon_max_st,lat_min_st,lat_max_st], crs=ccrs.PlateCarree())
    if Storm == 'IRMA':
        lat_min_ed = 17
        lat_max_ed = 20
        lon_min_ed = -47
        lon_max_ed = -43
    else:
        d_domain = identify_domain( 'output')
        lat_min_ed = d_domain['min_lat']-0.5
        lon_min_ed = d_domain['min_lon']-0.5
        lat_max_ed = d_domain['max_lat']+0.5
        lon_max_ed = d_domain['max_lon']+0.5
    ax[1].set_extent([lon_min_ed,lon_max_ed,lat_min_ed,lat_max_ed], crs=ccrs.PlateCarree())

    # Normalize the variable for colormap
    if Storm == 'IRMA':
        norm = plt.Normalize(950, 1020)
    else:
        norm = plt.Normalize(1008, 1015)
    cmap= {'WSM6':cm.Reds_r,'THO':cm.Blues_r}

    # Scatter ensemble
    for imp in MP:
        # loop thru ids
        for im in d_st[imp].keys():
            lon_st = d_st[imp][im][0]
            lon_ed = d_ed[imp][im][0]
            lat_st = d_st[imp][im][1]
            lat_ed = d_ed[imp][im][1]
            color_st = cmap[imp](norm(d_st[imp][im][2]))
            color_ed = cmap[imp](norm(d_ed[imp][im][2]))
            # scatter the start points
            ax[0].plot(lon_st,lat_st, 's', color=color_st, markersize=2, transform=ccrs.PlateCarree())
            # scatter the end points
            ax[1].plot(lon_ed,lat_ed, 's', color=color_ed, markersize=2, transform=ccrs.PlateCarree())
    
    # Scatter TCvital
    tc_lon_st, tc_lat_st, tc_slp_st = UD.read_TCvitals(small_dir, Storm, start_time_str[Storm])
    ax[0].plot(tc_lon_st,tc_lat_st, '*', color='black', markersize=5, transform=ccrs.PlateCarree())
    tc_lon_ed, tc_lat_ed, tc_slp_ed = UD.read_TCvitals(small_dir, Storm, end_time_str[Storm])
    ax[1].plot(tc_lon_ed,tc_lat_ed, '*', color='black', markersize=5, transform=ccrs.PlateCarree())
    
    fig.text(0.5,0.85,'*: TCvital',ha='center',va='center',rotation='horizontal',fontsize=11)

    # Add a colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap['THO'], norm=norm)
    sm.set_array([])  # Dummy array for colorbar
    cax_above = fig.add_axes([0.15, 0.1, 0.7, 0.02])  # [left, bottom, width, height]
    if var_name == 'slp':
        label_name = 'THO: Minimum SLP (hPa)'
    elif var_name == 'PSFC':
        label_name = 'THO: Minimum PSFC (hPa)'
    fig.colorbar(sm, cax=cax_above, orientation='horizontal',label=label_name, extend='both')

    sm = plt.cm.ScalarMappable(cmap=cmap['WSM6'], norm=norm)
    sm.set_array([])  # Dummy array for colorbar
    cax_below = fig.add_axes([0.15, 0.15, 0.7, 0.02])  # [left, bottom, width, height]
    cbar2 = fig.colorbar(sm, cax=cax_below, orientation='horizontal',label='WSM6', extend='both')
    cbar2.set_ticks([])

    #subplot title
    font = {'size':12,}
    ax[0].set_title( 'Start Point', font, fontweight='bold')
    ax[1].set_title( 'End Point', font, fontweight='bold')

    #title for all
    tt_name = Storm+': movement over 12-hr spinup'
    fig.suptitle(tt_name, fontsize=15)

    # Axis labels
    def axis_labels( lon_min,lat_min,lon_max,lat_max,gl):
        if Storm == 'IRMA':
            lon_ticks = list(range(math.ceil(lon_min)-1, math.ceil(lon_max)+1,1))

        else:
            lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
        lat_ticks = list(range(math.ceil(lat_min), math.ceil(lat_max),1))
        gl.top_labels = False
        gl.bottom_labels = True
        gl.left_labels = True
        gl.right_labels = False
        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10}
        gl.ylabel_style = {'size': 12}

    # plot axis labels
    gl = ax[0].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    axis_labels( lon_min_st,lat_min_st,lon_max_st,lat_max_st, gl )
    gl = ax[1].gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    axis_labels( lon_min_ed,lat_min_ed,lon_max_ed,lat_max_ed, gl )

    if var_name == 'slp':
        save_des = small_dir+Storm+'/'+Exper_name[MP[0]]+'/Vis_analyze/Model/Ens_MSLP_time/Ens_spinup_MSLP_geo_compare.png'
    elif var_name == 'PSFC':
        save_des = small_dir+Storm+'/'+Exper_name[MP[0]]+'/Vis_analyze/Model/Ens_MSLP_time/Ens_spinup_MPSFC_geo_compare.png'
    plt.savefig( save_des )
    print( 'Saving the figure: ', save_des )
    plt.close()



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'#'/expanse/lustre/scratch/zuy121/temp_project/Pro2_PSU_MW/' #'/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'#'/expanse/lustre/projects/pen116/zuy121/Pro2_PSU_MW/'  #'/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storm = 'IRMA'
    DA = 'CONV'
    MP = ['THO','WSM6']
    #fort_v = ['obs_type','lat','lon','obs']

    # observation type
    #obs_assimilated = True
    #if obs_assimilated:
    #    obs_type = 'slp' # Radiance

    # model variable
    if Storm == 'HARVEY':
        model_v = ['slp',]
    else:
        model_v = [ 'PSFC',]#'QSNOW','QCLOUD','QRAIN','QICE','QGRAUP']

    # time
    start_time_str = {'HARVEY':'201708220000','IRMA':'201709021200','JOSE':'201709041200','MARIA':'201709151200'}
    end_time_str =  {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}

    key = 'end' # 'end'

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    # if calculate data
    calculate_ens_data = False
    If_save = True

    # plot
    if_plot_geo = True

    # -------------------------------------------------------
    Exper_name = {}
    for imp in MP:
        Exper_name[imp] = UD.generate_one_name( Storm,DA,imp )

    # Calculate more data
    if calculate_ens_data:
        for imp in MP:
            wrf_dir = big_dir+Storm+'/'+Exper_name[imp]+'/fc/'+start_time_str[Storm]+'/'
            for var_name in model_v:
                if var_name == 'PSFC':
                    print('Finding the minimum PSFC for each ens member...')
                    identify_mpsfc_ens( wrf_dir, key )   
                elif var_name == 'slp':
                    print('Finding the mslp for each ens member...')
                    identify_mslp_ens( wrf_dir, key )

    # Plot
    for var_name in model_v:
        d_st, d_ed = Read_vars_Ens( var_name ) 

        if if_plot_geo:
            plot_geo_compare( d_st, d_ed )
            #plot_geo_move( d_st, d_ed )




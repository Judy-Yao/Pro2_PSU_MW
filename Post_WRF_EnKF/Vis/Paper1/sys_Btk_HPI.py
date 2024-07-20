
import numpy as np
import matplotlib
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import time
import matplotlib.dates as mdates
from matplotlib import pyplot
import matplotlib.lines as mlines
from matplotlib.colors import Normalize
from datetime import datetime, timedelta
from shapely.geometry import Polygon
import matplotlib.patches as patches

import Util_data as UD



# layout
# rows: storms
# columns: track, HPI
def plot_btk():

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=300)
    outer_grid = fig.add_gridspec(ncols=1,nrows=4,top=0.95,left=0.1,)

    # Create a colormap for track; marker colors for Vmax
    cmap = cm.get_cmap('nipy_spectral')
    norm = Normalize(vmin=0, vmax=90)

    ax = {}
    ax_mslp = {}
    for ist in Storms:
        ax[ist] = {}
        ir = Storms.index( ist )
        # gridspec inside gridspec
        inner_grid = outer_grid[ir].subgridspec(1, 2)
        ax[ist]['track'] = fig.add_subplot( inner_grid[0], projection=ccrs.PlateCarree())
        d_ist = domain[ist]
        ax[ist]['track'].set_extent( [d_ist['lon_min']-0.01,d_ist['lon_max']+0.01,d_ist['lat_min']-0.01,d_ist['lat_max']+0.01],crs=ccrs.PlateCarree())
        ax[ist]['track'].coastlines( resolution='10m', color='black',linewidth=0.5 )
        if ist =='HARVEY':
            ax[ist]['track'].set_aspect(0.7)
        elif ist == 'IRMA':
            ax[ist]['track'].set_aspect(2.8)
        elif ist == 'MARIA':
            ax[ist]['track'].set_aspect(2.5)
        elif ist == 'JOSE':
            ax[ist]['track'].set_aspect(2.3)
        ax[ist]['its'] = fig.add_subplot( inner_grid[1] )

        # Find onset of RI
        time = d_btk[ist]['time']
        lon = d_btk[ist]['lon']
        lat = d_btk[ist]['lat']
        lon_ri = lon[time == onse_RI[ist]]
        lat_ri = lat[time == onse_RI[ist]]
        
        # Plot Track
        ax[ist]['track'].plot(d_btk[ist]['lon'],d_btk[ist]['lat'],color='gray',linewidth=1.5,linestyle='-',transform=ccrs.PlateCarree())
        # colorcode markers
        colors = cmap(norm(d_btk[ist]['max_ws']))
        max_marker = 5
        min_marker = 2.5
        slope = (max_marker -min_marker)/(90-0)
        for i in range(len(d_btk[ist]['max_ws'])):
            color = cmap(norm(d_btk[ist]['max_ws'][i]))
            if time[i] == onse_RI[ist]:
                marker = '+'
                color = 'black'
                markersize = 5
            else:
                marker = 'o'
                markersize=min_marker+slope*d_btk[ist]['max_ws'][i]
            ax[ist]['track'].plot(d_btk[ist]['lon'][i],d_btk[ist]['lat'][i],marker=marker,
                markersize=markersize,markeredgecolor=color,markerfacecolor=color,transform=ccrs.PlateCarree())
        # Mark EnKF window
        da_st = EnKF_wd[ist]['st']
        da_end = EnKF_wd[ist]['end']
        lon_st = lon[time == da_st]
        lon_end = lon[time == da_end]
        lat_st = lat[time == da_st]
        lat_end = lat[time == da_end]
        # create a Polygon from the corners
        polygon = Polygon([(lon_st,lat_st),(lon_end,lat_st),(lon_end,lat_end),(lon_st,lat_end),(lon_st,lat_st)])
        mid_top_x = (lon_end + lon_st) / 2
        # Add the polygon to the map
        ax[ist]['track'].add_geometries([polygon], ccrs.PlateCarree(), linewidth=2, edgecolor='red', facecolor='none')
        # Annotate above the polygon
        if ist == 'IRMA':
            mid_top_y = lat_st
        else:
            mid_top_y = lat_end
        ax[ist]['track'].annotate('EnKF', xy=(mid_top_x, mid_top_y), xytext=(mid_top_x, mid_top_y+0.5),
            color='black', fontweight='bold',fontsize=8, ha='center', transform=ccrs.PlateCarree())

        # Plot intensity
        ax_mslp[ist] = ax[ist]['its'].twinx()
        dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in d_btk[ist]['time']]
        # Saffir-Simpson scale
        ax[ist]['its'].fill_between([dates[0], dates[-1]], 0, 90, color='#A12830')
        ax[ist]['its'].text(dates[-5], 80, 'CAT 5', fontsize=6, fontweight='bold', color='blue')
        ax[ist]['its'].fill_between([dates[0], dates[-1]], 0, 70, color='#DB1F2A')
        ax[ist]['its'].text(dates[-5], 64, 'CAT 4', fontsize=6, fontweight='bold', color='blue')
        ax[ist]['its'].fill_between([dates[0], dates[-1]], 0, 58, color='#DB4E4E')
        ax[ist]['its'].text(dates[-5], 53.5, 'CAT 3', fontsize=6, fontweight='bold', color='blue')
        ax[ist]['its'].fill_between([dates[0], dates[-1]], 0, 49, color='#E26E6E')
        ax[ist]['its'].text(dates[-5], 45.5, 'CAT 2', fontsize=6, fontweight='bold', color='blue')
        ax[ist]['its'].fill_between([dates[0], dates[-1]], 0, 42, color='#F1B3B3')
        ax[ist]['its'].text(dates[-5], 37, 'CAT 1', fontsize=6, fontweight='bold', color='blue')
        ax[ist]['its'].fill_between([dates[0], dates[-1]], 0, 32, color='#688FAD')
        ax[ist]['its'].text(dates[-6], 24.5, 'TC Storm', fontsize=6, fontweight='bold', color='blue')
        ax[ist]['its'].fill_between([dates[0], dates[-1]], 0, 17, color='#9FC1D3')
        ax[ist]['its'].text(dates[-8], 8.5, 'TC Depression', fontsize=6, fontweight='bold', color='blue')
        # plot Vmax
        ax[ist]['its'].plot(dates,d_btk[ist]['max_ws'],color='black',linewidth=2.5,linestyle='--',label='Vmax',zorder=1)
        ax[ist]['its'].set_ylabel('Vmax($\mathregular{ms^{-1}}$)',fontsize=7)
        # plot MSLP
        ax_mslp[ist].plot(dates,d_btk[ist]['min_slp'],color='black',linewidth=2.5,linestyle='-',marker=None,label='MSLP',zorder=2)
        ax_mslp[ist].set_ylabel( 'MSLP(hPa)',fontsize=7)
        # Mark onset of RI
        time_ri = datetime.strptime(onse_RI[ist],"%Y%m%d%H%M")
        ax[ist]['its'].scatter(time_ri,d_btk[ist]['max_ws'][time==onse_RI[ist]],s=40,color='white',marker='+',zorder=3)
        # Mark EnKF window
        da_st = datetime.strptime(da_st,"%Y%m%d%H%M")
        da_end = datetime.strptime(da_end,"%Y%m%d%H%M")
        corners = [(da_st,0),(da_st,90),(da_end,90),(da_end,0)] 
        # convert datetime objects to numerical format for plotting
        corners_num = [(mdates.date2num(c[0]), c[1]) for c in corners]
        polygon = patches.Polygon(corners_num,closed=True,linewidth=2,edgecolor='black',facecolor='gray',alpha=0.3)
        ax[ist]['its'].add_patch(polygon)
        ax[ist]['its'].text(dates[3], 25, 'EnKF', fontsize=8, fontweight='bold', color='black')

    # Add colorbar
    cbar_ax = fig.add_axes([0.10, 0.07, 0.25, 0.01])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks([0, 30, 60, 90])
    cbar.set_ticklabels(['0','30','60','90'])
    fig.text(0.42,0.075,'Vmax (m $\mathregular{s^{-1}}$)', fontsize=9, ha='center', va='center',rotation='horizontal')

    # get the position of the subplot
    pos = ax[Storms[-1]]['its'].get_position()
    line_y = ax[Storms[-1]]['its'].get_position().y0-0.04
    line = plt.Line2D([pos.x0+0.02, pos.x0+0.09],[line_y, line_y],transform=fig.transFigure,color='black',linewidth=2.5,linestyle='-')
    fig.lines.append(line)
    fig.text(pos.x0+0.14,line_y,'MSLP', fontsize=9, ha='center', va='center')

    line = plt.Line2D([pos.x1-0.15, pos.x1-0.08],[line_y, line_y],transform=fig.transFigure,color='black',linewidth=2.5,linestyle='--')
    fig.lines.append(line)
    fig.text(pos.x1-0.03,line_y,'Vmax', fontsize=9, ha='center', va='center')


    # Set Minor ticks/labels for track subplot
    for ist in Storms:
        d_ist = domain[ist]
        lon_ticks =  np.arange( d_ist['lon_min'],d_ist['lon_max']+1, 1)
        lat_ticks = np.arange( d_ist['lat_min'],d_ist['lat_max']+1, 1)
        gl = ax[ist]['track'].gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=0.8, color='gray', alpha=0.3, linestyle='--')
        gl.left_labels = False
        gl.bottom_labels = False
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.ylocator = mticker.FixedLocator(lat_ticks)

    # Set Major ticks/labels for track subplot
    for ist in Storms:
        d_ist = domain[ist]
        lon_ticks =  np.arange( d_ist['lon_min'],d_ist['lon_max']+5, 5)
        lat_ticks = np.arange( d_ist['lat_min'],d_ist['lat_max']+5, 5)
        gl = ax[ist]['track'].gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=1, color='gray', alpha=0.5, linestyle='-')
        gl.left_labels = True
        gl.bottom_labels = True
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.yformatter = LATITUDE_FORMATTER
        gl.xformatter = LONGITUDE_FORMATTER
        gl.xlabel_style = {'size': 8}
        gl.ylabel_style = {'size': 8}

    # Set ticks/other attributes for intensity subplots
    for ist in Storms:
        date_form = mdates.DateFormatter("%m-%d")
        Btk_start = start_time_str[ist]
        Btk_end = end_time_str[ist]
        ax[ist]['its'].set_xlim([datetime(int(Btk_start[0:4]),int(Btk_start[4:6]), int(Btk_start[6:8]), int(Btk_start[8:10])), datetime(int(Btk_start[0:4]),int(Btk_end[4:6]), int(Btk_end[6:8]), int(Btk_end[8:10]))])
        ax[ist]['its'].xaxis.set_major_locator(mdates.DayLocator())
        ax[ist]['its'].xaxis.set_major_formatter(date_form)
        ax[ist]['its'].tick_params(axis='x', labelrotation=20,labelsize=6)
        ax_mslp[ist].set_ylim([900,1025])     #([940, 1015])
        ax[ist]['its'].set_ylim([0,90])        
   

    # Add storm info
    for ist in Storms:
        if Storms.index(ist) == 0:
            fig.text(0.02,0.86,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 1:
            fig.text(0.02,0.64,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 2:
            fig.text(0.02,0.42,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 3:
            fig.text(0.02,0.19,ist, fontsize=12, ha='center', va='center',rotation='vertical')


    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_Btk_storms.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    start_time_str = {'HARVEY':'201708220000','IRMA':'201709021200','JOSE':'201709041200','MARIA':'201709151200'}
    end_time_str = {'HARVEY':'201708270000','IRMA':'201709080000','JOSE':'201709100000','MARIA':'201709210000'}
    onse_RI = {'HARVEY':'201708231200','IRMA':'201709040600','JOSE':'201709060600','MARIA':'201709170000'}
    EnKF_wd = {'HARVEY':{'st':'201708221200','end':'201708231200'},
            'IRMA':{'st':'201709030000','end':'201709040000'},
            'JOSE':{'st':'201709050000','end':'201709060000'},
            'MARIA':{'st':'201709160000','end':'201709170000'},
    } #window

    domain = {'HARVEY':{'lon_min':-100,'lon_max':-85,'lat_min':15,'lat_max':30},
            'IRMA':{'lon_min':-75,'lon_max':-40,'lat_min':15,'lat_max':23},
            'JOSE': {'lon_min':-70,'lon_max':-30,'lat_min':9,'lat_max':20},
            'MARIA':{'lon_min':-75,'lon_max':-30,'lat_min':9,'lat_max':20}}
    hour_step = 6 # 6 hours
    #------------------------------------

    # Read best track
    d_btk = {}
    for ist in Storms:
        d_btk[ist] = UD.btk_in_duration(small_dir,ist,start_time_str[ist],end_time_str[ist],hour_step)

    # Plot
    plot_btk()


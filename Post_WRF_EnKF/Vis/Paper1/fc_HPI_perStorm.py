import os,fnmatch # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
from netCDF4 import Dataset
from wrf import getvar
import math
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
import copy

import Util_data as UD
import Track

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['lines.linewidth'] = 2.5#1.5
matplotlib.rcParams['lines.markersize'] = 2.5
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['font.size'] = 8

# Initialized forecast times
def fc_iniT( Storm ):
    if Storm == 'HARVEY':
        fc_init = ['201708221800','201708230000','201708230600','201708231200']
    elif Storm == 'IRMA':
        fc_init = ['201709030600','201709031200','201709031800','201709040000']
    elif Storm == 'JOSE':
        fc_init = ['201709050600','201709051200','201709051800','201709060000']
    elif Storm == 'MARIA':
        fc_init = ['201709160600','201709161200','201709161800','201709170000']
    else:
        raise ValueError('Storm does not exist!')
    return fc_init

# Best track: start and end of the period
def btk_time( Storm ):
    if Storm == 'HARVEY':
        Btk_start = '201708230000' # '201709161800' #'201709030600'
        Btk_end = '201708270000' # '201709210000' #'201709090000'
    elif Storm == 'IRMA':
        Btk_start = '201709030000'
        Btk_end = '201709080000' #201709080000
    elif Storm == 'MARIA':
        Btk_start = '2017091600'#'201709160000'
        Btk_end = '201709210000'
    elif Storm == 'JOSE':
        Btk_start = '201709050000'
        Btk_end = '201709100000'
    else:
        raise ValueError('Storm does not exist!')
    return Btk_start, Btk_end


def plot_one_hpi( ax0, ax1, ax2,  state, color, line, line_width, label, steps=6, btk=None):

    # Best track
    if label == 'Best track':
        # Read data
        times = state['time']
        lon = state['lon']
        lat = state['lat']
        x_mslp = state['mslp']
        x_vmax = state['vmax']

        # track
        ax0.plot(lon, lat, marker='o',markersize=3, color=color,linewidth=3, label=label, linestyle=line, transform=ccrs.PlateCarree())
        dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]
        for it in times:
            if it[8:10] == '00':
                boolean_compare = [ it  == eachT for eachT in times ]
                idx = int( np.where(boolean_compare)[0] )
                ax0.scatter( lon[idx], lat[idx],s=5, marker='o',edgecolor="white",transform=ccrs.PlateCarree())
                if lon[idx]>lon_max or lon[idx]<lon_min or lat[idx]<lat_min or lat[idx]>lat_max:
                    continue
                elif Storm == 'JOSE' and it == '201709100000':
                    continue
                elif Storm == 'MARIA' and it == '201709210000':
                    continue
                else:
                    ax0.annotate(it[6:8], xy=(lon[idx], lat[idx]), color=color, xycoords='data', transform=ccrs.PlateCarree())
        # intensity
        ax1.plot_date(dates, x_mslp, color, label=label, linewidth=3)
        ax2.plot_date(dates, x_vmax, color, label=label, linewidth=3)

    # Model simulation
    else:
       # Read data (and process it if it is special)
        Diff_start_next6x = state['Diff_start_next6x']
        if Diff_start_next6x != 6 and Diff_start_next6x != 0: # it means time_to_move in namelist is not equal to 0
            # 
            start_next6x = datetime.strptime(state['time'][0],"%Y%m%d%H%M") + timedelta(hours=Diff_start_next6x)
            boolean_compare = [ start_next6x.strftime("%Y%m%d%H%M") == time_str for time_str in state['time'][:] ]
            idx_next6x = int( np.where(boolean_compare)[0] )

            times = state['time'][idx_next6x::steps]
            lon = state['lon'][idx_next6x::steps]
            lat = state['lat'][idx_next6x::steps]
            x_mslp = state['mslp'][idx_next6x::steps]
            x_vmax = state['vmax'][idx_next6x::steps]
            # Approximate the HPI at initialization time with the 1st record in ATCF_rsl.error.0000
            times = np.insert(times,0,state['time'][0])
            lon = np.insert(lon,0,state['lon'][0])
            lat = np.insert(lat,0,state['lat'][0])
            x_mslp = np.insert(x_mslp,0,state['mslp'][0])
            x_vmax = np.insert(x_vmax,0,state['vmax'][0])
        else:
            times = state['time'][::steps]
            lon = state['lon'][::steps]
            lat = state['lat'][::steps]
            x_mslp = state['mslp'][::steps]
            x_vmax = state['vmax'][::steps]

        # track        
        ax0.plot(lon, lat, marker='o', markersize=3, color=color,linewidth=line_width, label=label, linestyle=line, transform=ccrs.PlateCarree())
        for it in times:
            if it[8:10] == '00':
                boolean_compare = [ it  == eachT for eachT in times ]
                idx = int( np.where(boolean_compare)[0] )
                ax0.scatter( lon[idx], lat[idx],s=5, marker='o',edgecolor="white",transform=ccrs.PlateCarree())
                if lon[idx]>lon_max or lon[idx]<lon_min or lat[idx]<lat_min or lat[idx]>lat_max:
                    continue
                elif Storm == 'JOSE' and it == '201709100000':
                    continue
                elif Storm == 'MARIA' and it == '201709210000':
                    continue
                else:
                    ax0.annotate(it[6:8], xy=(lon[idx], lat[idx]), color=color, xycoords='data', transform=ccrs.PlateCarree())
                #ax0.text((state['lon'][idx], state['lat'][idx],it[6:8], fontsize=4,transform=ccrs.PlateCarree())

        # intensity
        dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]
        # MSLP
        ax1.plot_date(dates, x_mslp, color=color,  linestyle=line, linewidth=line_width,)
        # Vmax
        ax2.plot_date(dates, x_vmax, color=color, linestyle=line, linewidth=line_width,)


# Plot deterministic forecasts
def plot_hpi_df( ):

    num_exp = 3

    # Set the end point of the period to investigate
    if Storm == 'HARVEY':
        DF_model_end  = '201708270000'
    elif Storm == 'IRMA':
        DF_model_end  = '201709080000'
    elif Storm == 'MARIA':
        DF_model_end  = '201709210000'
    elif Storm == 'JOSE':
        DF_model_end  = '201709100000'
    else:
        DF_model_end  = None

    # Customize color maps
    if distinct_colors:
        Color1 = ["#c23728","#e14b31","#de6e56","#e1a692","#786028","#a57c1b","#d2980d","#ffb400","#503f3f","#6d4b4b","#a86464","#e27c7c"] #redish
        Color2 = ["#115f9a", "#1984c5", "#22a7f0", "#48b5c4", "#48446e", "#5e569b", "#776bcd", "#9080ff","#3c4e4b", "#466964", "#599e94", "#6cd4c5"] #blueish
    else:
        red_cm = cm.Reds
        blue_cm = cm.Blues
        num_colors = 14
        discretize_red = ListedColormap(red_cm(np.linspace(0,1,num_colors)))
        discretize_blue = ListedColormap(blue_cm(np.linspace(0,1,num_colors)))
        Color1 = discretize_red.colors[5::2] #discretize_red.colors[5:] 
        Color2 = discretize_blue.colors[5::2]#discretize_blue.colors[5:]
    Color_set = {'c0':Color1, 'c1':Color2}
    print( Color1 )
    Ana_color = ['#eea990','#748b97'] #'#748b97'
    Line_types = ['-','-']

    # Customize labels ###### Chnage it every time !!!!!!!!!!!!!!! 
    Labels = ['WSM6:','THO:']
    #Ana_labels = ['Ref Ans','MD Ans' ]

    # Best-track data
    Btk_start, Btk_end = btk_time( Storm )
    best_track = UD.btk_in_duration(small_dir, Storm, Btk_start, Btk_end, hour_step=6)

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200)
    widths = [2.78,1.79,1.79]
    #heights = [2.8,2.8,2.8]
    outer_grid = fig.add_gridspec(ncols=1,nrows=num_exp,bottom=0.03,top=0.96,left=0.10,right=0.98,wspace=0.05,hspace=0.03)

    ax = {}
    for ida in DA:
        ax[ida] = {}
        ir = DA.index( ida )
        # gridspec inside gridspec
        inner_grid = outer_grid[ir].subgridspec(1, 3, width_ratios=widths)
        ax[ida]['ax0'] = fig.add_subplot( inner_grid[0], projection=ccrs.PlateCarree())
        ax[ida]['ax0'].set_extent( domain_range,crs=ccrs.PlateCarree())
        ax[ida]['ax0'].coastlines( resolution='10m', color='black',linewidth=0.5 )
        if Storm =='HARVEY':
            ax[ida]['ax0'].set_aspect(1.2)
        elif Storm == 'IRMA':
            ax[ida]['ax0'].set_aspect(3.5)
        elif Storm == 'MARIA':
            ax[ida]['ax0'].set_aspect(3)
        elif Storm == 'JOSE':
            ax[ida]['ax0'].set_aspect(3.5)
        ax[ida]['ax1'] = fig.add_subplot( inner_grid[1] )
        ax[ida]['ax2'] = fig.add_subplot( inner_grid[2] )
        # add Saffir-Simpson Scale
        times = best_track['time']
        dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]
        ax[ida]['ax2'].fill_between([dates[0], dates[-1]], 0, 80, color='#71797E',alpha=0.5) # steel gray
        ax[ida]['ax2'].text(dates[1], 74, 'CAT 5', fontsize=6, fontweight='bold', color='white')
        ax[ida]['ax2'].fill_between([dates[0], dates[-1]], 0, 70, color='#6082B6',alpha=0.5) # Glaucous
        ax[ida]['ax2'].text(dates[1], 64, 'CAT 4', fontsize=6, fontweight='bold', color='white')
        ax[ida]['ax2'].fill_between([dates[0], dates[-1]], 0, 58, color='#8A9A5B',alpha=0.5) # sage green
        ax[ida]['ax2'].text(dates[1], 53.5, 'CAT 3', fontsize=6, fontweight='bold', color='white')
        ax[ida]['ax2'].fill_between([dates[0], dates[-1]], 0, 49, color='#A9A9A9',alpha=0.5) # dark gray
        ax[ida]['ax2'].text(dates[1], 45.5, 'CAT 2', fontsize=6, fontweight='bold', color='white')
        ax[ida]['ax2'].fill_between([dates[0], dates[-1]], 0, 42, color='#7393B3',alpha=0.5) # blue gray
        ax[ida]['ax2'].text(dates[1], 37, 'CAT 1', fontsize=6, fontweight='bold', color='white')
        ax[ida]['ax2'].fill_between([dates[0], dates[-1]], 0, 32, color='#B2BEB5',alpha=0.5) # ash gray
        ax[ida]['ax2'].text(dates[1], 24.5, 'TS', fontsize=6, fontweight='bold', color='white')
        ax[ida]['ax2'].fill_between([dates[0], dates[-1]], 0, 17, color='#E5E4E2',alpha=0.5) # Platinum
        ax[ida]['ax2'].text(dates[1], 11.5, 'TD', fontsize=6, fontweight='bold', color='white')
        # plot best track
        plot_one_hpi( ax[ida]['ax0'], ax[ida]['ax1'], ax[ida]['ax2'], best_track,  'black', '-', 4, 'Best track')

    # Plot simulaitons
    fc_init = fc_iniT( Storm )
    for imp in MP:
        iCtg = MP.index(imp) #category
        for ida in DA:
            for it in range(len(fc_init)):
                #print( big_dir+'/'+Storm+'/'+Exper_names[imp][ida]+'/wrf_df/'+fc_init[it] )
                if not os.path.exists( big_dir+'/'+Storm+'/'+Exper_names[imp][ida]+'/wrf_df/'+fc_init[it] ):
                    print('Forecast does not exist!')
                    continue
                # plot
                HPI_model = Track.read_rsl_error(Storm,Exper_names[imp][ida],big_dir+'/'+Storm+'/'+Exper_names[imp][ida]+'/wrf_df/'+fc_init[it],fc_init[it], DF_model_end) 
                plot_one_hpi( ax[ida]['ax0'], ax[ida]['ax1'], ax[ida]['ax2'], HPI_model, Color_set['c'+str(iCtg)][it], Line_types[iCtg], 1.5, Labels[iCtg]+fc_init[it], steps=6 )

    # Manullay control Legend
    if Storm == 'HARVEY':
        Loc = 'lower left'
    else:
        Loc = 'upper right'

    fc_times = copy.deepcopy( fc_init )
    lgd_1 = [MP[0]+':'+it for it in fc_times]
    lgd_1 = ['Best Track'] + lgd_1
    lines_DA0 = ax[DA[1]]['ax0'].get_lines()
    legend1 = ax[DA[0]]['ax0'].legend([lines_DA0[i] for i in [0,1,2,3,4]], lgd_1,fontsize='7',loc=Loc)

    lgd_2 = [MP[1]+':'+it for it in fc_times]
    lines_DA1 = ax[DA[1]]['ax0'].get_lines() 
    #if Storm == 'HARVEY':
    #    pass
    #    legend2 = ax[DA[1]]['ax1'].legend([lines_DA1[i] for i in [5,6,7,8]], lgd_2,fontsize='7',loc=Loc)
    #else:
    legend2 = ax[DA[1]]['ax0'].legend([lines_DA1[i] for i in [5,6,7,8]], lgd_2,fontsize='7',loc=Loc)

    # Set Minor ticks/labels for track subplot
    lon_ticks =  np.arange(lon_min, lon_max+1, 1)
    lat_ticks = np.arange(lat_min, lat_max+1, 1)
    for ida in DA:
        gl = ax[ida]['ax0'].gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=0.8, color='gray', alpha=0.3, linestyle='--')
        gl.left_labels = False
        gl.bottom_labels = False
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.ylocator = mticker.FixedLocator(lat_ticks)

    # Set Major ticks/labels for track subplot
    lon_ticks =  np.arange(lon_min, lon_max+5, 5)
    lat_ticks = np.arange(lat_min, lat_max+5, 5)
    #lon_ticks = list(range(math.ceil(domain_range[0])-2, math.ceil(domain_range[1])+3, 4))
    #lat_ticks = list(range(math.ceil(domain_range[2])-2, math.ceil(domain_range[3])+2, 2))
    for ida in DA:
        gl = ax[ida]['ax0'].gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=1, color='gray', alpha=0.5, linestyle='-')
        gl.left_labels = True
        if DA.index( ida ) == num_exp-1:
            gl.bottom_labels = True
        else:
            gl.bottom_labels = False
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.yformatter = LATITUDE_FORMATTER
        gl.xformatter = LONGITUDE_FORMATTER
        gl.xlabel_style = {'size': 8}
        gl.ylabel_style = {'size': 8}


    # Set ticks/other attributes for intensity subplots
    for ida in DA:
        date_form = mdates.DateFormatter("%m-%d")
        ax[ida]['ax1'].set_xlim([datetime(int(Btk_start[0:4]),int(Btk_start[4:6]), int(Btk_start[6:8]), int(Btk_start[8:10])), datetime(int(Btk_start[0:4]),int(Btk_end[4:6]), int(Btk_end[6:8]), int(Btk_end[8:10]))])
        ax[ida]['ax2'].set_xlim([datetime(int(Btk_start[0:4]),int(Btk_start[4:6]), int(Btk_start[6:8]), int(Btk_start[8:10])), datetime(int(Btk_start[0:4]),int(Btk_end[4:6]), int(Btk_end[6:8]), int(Btk_end[8:10]))])
        ax[ida]['ax1'].xaxis.set_major_locator(mdates.DayLocator())
        ax[ida]['ax2'].xaxis.set_major_locator(mdates.DayLocator())
        ax[ida]['ax1'].xaxis.set_major_formatter(date_form)
        ax[ida]['ax2'].xaxis.set_major_formatter(date_form)
        if DA.index( ida ) == num_exp-1:
            ax[ida]['ax1'].tick_params(axis='x', labelrotation=20,labelsize=8)
            ax[ida]['ax1'].tick_params(axis='x', labelrotation=20,labelsize=8)
            ax[ida]['ax2'].tick_params(axis='x', labelrotation=20,labelsize=8)
        else:
            ax[ida]['ax1'].set_xticklabels([])  # Set x ticks as empty
            ax[ida]['ax2'].set_xticklabels([])  # Set x ticks as empty
        ax[ida]['ax1'].set_ylim([900,1020])     #([940, 1015])
        ax[ida]['ax2'].set_ylim([10,80])   #([10,60])
        ax[ida]['ax1'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')
        ax[ida]['ax2'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    # Set y label
    for ida in DA:
        if DA.index(ida) == 0:
            fig.text(0.02,0.82,ida, fontsize=12, ha='center', va='center',rotation='vertical')
        elif DA.index(ida) == 1:
            fig.text(0.02,0.50,ida, fontsize=12, ha='center', va='center',rotation='vertical')
        elif DA.index(ida) == 2:
            fig.text(0.02,0.18,ida, fontsize=12, ha='center', va='center',rotation='vertical')

    # Set panel labels
    #for ida in DA:
    #    if DA.index(ida) == 0:
    #        fig.text(0.41,0.67,'(a1)', fontsize=12, ha='center', va='center')
    #        fig.text(0.68,0.67,'(a2)', fontsize=12, ha='center', va='center')
    #        fig.text(0.95,0.67,'(a3)', fontsize=12, ha='center', va='center')
    #    elif DA.index(ida) == 1:
    #        fig.text(0.41,0.355,'(b1)', fontsize=12, ha='center', va='center')
    #        fig.text(0.68,0.355,'(b2)', fontsize=12, ha='center', va='center')
    #        fig.text(0.95,0.355,'(b3)', fontsize=12, ha='center', va='center')
    #    elif DA.index(ida) == 2:
    #        fig.text(0.41,0.044,'(c1)', fontsize=12, ha='center', va='center')
    #        fig.text(0.68,0.044,'(c2)', fontsize=12, ha='center', va='center')
    #        fig.text(0.95,0.044,'(c3)', fontsize=12, ha='center', va='center')

    # Set titles
    ax[DA[0]]['ax0'].set_title( 'Track',fontsize = 12 )
    ax[DA[0]]['ax1'].set_title( 'MSLP (hPa)',fontsize = 12 )
    ax[DA[0]]['ax2'].set_title( 'Vmax (m $\mathregular{s^{-1}}$)',fontsize = 12 )

    # Save figure
    des_name = small_dir+Storm+'/Paper1/'+Storm+'_fc_HPI.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )




if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/Clean_results/'

    # Configuration
    Storm = 'IRMA'
    MP = ['WSM6','THO']
    DA = ['CONV','IR','MW']

    distinct_colors = False
    Plot_analyses = False # Feature that plots the analyses of an experiment
    Plot_HPI = True
    Plot_abs_error = False
    #fc_run_hrs = 60

    # Generate experiment names
    Exper_names = {}
    for imp in MP:
        Exper_names[imp] = {}
        for ida in DA:
            Exper_names[imp][ida] = UD.generate_one_name( Storm,ida,imp )

    # Pre-set the domain for DF forecast
    if Storm == 'HARVEY':
        lon_min = -102
        lon_max = -87
        lat_min = 16
        lat_max = 31
    elif Storm == 'IRMA':
        lon_min = -73
        lon_max = -43
        lat_min = 14
        lat_max = 24
    elif Storm == 'MARIA':
        lon_min = -68
        lon_max = -43#-45
        lat_min = 9.3# 
        lat_max = 19.3
    elif Storm == 'JOSE':
        lon_min = -65
        lon_max = -35#-45
        lat_min = 9.5
        lat_max = 19.5
    domain_range = [lon_min-0.01, lon_max+0.01, lat_min-0.01, lat_max+0.01]

    if Plot_HPI:
        plot_hpi_df( )

    if Plot_abs_error:
        start_time=time.process_time()
        plot_abs_err( Config )
        end_time = time.process_time()
        print('time needed: ', end_time-start_time, ' seconds')



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

# Specify time range for model data
def t_end_model( Storm ):

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

    return DF_model_end

# Best track: start and end of the period
def t_range_btk( Storm ):
    if Storm == 'HARVEY':
        Btk_start = '201708221800' 
        Btk_end = '201708270000' 
    elif Storm == 'IRMA':
        Btk_start = '201709030600'
        Btk_end = '201709080000' 
    elif Storm == 'MARIA':
        Btk_start = '201709160600'
        Btk_end = '201709210000'
    elif Storm == 'JOSE':
        Btk_start = '201709050600'
        Btk_end = '201709100000'
    else:
        raise ValueError('Storm does not exist!')
    return Btk_start, Btk_end

# Add one forecast for one var
def add_one_data( ax0,state,line_color,line_type,line_width,label,steps=6 ):

    # --- Extract varibales 
    # Best track
    if label == 'btk':
        times = state['time']
        x_mslp = state['mslp']
    else:
       # Read data (and process it if it is special)
        Diff_start_next6x = state['Diff_start_next6x']
        if Diff_start_next6x != 6 and Diff_start_next6x != 0: # it means time_to_move in namelist is not equal to 0
            # 
            start_next6x = datetime.strptime(state['time'][0],"%Y%m%d%H%M") + timedelta(hours=Diff_start_next6x)
            boolean_compare = [ start_next6x.strftime("%Y%m%d%H%M") == time_str for time_str in state['time'][:] ]
            idx_next6x = int( np.where(boolean_compare)[0] )

            times = state['time'][idx_next6x::steps]
            x_mslp = state['mslp'][idx_next6x::steps]
            # Approximate the HPI at initialization time with the 1st record in ATCF_rsl.error.0000
            times = np.insert(times,0,state['time'][0])
            x_mslp = np.insert(x_mslp,0,state['mslp'][0])
        else:
            times = state['time'][::steps]
            x_mslp = state['mslp'][::steps]

    # --- Plot
    dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]
    ax0.plot_date(dates, x_mslp, color=line_color, linestyle=line_type, label=label, linewidth=line_width)


def plot_fc_mslp( ):

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=300) # standard: 6.5,8.5
    grids = fig.add_gridspec(ncols=4,nrows=3,hspace=0.04,wspace=0.03,top=0.93,bottom=0.10,left=0.1,right=0.95)
    ax = {}
    for ida in DA:
        ax[ida] = {}
        for ist in Storms:
            ax[ida][ist] = fig.add_subplot( grids[DA.index(ida),Storms.index(ist)] )

     # Customize color maps
    red_cm = cm.Reds
    blue_cm = cm.Blues
    num_colors = 14
    discretize_red = ListedColormap(red_cm(np.linspace(0,1,num_colors)))
    discretize_blue = ListedColormap(blue_cm(np.linspace(0,1,num_colors)))
    Color1 = discretize_red.colors[3::2] #discretize_red.colors[5:] 
    Color2 = discretize_blue.colors[3::2]#discretize_blue.colors[5:]
    Color_set = {'WSM6':Color1, 'THO':Color2}

    Line_types = ['-','-']

    # Plot
    for ida in DA: # loop thru rows
        for ist in Storms: # loop thru storm names
            fc_init = fc_iniT( ist )
            # plot best track
            add_one_data( ax[ida][ist],btk[ist],'black','-',3,'btk',steps=6 )
            # plot forecasts
            for imp in MP:
                iCtg = MP.index(imp) #category
                for it in range(len(fc_init)):
                    add_one_data( ax[ida][ist],fc[ida][ist][imp][fc_init[it]],Color_set[imp][it], Line_types[iCtg], 1.5, imp, steps=6)   

    # Legend
    for ist in Storms:

        # specify location
        if ist == 'IRMA':
            Loc = 'upper right'
        else:
            Loc = 'lower left'

        fc_init = fc_iniT( ist )

        lgd_0 = ['Best Track',] 
        lines_0 = ax[DA[0]][ist].get_lines()
        legend0 = ax[DA[0]][ist].legend([lines_0[i] for i in [0,]], lgd_0,fontsize='7',loc=Loc)

        lgd_1 = [MP[0]+':'+it for it in fc_init]
        lines_1 = ax[DA[1]][ist].get_lines()
        legend1 = ax[DA[1]][ist].legend([lines_1[i] for i in [1,2,3,4]], lgd_1,fontsize='6',loc=Loc)

        lgd_2 = [MP[1]+':'+it for it in fc_init]
        lines_2 = ax[DA[2]][ist].get_lines()
        legend2 = ax[DA[2]][ist].legend([lines_2[i] for i in [5,6,7,8]], lgd_2,fontsize='6',loc=Loc)


    # Set ticks/other attributes for intensity subplots
    for ist in Storms:
        Btk_start, Btk_end = t_range_btk( ist )
        for ida in DA:
            # x axis
            date_form = mdates.DateFormatter("%m-%d")
            ax[ida][ist].set_xlim([datetime(int(Btk_start[0:4]),int(Btk_start[4:6]), int(Btk_start[6:8]), int(Btk_start[8:10])), datetime(int(Btk_start[0:4]),int(Btk_end[4:6]), int(Btk_end[6:8]), int(Btk_end[8:10]))])
            ax[ida][ist].xaxis.set_major_locator(mdates.DayLocator())
            ax[ida][ist].xaxis.set_major_formatter(date_form)
            if DA.index( ida ) == len(DA)-1:
                ax[ida][ist].tick_params(axis='x', labelrotation=20,labelsize=8)
                ax[ida][ist].tick_params(axis='x', labelrotation=20,labelsize=8)
            else:
                ax[ida][ist].set_xticklabels([])  # Set x ticks as empty
            # y axis
            ax[ida][ist].set_ylim([900,1020])     #([940, 1015])
            if Storms.index( ist ) != 0:
                ax[ida][ist].set_yticklabels([])
            ax[ida][ist].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

   # Set y label
    for ida in DA:
        if DA.index(ida) == 0:
            fig.text(0.03,0.80,ida+': MSLP (hPa)', fontsize=12, ha='center', va='center',rotation='vertical')
        elif DA.index(ida) == 1:
            fig.text(0.03,0.52,ida+': MSLP (hPa)', fontsize=12, ha='center', va='center',rotation='vertical')
        elif DA.index(ida) == 2:
            fig.text(0.03,0.24,ida+': MSLP (hPa)', fontsize=12, ha='center', va='center',rotation='vertical')

    # Set titles
    for ist in Storms:
        ax[DA[0]][ist].set_title( ist,fontsize = 12 )


    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/fc_mslp_allStorms.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/Clean_results/'

    #--------Configuration------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    DA = ['CONV','IR','MW']
    MP = ['WSM6','THO'] #

    var = 'mslp' # track, vmax
    #------------------------------------
    wrf_dir = big_dir

    if_btk = True

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = {}
        for imp in MP:
            Exper_names[istorm][imp] = {}
            for ida in DA:
                Exper_names[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # Read best-track
    if if_btk:
        btk = {}
        for ist in Storms:
            Btk_start, Btk_end = t_range_btk( ist )
            btk[ist] = UD.btk_in_duration(small_dir, ist, Btk_start, Btk_end, hour_step=6)

    # Read forecasts
    fc = {}
    for ida in DA:
        fc[ida] = {}
        for ist in Storms:
            fc_init = fc_iniT( ist )
            fc[ida][ist] = {}
            for imp in MP:
                fc[ida][ist][imp] = {}
                for it in range(len(fc_init)):
                    fc_dir = big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida]+'/wrf_df/'+fc_init[it]
                    if not os.path.exists( fc_dir ):
                        print('Forecast does not exist!')
                        continue
                    # read forecast launched at a time
                    fc[ida][ist][imp][fc_init[it]] = Track.read_rsl_error(ist,Exper_names[ist][imp][ida],fc_dir,fc_init[it], t_end_model(ist) )

    
    # Plot
    if var == 'mslp':
        plot_fc_mslp( )



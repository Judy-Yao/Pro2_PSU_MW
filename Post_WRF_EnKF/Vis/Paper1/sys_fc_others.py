import os,fnmatch # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
#from wrf import getvar
import math
import scipy as sp
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
import numpy.ma as ma
from scipy import interpolate
import copy

import EnKF_Vmax as Vx #Vmax
import EnKF_minSlp_track as SC #StormCenter
import Util_data as UD

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['lines.linewidth'] = 1.5#1.5
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
def t_range_model( Storm ):

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

# plot: start and end of the plotting period$
def plot_time( Storm ):
    if Storm == 'HARVEY':
        plot_start = '201708221800' # '201709161800' #'201709030600'$
        plot_end = '201708270000' # '201709210000' #'201709090000'$
    elif Storm == 'IRMA':
        plot_start = '201709030000'
        plot_end = '201709080000' #201709080000$
    elif Storm == 'MARIA':
        plot_start = '2017091600'#'201709160000'$
        plot_end = '201709210000'
    elif Storm == 'JOSE':
        plot_start = '201709050000'
        plot_end = '201709100000'
    else:
        raise ValueError('Storm does not exist!')
    return plot_start, plot_end

# Customize color maps
def cst_color():
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
    return Color_set

# Read fields for a forecast
def read_forecast(Storm,exp_dir):

    times = []
    d_fc = {var:[] for var in varnames} 

    # list hourly wrfout
    wrfouts = sorted(glob.glob( exp_dir+'/wrfout_d03_*00' ) )

    for wrfout in wrfouts:
        # get the directory and file name
        dire, filename = os.path.split( wrfout )
        times.append( filename[11:26].replace('-', '').replace('_', '').replace(':', '') )

        for var in varnames:
            if var == 'Precip':
                ncdir = nc.Dataset( wrfout )
                ppb = ncdir.variables['RAINNC'][0,:,:]
                # Make precipitation related to convection stand out
                #mask = ppb <= 0.5
                #ppb_masked = ma.masked_array(ppb, mask=mask)
                #x_pp.append( np.ma.mean(ppb) )

            d_fc[var].append( (np.sum(ppb)/1e6) )  # /1e6
    
    d_fc[var] = np.array( d_fc[var] )

    # add times
    d_fc['time'] = np.array( times )

    return d_fc


def plot_one( ax, var, state, color, line, line_width, label, steps=1 ):

    times = state['time'][::steps]  
    dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]
    ax.plot(dates,state[var],color=color,linestyle=line,linewidth=line_width) # instead of plot_date


# simply test 
def plot_single_storm( Storm, var):

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200)
    grids = fig.add_gridspec(ncols=1,nrows=3,bottom=0.03,top=0.96,left=0.10,right=0.98,wspace=0.05,hspace=0.03)

    ax = {}
    for ida in DA:
        ax[ida] = {}
        ir = DA.index( ida )
        ax[ida] = fig.add_subplot( grids[ir] )

    # Customize color 
    Color_set = cst_color() 
    
    Line_types = ['-','-']

    # Customize labels ###### Chnage it every time !!!!!!!!!!!!!!! 
    Labels = ['WSM6:','THO:']

    # Plot simulations
    fc_init = fc_iniT( Storm )
    fc_end = t_range_model( Storm )

    for imp in MP:
        iCtg = MP.index(imp) #category
        for ida in DA:
            for it in range(len(fc_init)):
                print(Storm+': '+ida+'_'+imp+' '+fc_init[it])
                if not os.path.exists( big_dir+'/'+Storm+'/'+Exper_names[Storm][imp][ida]+'/wrf_df/'+fc_init[it] ):
                    print('Forecast does not exist!')
                    continue
                # plot
                d_fc = read_forecast(Storm,big_dir+'/'+Storm+'/'+Exper_names[Storm][imp][ida]+'/wrf_df/'+fc_init[it]) #,fc_init[it], DF_model_end)
                plot_one( ax[ida], var, d_fc, Color_set['c'+str(iCtg)][it], Line_types[iCtg], 1.5, Labels[iCtg]+fc_init[it], steps=1 )

    # Set ticks/other attributes for intensity subplots
    for ida in DA:
        date_form = mdates.DateFormatter("%m-%d")
        ax[ida].set_xlim([datetime(int(fc_init[0][0:4]),int(fc_init[0][4:6]), int(fc_init[0][6:8]), int(fc_init[0][8:10])), datetime(int(fc_end[0:4]),int(fc_end[4:6]), int(fc_end[6:8]), int(fc_end[8:10]))])
        ax[ida].xaxis.set_major_locator(mdates.DayLocator())
        ax[ida].xaxis.set_major_formatter(date_form)
        ax[ida].tick_params(axis='x', labelrotation=20,labelsize=8)
        #ax[ida].set_ylim([900,1020])     #([940, 1015])
        ax[ida].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')


    # Save figure
    des_name = small_dir+Storm+'/Paper1/'+Storm+'_fc_precip.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


def plot_4by3(  ):


    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=4,top=0.93,left=0.12,hspace=0.11)

    ax = {}
    for ist in Storms:
        ax[ist] = {}
        ir = Storms.index(ist)
        ist_grids = outer_grids[ir].subgridspec(1, 3, wspace=0.05)
        for ida in DA:
            ax[ist][ida] = fig.add_subplot( ist_grids[DA.index(ida)] )

    # Customization 
    Color_set = cst_color()
    Line_types = ['-','-']
    Labels = ['WSM6:','THO:']
    # Plot simulations
    for ist in Storms:
        # loop thru each storm
        fc_init = fc_iniT( ist )
        fc_end = t_range_model( ist ) 
        # loop thru each DA
        for ida in DA:
            for imp in MP:
                iCtg = MP.index(imp) #category
                for it in range(len(fc_init)):
                    print(ist+': '+ida+'_'+imp+' '+fc_init[it])
                    if not os.path.exists( big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida]+'/wrf_df/'+fc_init[it] ):
                        print('Forecast does not exist!')
                        continue
                    # plot
                    d_fc = read_forecast(ist,big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida]+'/wrf_df/'+fc_init[it]) #,fc_init[it], DF_model_end)
                    plot_one( ax[ist][ida], var, d_fc, Color_set['c'+str(iCtg)][it], Line_types[iCtg], 1.5, Labels[iCtg]+fc_init[it], steps=1 )


    # Manullay control Legend
    Loc = 'upper left'

    for ist in Storms:
        fc_init = fc_iniT( ist )
        fc_times = copy.deepcopy( fc_init )
        lgd_1 = [MP[0]+':'+it for it in fc_times]
        lines_DA0 = ax[ist]['IR'].get_lines()
        legend1 = ax[ist]['IR'].legend([lines_DA0[i] for i in [0,]], lgd_1,fontsize='7',loc=Loc)
        legend1.set_alpha( 0.5 )

        lgd_2 = [MP[1]+':'+it for it in fc_times]
        lines_DA1 = ax[ist]['IR+MW'].get_lines() 
        legend2 = ax[ist]['IR+MW'].legend([lines_DA1[i] for i in [1,]], lgd_2,fontsize='7',loc=Loc)
        legend2.set_alpha( 0.5 )

    # Set ticks/other attributes 
    date_form = mdates.DateFormatter("%m-%d")
    for ist in Storms:
        # loop thru each stormplot_time( Storm )
        plot_st, plot_end = plot_time( ist )
        # loop thru each DA
        for ida in DA:
            # x axis
            ax[ist][ida].set_xlim([datetime(int(plot_st[0:4]),int(plot_st[4:6]), int(plot_st[6:8]), int(plot_st[8:10])), datetime(int(plot_end[0:4]),int(plot_end[4:6]), int(plot_end[6:8]), int(plot_end[8:10]))])
            # Customize x-axis tick labels precisely
            # Define custom tick locations
            if DA.index( ida ) == 1:
                if ist == 'HARVEY':
                    plot_t = ['201708240000','201708250000','201708260000']
                elif ist == 'JOSE':
                    plot_t = ['201709060000','201709070000','201709080000','201709090000']
                elif ist == 'IRMA':
                    plot_t = ['201709040000','201709050000','201709060000','201709070000']
                elif ist == 'MARIA':
                    plot_t = ['201709170000','201709180000','201709190000','201709200000']
                else:
                    pass
                tick_locations = [datetime.strptime(i,"%Y%m%d%H%M") for i in plot_t]
                # Define custom tick labelsi
                tick_labels = [date.strftime('%m-%d') for date in tick_locations]
                # Set x labels
                ax[ist][ida].set_xticks(tick_locations)
                ax[ist][ida].set_xticklabels(tick_labels) 
            else:
                ax[ist][ida].xaxis.set_major_locator(mdates.DayLocator())
                ax[ist][ida].xaxis.set_major_formatter(date_form)
            ax[ist][ida].tick_params(axis='x',labelsize=8)
            # y axis
            ax[ist][ida].set_ylim([0,7]) # X*1e6
            ax[ist][ida].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')
            if DA.index( ida ) != 0:
                ax[ist][ida].set_yticklabels([]) 


    # Add DA name
    for ida in DA:
        if DA.index(ida) == 0:
            fig.text(0.25,0.95,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 1:
            fig.text(0.51,0.95,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 2:
            fig.text(0.79,0.95,ida, fontsize=12, ha='center', va='center')

    # Add y label
    fig.text(0.03,0.51,'Area total grid scale precipitation ($10^6$ mm) ', fontsize=12, ha='center', va='center',rotation='vertical') 
    fig.text(0.08,0.85,'HARVEY', fontsize=10, ha='center', va='center',rotation='vertical') 
    fig.text(0.08,0.63,'JOSE', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.08,0.42,'IRMA', fontsize=10, ha='center', va='center' ,rotation='vertical')
    fig.text(0.08,0.21,'MARIA', fontsize=10, ha='center', va='center',rotation='vertical')


    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_precip.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


def plot_4by2(  ):

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=4,top=0.93,right=0.97,hspace=0.11)

    ax = {}
    for ist in Storms:
        ax[ist] = {}
        ir = Storms.index(ist)
        ist_grids = outer_grids[ir].subgridspec(1, 2, wspace=0.05)
        for imp in MP:
            ax[ist][imp] = fig.add_subplot( ist_grids[MP.index(imp)] )

    # Customization 
    Color_set = cst_color()
    Line_types = {'CONV':'-','IR':'--','IR+MW':(0, (1, 1))}

    # Plot simulations
    for ist in Storms:
        # loop thru each storm
        fc_init = fc_iniT( ist )
        fc_end = t_range_model( ist ) 
        # loop thru each MP
        for imp in MP:
            iCtg = MP.index(imp) #category
            for ida in DA[1:]:
                for it in range(len(fc_init)):
                    print(ist+': '+ida+'_'+imp+' '+fc_init[it])
                    if not os.path.exists( big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida]+'/wrf_df/'+fc_init[it] ):
                        print('Forecast does not exist!')
                        continue
                    # plot
                    d_fc = read_forecast(ist,big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida]+'/wrf_df/'+fc_init[it]) #,fc_init[it], DF_model_end)
                    plot_one( ax[ist][imp], var, d_fc, Color_set['c'+str(iCtg)][it], Line_types[ida], 2, fc_init[it], steps=1 )

    # Legends
    for ist in Storms:
        fc_init = fc_iniT( ist )
        fc_times = copy.deepcopy( fc_init )

        lines = ax[ist]['WSM6'].get_lines()
        lgd_0 = [it for it in fc_times]
        legend0 = ax[ist]['WSM6'].legend([lines[i] for i in [0,]], lgd_0,fontsize='7',loc='upper left')
        legend0.set_alpha( 0.5 )
        # Add the first legend manually to the current Axes
        ax[ist]['WSM6'].add_artist(legend0)
        
        #lgd_1 = ['CONV','IR','IR+MW']
        #legend1 = ax[ist]['WSM6'].legend([lines[i] for i in [0,1,2]], lgd_1,fontsize='7',loc='lower right',ncol=3)
        #legend1.set_alpha( 0.5 )

        lines = ax[ist]['THO'].get_lines()
        lgd_0 = [it for it in fc_times]
        legend0 = ax[ist]['THO'].legend([lines[i] for i in [0,]], lgd_0,fontsize='7',loc='upper left')
        legend0.set_alpha( 0.5 )

    # Create proxy artists for the second legend
    proxy_artist1 = plt.Line2D((0, 1), (0, 0), color=Color_set['c0'][0], linestyle='-')
    proxy_artist2 = plt.Line2D((0, 1), (0, 0), color=Color_set['c0'][0], linestyle='--')
    proxy_artist3 = plt.Line2D((0, 1), (0, 0), color=Color_set['c0'][0], linestyle=':')
    ax['HARVEY']['WSM6'].legend([proxy_artist1, proxy_artist2, proxy_artist3], ['CONV', 'IR', 'IR+MW'],fontsize='7',loc='lower right',ncol=3)

    # Set ticks/other attributes
    date_form = mdates.DateFormatter("%m-%d")
    for ist in Storms:
        # loop thru each stormplot_time( Storm )
        plot_st, plot_end = plot_time( ist )
        for imp in MP:
            # x axis
            ax[ist][imp].set_xlim([datetime(int(plot_st[0:4]),int(plot_st[4:6]), int(plot_st[6:8]), int(plot_st[8:10])), datetime(int(plot_end[0:4]),int(plot_end[4:6]), int(plot_end[6:8]), int(plot_end[8:10]))])
            if MP.index( imp ) == 0:
                if ist == 'HARVEY':
                    plot_t = ['201708230000','201708240000','201708250000','201708260000']
                elif ist == 'JOSE':
                    plot_t = ['201709050000','201709060000','201709070000','201709080000','201709090000']
                elif ist == 'IRMA':
                    plot_t = ['201709030000','201709040000','201709050000','201709060000','201709070000']
                elif ist == 'MARIA':
                    plot_t = ['201709160000','201709170000','201709180000','201709190000','201709200000']
                else:
                    pass
                tick_locations = [datetime.strptime(i,"%Y%m%d%H%M") for i in plot_t]
                # Define custom tick labelsi
                tick_labels = [date.strftime('%m-%d') for date in tick_locations]
                # Set x labels
                ax[ist][imp].set_xticks(tick_locations)
                ax[ist][imp].set_xticklabels(tick_labels) 
            else:
                # Customize x-axis tick labels precisely
                ax[ist][imp].xaxis.set_major_locator(mdates.DayLocator())
                ax[ist][imp].xaxis.set_major_formatter(date_form)
                ax[ist][imp].tick_params(axis='x',labelsize=8)
            # y axis
            ax[ist][imp].set_yscale('log')
            #ax[ist][imp].set_ylim([0,7]) # X*1e6
            #if MP.index(imp) != 0:
            #     ax[ist][imp].set_yticklabels([])
            #ax[ist][imp].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    # Add y label
    fig.text(0.03,0.51,'Area total grid scale precipitation ($10^6$ mm) ', fontsize=12, ha='center', va='center',rotation='vertical')
    fig.text(0.08,0.85,'HARVEY', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.08,0.63,'JOSE', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.08,0.42,'IRMA', fontsize=10, ha='center', va='center' ,rotation='vertical')
    fig.text(0.08,0.21,'MARIA', fontsize=10, ha='center', va='center',rotation='vertical')

    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_precip.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )



# ------------------------------------------------------------------------------------------------------
#            Operation: Plot Differences
# ------------------------------------------------------------------------------------------------------

# Add a smaller subplot in a subplot.
# See: https://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
def add_subplot_axes(ax,rect):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax


def plot_one_diff( ist, imp, conv_each, ax, var, state, color, line, line_width, label, steps=1 ):

    times = state['time'][::steps]
    dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]
    if ist == 'HARVEY':
        ax.plot(dates,state[var]-conv_each[var][::3],color=color,linestyle=line,linewidth=line_width) # instead of plot_date
    else:
        ax.plot(dates,state[var]-conv_each[var],color=color,linestyle=line,linewidth=line_width) # instead of plot_date

def plot_4by2_diff(  ):

    # Read conv data
    conv = {}
    for ist in Storms:
        fc_init = fc_iniT( ist )
        conv[ist] = {}
        for imp in MP:  
            conv[ist][imp] = {}
            for it in fc_init:
                print('Reading CONV experiment: '+imp+' '+it)
                conv[ist][imp][it] = read_forecast(ist,big_dir+'/'+ist+'/'+Exper_names[ist][imp]['CONV']+'/wrf_df/'+it) 

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=4,top=0.968,right=0.95,hspace=0.11)

    ax = {}
    for ist in Storms:
        ax[ist] = {}
        ir = Storms.index(ist)
        ist_grids = outer_grids[ir].subgridspec(1, 2, wspace=0.02)
        for imp in MP:
            ax[ist][imp] = fig.add_subplot( ist_grids[MP.index(imp)] )

    # Customization 
    Color_set = cst_color()
    Line_types = {'IR':'-','IR+MW':(0, (1, 1))}

    # Plot simulations
    subax = add_subplot_axes(ax['JOSE']['THO'],[0.1,0.05,0.35,0.25]) # small subplot

    for ist in Storms:
        # loop thru each storm
        fc_init = fc_iniT( ist )
        fc_end = t_range_model( ist )
        # loop thru each MP
        for imp in MP:
            iCtg = MP.index(imp) #category
            for ida in DA[1:]:
                for it in range(len(fc_init)):
                    print(ist+': '+ida+'_'+imp+' '+fc_init[it])
                    if not os.path.exists( big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida]+'/wrf_df/'+fc_init[it] ):
                        print('Forecast does not exist!')
                        continue
                    # plot
                    d_fc = read_forecast(ist,big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida]+'/wrf_df/'+fc_init[it]) #,fc_init[it], DF_model_end)
                    plot_one_diff( ist,imp, conv[ist][imp][fc_init[it]], ax[ist][imp], var, d_fc, Color_set['c'+str(iCtg)][it], Line_types[ida], 2, fc_init[it], steps=1 )
                    if ist == 'JOSE' and imp == "THO":
                        plot_one_diff( ist,imp, conv[ist][imp][fc_init[it]], subax, var, d_fc, Color_set['c'+str(iCtg)][it], Line_types[ida], 1, fc_init[it], steps=1 )
                        subax.set_yticks( [-5,-2.5,0] )
                        subax.set_ylim([-5,1])
                        subax.set_yticklabels(['-5.0','-2.5','0.0'],fontsize=6)
                        subax.set_xticklabels([''])

    # Legends
    for ist in Storms:
        fc_init = fc_iniT( ist )
        fc_times = copy.deepcopy( fc_init )

        lines = ax[ist]['WSM6'].get_lines()
        lgd_0 = [it for it in fc_times]
        legend0 = ax[ist]['WSM6'].legend([lines[i] for i in [0,1,2,3]], lgd_0,fontsize='7',loc='lower center',ncol=2)
        legend0.set_alpha( 0.5 )
        # Add the first legend manually to the current Axes
        ax[ist]['WSM6'].add_artist(legend0)

        lines = ax[ist]['THO'].get_lines()
        lgd_0 = [it for it in fc_times]
        legend0 = ax[ist]['THO'].legend([lines[i] for i in [0,1,2,3]], lgd_0,fontsize='7',loc='upper center',ncol=2)
        legend0.set_alpha( 0.5 )

    # Emphasize y =0 line
    for ist in Storms:
        for imp in MP:
            ax[ist][imp].axhline(y=0, color='gray', linestyle='-',linewidth=1.5,alpha=0.5)

    # Create proxy artists for the second legend
    proxy_artist1 = plt.Line2D((0, 1), (0, 0), color=Color_set['c0'][0], linestyle='-')
    proxy_artist2 = plt.Line2D((0, 1), (0, 0), color=Color_set['c0'][0], linestyle=':')
    ax['HARVEY']['WSM6'].legend([proxy_artist1, proxy_artist2], ['IR', 'IR+MW'],fontsize='8',loc='upper center',ncol=3)

    # Set ticks/other attributes
    date_form = mdates.DateFormatter("%m-%d")
    for ist in Storms:
        # loop thru each stormplot_time( Storm )
        plot_st, plot_end = plot_time( ist )
        for imp in MP:
            # x axis
            ax[ist][imp].set_xlim([datetime(int(plot_st[0:4]),int(plot_st[4:6]), int(plot_st[6:8]), int(plot_st[8:10])), datetime(int(plot_end[0:4]),int(plot_end[4:6]), int(plot_end[6:8]), int(plot_end[8:10]))])
            if MP.index( imp ) == 0:
                if ist == 'HARVEY':
                    plot_t = ['201708230000','201708240000','201708250000','201708260000']
                elif ist == 'JOSE':
                    plot_t = ['201709050000','201709060000','201709070000','201709080000','201709090000']
                elif ist == 'IRMA':
                    plot_t = ['201709030000','201709040000','201709050000','201709060000','201709070000']
                elif ist == 'MARIA':
                    plot_t = ['201709160000','201709170000','201709180000','201709190000','201709200000']
                else:
                    pass
                tick_locations = [datetime.strptime(i,"%Y%m%d%H%M") for i in plot_t]
                # Define custom tick labelsi
                tick_labels = [date.strftime('%m-%d') for date in tick_locations]
                # Set x labels
                ax[ist][imp].set_xticks(tick_locations)
                ax[ist][imp].set_xticklabels(tick_labels)
            else:
                # Customize x-axis tick labels precisely
                ax[ist][imp].xaxis.set_major_locator(mdates.DayLocator())
                ax[ist][imp].xaxis.set_major_formatter(date_form)
                ax[ist][imp].tick_params(axis='x',labelsize=8)
            # y axis
            #ax[ist][imp].set_yscale('log')
            ax[ist][imp].set_ylim([-1.75,1.0]) # X*1e6
            if MP.index(imp) != 0:
                 ax[ist][imp].set_yticklabels([])
            ax[ist][imp].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    # Add y label
    fig.text(0.04,0.87,'HARVEY', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.04,0.65,'JOSE', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.04,0.44,'IRMA', fontsize=10, ha='center', va='center' ,rotation='vertical')
    fig.text(0.04,0.21,'MARIA', fontsize=10, ha='center', va='center',rotation='vertical')

    fig.text(0.55,0.05,'Difference of area total grid scale precipitation ($10^6$ mm):\n Experiment - CONV', fontsize=10, ha='center', va='center')

    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_precip_diff.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY','JOSE','IRMA','MARIA']
    DA = ['CONV','IR','IR+MW']
    MP = ['WSM6','THO'] #

    varnames = ['Precip',]

    # if operate over the same number of samples for all forecasts
    sameNum_sample = False
    if sameNum_sample:
        fc_run_hrs = 60

    single_storm = False
    distinct_colors = False
    
    #------------------------------------
    wrf_dir = big_dir

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = {}
        for imp in MP:
            Exper_names[istorm][imp] = {}
            for ida in DA:
                Exper_names[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # Plot
    start_time=time.process_time()

    for var in varnames:

        if single_storm:
            plot_single_storm( Storms[0],var )
        else:
            plot_4by2_diff()    

    end_time = time.process_time()
    print('time needed: ', end_time-start_time, ' seconds')






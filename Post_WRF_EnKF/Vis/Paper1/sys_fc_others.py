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
import sys_EnKF_minSlp as SC #StormCenter
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

def storm_color():
    colorset = {'HARVEY':'#FFA500','IRMA':'#FF13EC','JOSE':'#0D13F6','MARIA':'#097969',}
    return colorset

def alpha_fc_init():
    alphas = np.linspace(0.2,1,4)
    return alphas

def DA_marker():
    marker_type= {'CONV':'P','IR':'o','IR+MW':'s'}
    return marker_type

# ------------------------------------------------------------------------------------------------------
#            Operation: Data Processing
# ------------------------------------------------------------------------------------------------------

# Only read forecasts from one initialization time
# Read field from forecasts listed with wildcard
# Without constraints on names
def read_listed_forecast(storm,exp_dir):

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

# Read forecasts from all initialization times
# Sample number are different across different forecasts 
# I.e., sample numbers are decided by the start of, end of the forecast time, and interval
def read_fc_eachInit( storm,exp_name,exper_inits,DF_model_end ):

    d_fc = {}
    # Loop over initialization time for each experiment
    for init_t in exper_inits:
        fc_dir = wrf_dir+storm+'/'+exp_name+'/wrf_df/'

        # initialize container
        d_fc[init_t] = {var:[] for var in varnames}
        # generate the forecast time series every 3 hours
        fc_times_str = UD.generate_times(init_t,DF_model_end,3)
        run_times = len(fc_times_str)

        # read thru each forecast
        for it in fc_times_str:
            wrfout = ( fc_dir + init_t +
                "/wrfout_d03_" +
                it[:4] + "-" +  # Year
                it[4:6] + "-" +  # Month
                it[6:8] + "_" +  # Day
                it[8:10] + ":" +  # Hour
                it[10:12] + ":00"  # Minute and seconds
            ) 
            # loop thru vars
            for var in varnames:
                if var == 'Precip':
                    ncdir = nc.Dataset( wrfout )
                    ppb = ncdir.variables['RAINNC'][0,:,:]
                    # Make precipitation related to convection stand out
                    #mask = ppb <= 0.5
                    #ppb_masked = ma.masked_array(ppb, mask=mask)
                    #x_pp.append( np.ma.mean(ppb) )

                d_fc[init_t][var].append( (np.sum(ppb)/1e6) )  # /1e6

        # convert to numpy array
        d_fc[init_t][var] = np.array( d_fc[init_t][var] )

    return d_fc


# Exper_vars: storms * microphysics * DA * fc_init 
# Store forecasts from all initialization times
def Read_allFCs():

    # Read vars with lead times of various length
    Exper_vars = {}
    for istorm in Storms:
        Exper_vars[istorm] = {}
        # Time range for model data
        Exper_initTimes = fc_iniT( istorm )
        DF_model_end = t_range_model( istorm )
        for imp in MP:
            Exper_vars[istorm][imp] = {}
            for ida in DA:
                Exper_vars[istorm][imp][ida] = {}
                iExper = Exper_names[istorm][imp][ida]
                if os.path.exists( wrf_dir+'/'+istorm+'/'+iExper ):
                    print('Reading variable value for '+istorm+' :'+iExper)
                    if sameNum_sample:
                        pass
                    else:
                        Exper_vars[istorm][imp][ida] = read_fc_eachInit(istorm,iExper,Exper_initTimes,DF_model_end)
                else:
                    Exper_vars[istorm][imp][ida] = None

    return Exper_vars

# Exper_vars: storms * microphysics * DA * fc_init 
def Time_accu_means():

    Exper_vars = Read_allFCs()

    # find the shortest forecast period
    fc_srt_len = {}
    for ist in Storms:
        fc_inits = fc_iniT( ist )
        fc_lens = [ len( Exper_vars[ist][MP[0]][DA[0]][fc_init][varnames[0]] ) for fc_init in fc_inits]
        fc_srt_len[ist] = np.min( fc_lens )

    # Average over storms
    Means = {}
    for istorm in Storms:
        fc_inits = fc_iniT( istorm )
        Means[istorm] = {}
        for imp in MP:
            Means[istorm][imp] = {}
            for ida in DA:
                Means[istorm][imp][ida] = {}
                for var in varnames:
                    # assemble values of forecasted attributes of len(fc_inits)
                    across_value = np.empty((1,fc_srt_len[istorm]))
                    for fc_init in fc_inits:
                        # only assemble values to the end of the shortest forecast period
                        fc_values = Exper_vars[istorm][imp][ida][fc_init][var][:fc_srt_len[istorm]]
                        accumulated_sums = [sum(fc_values[:i+1]) for i in range(len(fc_values))]
                        across_value = np.concatenate( (across_value,np.array(accumulated_sums).reshape(1,fc_srt_len[istorm])),axis=0)
                across_value = across_value[1:,:] # remove the first empty row
                # average over all forecasts for that forecast kind
                Means[istorm][imp][ida][var] = np.nanmean(across_value,axis=0)

    return fc_srt_len, Means 


def MEAN_overStorms():

    mean_sts = {}
    for imp in MP:
        mean_sts[imp] = {}
        for ida in DA:
            mean_sts[imp][ida] = {}
            for var in varnames:
                # initialize containers
                across_value = np.full( (len(Storms),max(fc_srt_len.values())), np.nan)
                # loop thru storms
                for ist in Storms:
                    across_value[Storms.index(ist),:fc_srt_len[ist]] = Means[ist][imp][ida][var]
                # average
                mean_sts[imp][ida][var] = np.nanmean(across_value,axis=0)

    return mean_sts

# ------------------------------------------------------------------------------------------------------
#            Operation: Plot RAINNC (accumulated grid scale precipitation)
# ------------------------------------------------------------------------------------------------------

# layout
# WSM6: time-accumulated rain, WSM6: time-accumulated rain diff
# THO: time-accumulated rain , THO: time-accumulated rain diff
def plot_2by2_Means():

    # Set up figure
    fig = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=2,top=0.93,right=0.96,hspace=0.03)

    ax = {}
    for imp in MP:
        ax[imp] = {}
        ir = MP.index(imp)
        ist_grids = outer_grids[ir].subgridspec(1, 2, wspace=0.2)
        ax[imp]['orig'] = fig.add_subplot( ist_grids[0] )
        ax[imp]['diff'] = fig.add_subplot( ist_grids[1] )

    # Customization
    colorset = storm_color()
    #marker_type= DA_marker()
    line_width = {'HARVEY':1.5,'IRMA':1.5,'JOSE':1.5,'MARIA':1.5, 'Mean':2}
    Line_types = {'CONV':'-','IR':'--','MW':':'}
    alphas = {'CONV':0.5,'IR':0.7,'MW':0.9}

    # Plot simulations
    for imp in MP:
        for ist in Storms:
            lead_t = list(range(0,fc_srt_len[ist]))
            # Plot original value
            for ida in DA:
                ax[imp]['orig'].plot(lead_t,Means[ist][imp][ida]['Precip'],Line_types[ida],color=colorset[ist],linewidth=line_width[ist],alpha=alphas[ida])
            # plot Exper - CONV
            ax[imp]['diff'].axhline(y=0, color='gray', linestyle='-',linewidth=1.5,alpha=0.5)
            for ida in DA[1:]:
                ax[imp]['diff'].plot(lead_t,Means[ist][imp][ida]['Precip']-Means[ist][imp]['CONV']['Precip'],Line_types[ida],color=colorset[ist],linewidth=2,alpha=alphas[ida])

    # mean over storms
    mean_sts = MEAN_overStorms()
    lead_t = list(range(0,max(fc_srt_len.values())))
    for imp in MP:
        # Plot original value
        for ida in DA:
            ax[imp]['orig'].plot(lead_t,mean_sts[imp][ida]['Precip'],Line_types[ida],color='black',linewidth=2,alpha=alphas[ida])
        # Plot Exper - CONV
        #for ida in DA[1:]:
        #    ax[imp]['diff'].plot(lead_t,mean_sts[imp][ida]['Precip']-mean_sts[imp]['CONV']['Precip'],Line_types[ida],color='black',linewidth=2,alpha=alphas[ida])


    # Manully add legends
    # create proxy artists for the legend with different line widths and colors
    lgd_1 = Storms + ['Mean']
    legend_colors = list(colorset.values())
    legend_colors.append( 'black' )
    list_widths = list(line_width.values())
    proxy_artists = [plt.Line2D([0], [0], color=color, lw=lw) for color,lw in zip( legend_colors,list_widths )]
    #legend0.set_alpha( 0.5 )fig.text(0.5, 0.5, 'Time-accumulated differences of precipitation: Experiment - CONV', ha='center', va='center', rotation='vertical', fontsize=10)
    # Add the first legend manually to the current Axes
    ax['THO']['orig'].legend(proxy_artists,lgd_1,fontsize='8',loc='upper left',ncol=1)

    lines = ax['WSM6']['orig'].get_lines()
    lgd_0 = DA
    legend0 = ax['WSM6']['orig'].legend([lines[i] for i in [-3,-2,-1]], lgd_0,fontsize='8',loc='upper left')
    #legend0.set_alpha( 0.5 )
    # Add the first legend manually to the current Axes
    ax['WSM6']['orig'].add_artist(legend0)

    # Set axis attributes
    for imp in MP:
        # x axis
        ax[imp]['orig'].set_xlim( [-0.1,max(fc_srt_len.values())-1] )
        ax[imp]['diff'].set_xlim( [-0.1,max(fc_srt_len.values())-1] )
        ax[imp]['orig'].set_xticks([0,8,16,24,32] )
        ax[imp]['diff'].set_xticks([0,8,16,24,32] )
        if MP.index(imp) != 0:
            ax[imp]['orig'].set_xticklabels(['D0','D1','D2','D3','D4'])
            ax[imp]['diff'].set_xticklabels(['D0','D1','D2','D3','D4'])
            ax[imp]['orig'].set_xlabel('Forecast Time (days)')
            ax[imp]['diff'].set_xlabel('Forecast Time (days)')
        else:
            ax[imp]['orig'].set_xticklabels(['','','','',''])
            ax[imp]['diff'].set_xticklabels(['','','','',''])
        # y axis
        ax[imp]['orig'].set_ylim( [-0.1,85.1] )
        ax[imp]['diff'].set_ylim( [-30.1,10.1] )
        # y ticks
        diff_yticks = list(range(-30,10+1,5))
        ax[imp]['diff'].set_yticks( diff_yticks )
        ax[imp]['diff'].set_yticklabels( [str(it) for it in diff_yticks] )
        orig_yticks = list(range(0,85+1,10))
        ax[imp]['orig'].set_yticks( orig_yticks)
        ax[imp]['orig'].set_yticklabels( [str(it) for it in orig_yticks] )
        # y label
        #ax[imp]['orig'].set_ylabel(imp,fontsize=10)
        #ax[imp]['mslp'].set_ylabel('MSLP: MAE (hPa)')
        # grid lines 
        ax[imp]['orig'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')
        ax[imp]['diff'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    fig.text(0.03,0.73,'WSM6', fontsize=12, ha='center', va='center',rotation='vertical')
    fig.text(0.03,0.32,'THO', fontsize=12, ha='center', va='center',rotation='vertical')
    fig.text(0.07, 0.5, 'Time-accumulated Precip (TAP; $10^6$ mm)', ha='center', va='center', rotation='vertical', fontsize=10)
    fig.text(0.53, 0.5, 'Differences of TAP: Experiment - CONV ($10^6$ mm)', ha='center', va='center', rotation='vertical', fontsize=10)

    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_time_accumulated_fc_precip_withFCtime.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )




def plot_one( ist, ida, ax, var, state, color, line, line_width, label, steps=3 ):

    # Process data to be plotted
    if (ist == 'HARVEY') and (ida != 'CONV'): # Jerry's run. He saved forecasts every three hours. 
        times = state['time']
        v_to_plot = state[var]
    else:
        times = state['time'][::steps]
        v_to_plot = state[var][::steps]

    dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]

    if time_accumulated:
    # For each index i in the list, it calculates the sum of all elements from the start of the list up to (and including) index i. 
        accumulated_sums = [sum(v_to_plot[:i+1]) for i in range(len(v_to_plot))]
        ax.plot(dates,accumulated_sums,color=color,linestyle=line,linewidth=line_width)
    else:
        ax.plot(dates,v_to_plot,color=color,linestyle=line,linewidth=line_width) # instead of plot_date

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
                d_fc = read_listed_forecast(Storm,big_dir+'/'+Storm+'/'+Exper_names[Storm][imp][ida]+'/wrf_df/'+fc_init[it]) #,fc_init[it], DF_model_end)
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


# rows: storms; columns: DA methods
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
                    d_fc = read_listed_forecast(ist,big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida]+'/wrf_df/'+fc_init[it]) #,fc_init[it], DF_model_end)
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

# rows: storms; columns: microphysics schemes
def plot_4by2(  ):

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
    Line_types = {'CONV':'-','IR':'--','MW':(0, (1, 1))}

    # Plot simulations
    for ist in Storms:
        # loop thru each storm
        fc_init = fc_iniT( ist )
        fc_end = t_range_model( ist ) 
        # loop thru each MP
        for imp in MP:
            iCtg = MP.index(imp) #category
            for ida in DA:
                for it in range(len(fc_init)):
                    print(ist+': '+ida+'_'+imp+' '+fc_init[it])
                    if not os.path.exists( big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida]+'/wrf_df/'+fc_init[it] ):
                        print('Forecast does not exist!')
                        continue
                    # plot
                    d_fc = read_listed_forecast(ist,big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida]+'/wrf_df/'+fc_init[it]) #,fc_init[it], DF_model_end)
                    # Plot forecasted value every 3 hours (Jerry saved his HARVEY forecasts in this way.)
                    plot_one( ist, ida, ax[ist][imp], var, d_fc, Color_set['c'+str(iCtg)][it], Line_types[ida], 1.5, fc_init[it], steps=3 ) 

    # Legends
    for ist in Storms:
        fc_init = fc_iniT( ist )
        fc_times = copy.deepcopy( fc_init )

        lines = ax[ist]['WSM6'].get_lines()
        lgd_0 = [it for it in fc_times]
        legend0 = ax[ist]['WSM6'].legend([lines[i] for i in [0,1,2,3]], lgd_0,fontsize='7',loc='upper left')
        legend0.set_alpha( 0.5 )
        # Add the first legend manually to the current Axes
        ax[ist]['WSM6'].add_artist(legend0)
        
        lines = ax[ist]['THO'].get_lines()
        lgd_0 = [it for it in fc_times]
        legend0 = ax[ist]['THO'].legend([lines[i] for i in [0,1,2,3]], lgd_0,fontsize='7',loc='upper left')
        legend0.set_alpha( 0.5 )

    # Create proxy artists for the second legend
    proxy_artist1 = plt.Line2D((0, 1), (0, 0), color=Color_set['c0'][0], linestyle='-')
    proxy_artist2 = plt.Line2D((0, 1), (0, 0), color=Color_set['c0'][0], linestyle='--')
    proxy_artist3 = plt.Line2D((0, 1), (0, 0), color=Color_set['c0'][0], linestyle=':')
    ax['HARVEY']['WSM6'].legend([proxy_artist1, proxy_artist2, proxy_artist3], ['CONV', 'IR', 'IR+MW'],fontsize='7',loc='lower right')

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
            if time_accumulated:
                ax[ist][imp].set_ylim([0,120]) # X*1e6
            else:
                ax[ist][imp].set_ylim([0,7]) # X*1e6
            if MP.index(imp) != 0:
                 ax[ist][imp].set_yticklabels([])
            ax[ist][imp].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    # Add y label
    fig.text(0.06,0.87,'HARVEY', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.06,0.65,'JOSE', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.06,0.44,'IRMA', fontsize=10, ha='center', va='center' ,rotation='vertical')
    fig.text(0.06,0.21,'MARIA', fontsize=10, ha='center', va='center',rotation='vertical')

    if time_accumulated:
        fig.text(0.55,0.05,'Every 3 hours: time-accumulated area total grid scale precipitation ($10^6$ mm)', fontsize=10, ha='center', va='center')
    else:
        fig.text(0.55,0.05,'Every 3 hours: area total grid scale precipitation ($10^6$ mm)', fontsize=10, ha='center', va='center')

    # Save figure
    if time_accumulated:
        des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/Time_accu_sys_fc_precip.png'
    else:
        des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_precip.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

# layout: 
# top row: CONV (time-accumulated rain)
# bottom row: Exp - CONV (time-accumulated rain diff)
def plot_conv_diff():

    # Set up figure
    fig = plt.figure( figsize=(6.5,7),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=2,top=0.93,right=0.96,hspace=0.25)
    # set up axis
    ax = {}
    ax['orig'] = fig.add_subplot( outer_grids[0] )
    ax['diff'] = {}
    ist_grids = outer_grids[1].subgridspec(1, 2, wspace=0.2)
    for imp in MP:
        ir = MP.index(imp)
        ax['diff'][imp] = fig.add_subplot( ist_grids[ir] )

    # Customization
    colorset = storm_color()
    #marker_type= DA_marker()
    line_width = {'HARVEY':1.5,'IRMA':1.5,'JOSE':1.5,'MARIA':1.5, 'Mean':2}
    Line_types_DA = {'IR':'-','MW':'--'}
    Line_types_MP = {'WSM6':'-','THO':'--',}
    alphas = {'CONV':0.5,'IR':0.7,'MW':0.9}

    # Plot CONV simulations
    # mean over storms
    mean_sts = MEAN_overStorms()
    lead_t = list(range(0,max(fc_srt_len.values())))

    for imp in MP:
        # Plot means for each storm
        for ist in Storms:
            lead_t = list(range(0,fc_srt_len[ist]))
            ax['orig'].plot(lead_t,Means[ist][imp]['CONV']['Precip'],Line_types_MP[imp],color=colorset[ist],linewidth=line_width[ist],)
        # Plot means for either WSM6 or THO experiments
        ax['orig'].plot(lead_t,mean_sts[imp]['CONV']['Precip'],Line_types_MP[imp],color='black',linewidth=2)

    # Plot Exper - CONV
    for imp in MP:
        for ist in Storms:
            lead_t = list(range(0,fc_srt_len[ist]))
            ax['diff'][imp].axhline(y=0, color='gray', linestyle='-',linewidth=1.5,alpha=0.5)
            for ida in DA[1:]:
                ax['diff'][imp].plot(lead_t,Means[ist][imp][ida]['Precip']-Means[ist][imp]['CONV']['Precip'],Line_types_DA[ida],color=colorset[ist],linewidth=2,alpha=alphas[ida])

    # Manully add legends
    # First row
    lgd_00 = Storms + ['Mean']
    legend_colors = list(colorset.values())
    legend_colors.append( 'black' )
    list_widths = list(line_width.values())
    proxy_artists00 = [plt.Line2D([0], [0], color=color, lw=lw) for color,lw in zip( legend_colors,list_widths )]
    # Add the first legend manually to the current Axes
    legend00 = ax['orig'].legend(proxy_artists00,lgd_00,fontsize='10',loc='upper left')
    ax['orig'].add_artist(legend00)

    lgd_01 = MP
    list_types_MP = list(Line_types_MP.values())
    proxy_artists01 = [plt.Line2D([0], [0], color='black', linestyle=lt) for lt in list_types_MP]
    legend01 = ax['orig'].legend(proxy_artists01,lgd_01,fontsize='10',loc='lower right')
    ax['orig'].add_artist(legend01)

    # Second row
    lgd_10 = Storms
    legend_colors = list(colorset.values())
    proxy_artists10 = [plt.Line2D([0], [0], color=color ) for color in legend_colors ]
    legend10 = ax['diff'][MP[0]].legend(proxy_artists10,lgd_10,fontsize='8',loc='lower left',ncol=1)
    ax['diff'][MP[0]].add_artist(legend10)
    legend10 = ax['diff'][MP[1]].legend(proxy_artists10,lgd_10,fontsize='8',loc='lower left',ncol=1)
    ax['diff'][MP[1]].add_artist(legend10)

    lgd_11 = DA[1:]
    list_types = list(Line_types_DA.values())
    proxy_artists11 = [plt.Line2D([0], [0], color='black', linestyle=lt) for lt in list_types ]
    legend11 = ax['diff'][MP[0]].legend(proxy_artists11,lgd_11,fontsize='8',loc='upper left')
    ax['diff'][MP[0]].add_artist(legend11)
    legend11 = ax['diff'][MP[1]].legend(proxy_artists11,lgd_11,fontsize='8',loc='upper left')
    ax['diff'][MP[1]].add_artist(legend11)

    # Set axis attributes
    ax['orig'].set_title(' (a) Time-accumulated Precip (TAP) for CONV Exps',fontsize='12')
    ax['orig'].set_xlim( [-0.1,max(fc_srt_len.values())-1] )
    ax['orig'].set_xticks([0,8,16,24,32] )
    ax['orig'].set_xticklabels(['D0','D1','D2','D3','D4'])
    ax['orig'].set_xlabel('Forecast Time (days)',fontsize=9)
    ax['orig'].set_ylabel('CONV: Precip ($10^6$ mm)',fontsize=10)
    # y axis
    ax['orig'].set_ylim( [-0.1,85.1] )
    # y ticks
    orig_yticks = list(range(0,85+1,10))
    ax['orig'].set_yticks( orig_yticks)
    ax['orig'].set_yticklabels( [str(it) for it in orig_yticks],fontsize=9 )
    # grid lines
    ax['orig'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')
    for imp in MP:
        # x axis
        ax['diff'][imp].set_xlim( [-0.1,max(fc_srt_len.values())-1] )
        ax['diff'][imp].set_xticks([0,8,16,24,32] )
        ax['diff'][imp].set_xticklabels(['D0','D1','D2','D3','D4'])
        ax['diff'][imp].set_xlabel('Forecast Time (days)',fontsize=9)
        # y axis
        ax['diff'][imp].set_ylim( [-30.1,10.1] )
        # y ticks
        diff_yticks = list(range(-30,10+1,5))
        ax['diff'][imp].set_yticks( diff_yticks )
        ax['diff'][imp].set_yticklabels( [str(it) for it in diff_yticks],fontsize=9 )
        # y label
        ax['diff'][imp].set_ylabel(imp+' - CONV: Precip ($10^6$ mm)',fontsize=10)
        # grid lines
        ax['diff'][imp].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    fig.text(0.53, 0.5, '(b) Differences of TAP: DA Exps - CONV Exps', ha='center', va='center', fontsize=12)

    #fig.text(0.03,0.73,'WSM6', fontsize=12, ha='center', va='center',rotation='vertical')
    #fig.text(0.03,0.32,'THO', fontsize=12, ha='center', va='center',rotation='vertical')
    #fig.text(0.07, 0.5, 'Time-accumulated Precip (TAP; $10^6$ mm)', ha='center', va='center', rotation='vertical', fontsize=10)
    #fig.text(0.53, 0.5, 'Differences of TAP: Experiment - CONV ($10^6$ mm)', ha='center', va='center', rotation='vertical', fontsize=10)

    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_time_accumulated_fc_precip_withFCtime.png'
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
    height = rbox.height
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


def plot_one_diff( ist, ida, conv_each, ax, var, state, color, line, line_width, label, steps=3 ):

    # Process data to be plotted
    if ist == 'HARVEY': # Jerry's run. He saved forecasts every three hours. 
        conv_plot = conv_each[var][::steps] # Harvey conv_each is stored hourly
        
        times = state['time'] # Harvey other is stored every 3 hours
        state_plot = state[var] # Harvey other is stored every 3 hours
    else:
        conv_plot = conv_each[var][::steps] 
        times = state['time'][::steps]
        state_plot = state[var][::steps]

    dates = [datetime.strptime(i,"%Y%m%d%H%M") for i in times]

    if time_accumulated:
        v_to_plot = state_plot-conv_plot
        accumulated_sums = [sum(v_to_plot[:i+1]) for i in range(len(v_to_plot))]
        ax.plot(dates,accumulated_sums,color=color,linestyle=line,linewidth=line_width)
    else:
        v_to_plot = state_plot-conv_plot
        ax.plot(dates,v_to_plot,color=color,linestyle=line,linewidth=line_width) # instead of plot_date

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
                conv[ist][imp][it] = read_listed_forecast(ist,big_dir+'/'+ist+'/'+Exper_names[ist][imp]['CONV']+'/wrf_df/'+it) 

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
                    d_fc = read_listed_forecast(ist,big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida]+'/wrf_df/'+fc_init[it]) #,fc_init[it], DF_model_end)
                    plot_one_diff( ist,ida, conv[ist][imp][fc_init[it]], ax[ist][imp], var, d_fc, Color_set['c'+str(iCtg)][it], Line_types[ida], 2, fc_init[it], steps=3 )
                    if ist == 'JOSE' and imp == "THO":
                        plot_one_diff( ist,ida, conv[ist][imp][fc_init[it]], subax, var, d_fc, Color_set['c'+str(iCtg)][it], Line_types[ida], 1, fc_init[it], steps=3 )
                        if time_accumulated:
                            subax.set_yticks( [-90,-45,0] )
                            subax.set_ylim([-90,1])
                            subax.set_yticklabels(['-90','-45','0'],fontsize=6)
                            subax.set_xticklabels([''])
                        else:
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
            if time_accumulated:
                ax[ist][imp].set_ylim([-30,15]) # X*1e6
            else:
                ax[ist][imp].set_ylim([-1.75,1.0]) # X*1e6
            if MP.index(imp) != 0:
                 ax[ist][imp].set_yticklabels([])
            ax[ist][imp].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    # Add y label
    fig.text(0.04,0.87,'HARVEY', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.04,0.65,'JOSE', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.04,0.44,'IRMA', fontsize=10, ha='center', va='center' ,rotation='vertical')
    fig.text(0.04,0.21,'MARIA', fontsize=10, ha='center', va='center',rotation='vertical')

    if time_accumulated:
        fig.text(0.55,0.05,'Every 3 hours: time-accumulated difference \n of area total grid scale precipitation ($10^6$ mm):\n Experiment - CONV', fontsize=10, ha='center', va='center')
    else:
        fig.text(0.55,0.05,'Every 3 hours: Difference of area total grid scale precipitation ($10^6$ mm):\n Experiment - CONV', fontsize=10, ha='center', va='center')

    # Save figure
    if time_accumulated:
        des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/Time_accu_sys_fc_precip_diff.png'
    else:
        des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_precip_diff.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/Clean_results/'

    #--------Configuration------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    DA = ['CONV','IR','MW']
    MP = ['WSM6','THO'] #

    varnames = ['Precip',]

    # if operate over the same number of samples for all forecasts
    sameNum_sample = False
    if sameNum_sample:
        fc_run_hrs = 60

    single_storm = False
    multi_storms = False
    distinct_colors = False
    difference = False
    time_accumulated = False
    Mean_timeAccu_wrt_time = True

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

    # Calculate the time-accumulated values and average over storms
    # Fig should look like MAE_wrt_time
    if Mean_timeAccu_wrt_time:
        # calculate with same samples
        fc_srt_len, Means = Time_accu_means()
        #plot_2by2_Means()
        plot_conv_diff()

    # Plot
    start_time=time.process_time()

    for var in varnames:

        if single_storm:
            plot_single_storm( Storms[0],var )
        
        if multi_storms:
            #plot_4by3(  )
            plot_4by2(  )

        if difference:
            plot_4by2_diff()    

    end_time = time.process_time()
    print('time needed: ', end_time-start_time, ' seconds')






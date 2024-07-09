import os,fnmatch # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import pickle
import netCDF4 as nc
from wrf import getvar
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
import copy
import matplotlib.dates as mdates
from matplotlib import pyplot
import numpy.ma as ma
from scipy import interpolate

import EnKF_Vmax as Vx #Vmax
import EnKF_minSLP_track as SC #StormCenter
import Util_data as UD

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['lines.linewidth'] = 2.5#1.5
matplotlib.rcParams['lines.markersize'] = 2.5
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['font.size'] = 8


# Calculate mean absolute error
def mean_absolute_error(obs,model):
    return np.mean(np.abs(np.subtract(obs, model)))

# Generate time series
def generate_times( Storms, start_time_str, end_time_str, interval ):

    dict_times = {}
    for istorm in Storms:
        time_diff = datetime.strptime(end_time_str[istorm],"%Y%m%d%H%M") - datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+interval, interval))]
        dict_times[istorm] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times

# --------------------------------------------------------
#    HPI
# --------------------------------------------------------

# stack obs over storms
def stack_obs( d_obs,iv ):

    if iv == 'track':
        stacked = [[],[]] 
        for ist in Storms:
            stacked[0].extend( d_obs[ist][iv][0] ) # lon
            stacked[1].extend( d_obs[ist][iv][1] ) # lat
    else:
        stacked = []
        for ist in Storms:
            stacked.extend( d_obs[ist][iv] )

    return stacked

# stack simulations over storms
def stack_model( d_model ):

    stacked_model = {}
    for imp in MP:
        stacked_model[imp] = {}
        for ida in DA:
            stacked_model[imp][ida] = {}
            for iv in var:
                if iv == 'track':
                    stacked = [[],[]]
                    for ist in Storms:
                        stacked[0].extend( d_model[ist][imp][ida][iv][0] ) # lon
                        stacked[1].extend( d_model[ist][imp][ida][iv][1] ) # lat 
                else:
                    stacked = []
                    for ist in Storms:
                        stacked.extend( d_model[ist][imp][ida][iv] )
                # fill in value
                stacked_model[imp][ida][iv] = stacked
    return stacked_model


# bin the bias
def bin_bias( d_all_obs,d_all_model ):

    # calculate the bias
    bias = {}
    for imp in MP:
        bias[imp] = {}
        for ida in DA:
            bias[imp][ida] = {}
            for iv in var:
                if iv == 'track':
                    distance = []
                    # loop thru each DA time
                    for idx in range( len(d_all_obs[iv][0]) ):
                        lon_x = d_all_model[imp][ida][iv][0][idx]
                        lat_x = d_all_model[imp][ida][iv][1][idx]
                        lon_obs = d_all_obs[iv][0][idx]
                        lat_obs = d_all_obs[iv][1][idx]
                        dist, azi_btk = UD.twoP_inverse( lat_obs,lon_obs,lat_x,lon_x) #km
                        distance.append( dist )
                    bias[imp][ida][iv] = distance
                else:
                    bias[imp][ida][iv] = np.array(d_all_model[imp][ida][iv]) - np.array(d_all_obs[iv]) 

    # bin the bias
    num_bins = 20
    minDist = 0
    maxDist = 200 #km
    minMSLP = -30
    maxMSLP = 20
    minVmax = -25
    maxVmax = 25

    # Marginally-bin data
    hist_var = {}
    for imp in MP:
        hist_var[imp] = {}
        for ida in DA:
            hist_var[imp][ida] = {}
            for iv in var:
                hist_var[imp][ida][iv] = {}
                if iv == 'track':
                    hist, bin_edges = np.histogram(bias[imp][ida][iv],range=[minDist,maxDist],bins=num_bins,density=True)
                elif iv == 'minSLP':
                    hist, bin_edges = np.histogram(bias[imp][ida][iv],range=[minMSLP,maxMSLP],bins=num_bins,density=True)
                elif iv == 'Vmax':
                    hist, bin_edges = np.histogram(bias[imp][ida][iv],range=[minVmax,maxVmax],bins=num_bins,density=True)

                hist_var[imp][ida][iv]['hist'] = hist
                hist_var[imp][ida][iv]['bin_edges'] = bin_edges

    return hist_var

# --------------------------------------------------------
#    Plot HPI
# --------------------------------------------------------

def plot_hist_1by2():

    var_its = ['track','minSLP',]#'Vmax']

    # Set up figure
    fig = plt.figure( figsize=(6.5,4.25),dpi=200) # standard: 6.5,8.5
    grids = fig.add_gridspec(ncols=1,nrows=2,hspace=0.12)
    ax = {}
    for iv in var_its:
        ax[iv] = fig.add_subplot( grids[ var_its.index(iv) ] )

    # customization
    colors = {'WSM6': '#FF3333','THO':'#3333FF'}
    lines = {'CONV':'-','IR':'--','IR+MW':(0, (1, 1))}
    alphas = {'WSM6':1,'THO':0.8}

    # Plot
    for iv in var_its:
        for imp in MP:
            for ida in DA:
                print(iv)
                bin_edges = hist_var[imp][ida][iv]['bin_edges']
                x_axis = (bin_edges[:-1]+bin_edges[1:])/2
                ax[iv].plot(x_axis,hist_var[imp][ida][iv]['hist'],color=colors[imp],linestyle=lines[ida],linewidth='2',alpha=alphas[imp])
                ax[iv].grid(True,linewidth=1, color='gray', alpha=0.5, linestyle='-')

    # Legend
    lgd_1 = [MP[0]+':'+ida for ida in DA]
    lgd_2 = [MP[1]+':'+ida for ida in DA]
    lines_DA = ax[var_its[0]].get_lines()
    legend1 = ax[var_its[0]].legend([lines_DA[i] for i in [0,1,2,]], lgd_1,fontsize='8',loc='upper left')
    legend2 = ax[var_its[1]].legend([lines_DA[i] for i in [3,4,5,]], lgd_2,fontsize='8',loc='upper right') 


    # axes attributes
    for iv in var_its:
        ax[iv].axvline(x=0, color='gray',linewidth=1.5, alpha=0.5)
        ax[iv].set_xlim([np.amin(hist_var[MP[0]][DA[0]][iv]['bin_edges']),np.amax(hist_var[MP[0]][DA[0]][iv]['bin_edges'])])
        ax[iv].set_ylim([-0.001,0.25])
        if iv == 'Vmax':
            x_ticks = np.linspace(np.amin(hist_var[MP[0]][DA[0]][iv]['bin_edges']),np.amax(hist_var[MP[0]][DA[0]][iv]['bin_edges']),5)
            ax[iv].set_xticks( x_ticks )


    # Set titles
    for iv in var_its:
        if iv == 'track':
            ax[iv].set_ylabel( 'Track Bias (km)',fontsize = 12 )
        elif iv == 'minSLP':
            ax[iv].set_ylabel( 'MSLP Bias (hPa)',fontsize = 12 )
        elif iv == 'Vmax':
            ax[iv].set_ylabel( 'Vmax Bias (m $\mathregular{s^{-1}}$)',fontsize = 12 )

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/EnKF_HPI_hist_track.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )



# plot one 
def plot_one( ax0,ax1,state,color,line,line_width,label):
    
    # Read data
    x_min_slp = state['minSLP']
    x_max_ws = state['Vmax']
    # intensity
    ax0.plot(lead_t, x_min_slp, color, linestyle=line, label=label, linewidth=line_width)
    ax1.plot(lead_t, x_max_ws, color, linestyle=line, label=label, linewidth=line_width) 
        
def plot_EnKF_HPI():

    # Customize color
    colors = {'WSM6': '#FF3333','THO':'#3333FF'}
    
    # Customize linestyle
    lines = {'CONV':'-','IR':'--','IR+MW':(0, (1, 1))}  

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200)
    widths = [3.2,3.2]
    grids = fig.add_gridspec(ncols=2,nrows=4,bottom=0.05,top=0.96,left=0.10,right=0.98,wspace=0.08,hspace=0.06)
    #heights = [2.8,2.8,2.8]
    #outer_grid = fig.add_gridspec(ncols=1,nrows=4,bottom=0.05,top=0.96,left=0.10,right=0.98,wspace=0,hspace=0.05)

    ax = {}
    for ist in Storms:
        ax[ist] = {}
        ir = Storms.index( ist )
        # gridspec inside gridspec
        #inner_grid = outer_grid[ir].subgridspec(1, 2, )#width_ratios=widths)
        ax[ist]['ax0'] = fig.add_subplot( grids[ir,0] )
        ax[ist]['ax1'] = fig.add_subplot( grids[ir,1] )
        # plot TC vitals 
        if if_tcvital:
            plot_one( ax[ist]['ax0'], ax[ist]['ax1'],  d_tcvHr[ist],  'gray', '-', 3.5, 'TCvital')
        # plot btk
        if if_btk:
            plot_one( ax[ist]['ax0'], ax[ist]['ax1'],  d_btkHr[ist],  'black', '-', 3.5, 'TCvital')


    # Plot analyses
    for ist in Storms:
        for imp in MP:
            iCtg = MP.index(imp) #category
            for ida in DA:
                if not os.path.exists( big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida] ):
                    print('Experiment does not exist!')
                    continue
                # plot
                plot_one( ax[ist]['ax0'],ax[ist]['ax1'],d_model[ist][imp][ida],colors[imp],lines[ida],2,imp+':'+ida)

    # Manullay control Legend
    lgd = []
    for imp in MP:
        for ida in DA:
            lgd.append(imp+':'+ida)

    if if_btk:
        lgd = ['Best Track'] + lgd
    if if_tcvital:
        lgd = ['TCvitals'] + lgd
    lines_DA0 = ax[Storms[-1]]['ax0'].get_lines()
    legend = ax[Storms[0]]['ax0'].legend(lines_DA0,lgd,ncol=2,fontsize='7.5',loc='lower center')

    #lgd_1 = [MP[0]+':'+ida for ida in DA]
    #lgd_2 = [MP[1]+':'+ida for ida in DA]
    #lines_DA1 = ax[Storms[-1]]['ax1'].get_lines()
    #if if_btk and if_tcvital:
        #legend1 = ax[Storms[0]]['ax0'].legend([lines_DA0[i] for i in [0,1,2,3,4]], lgd_1,fontsize='7',loc='lower left')
        #legend2 = ax[Storms[0]]['ax1'].legend([lines_DA1[i] for i in [5,6,7,]], lgd_2,fontsize='7',loc='lower right')
    #else:
        #legend1 = ax[Storms[0]]['ax0'].legend([lines_DA0[i] for i in [0,1,2,3,]], lgd_1,fontsize='7',loc='lower left')
        #legend2 = ax[Storms[0]]['ax1'].legend([lines_DA1[i] for i in [4,5,6,]], lgd_2,fontsize='7',loc='lower right') 

    # Set ticks/other attributes
    for ist in Storms:
        ax[ist]['ax0'].set_xlim([0,cycles-1])
        ax[ist]['ax1'].set_xlim([0,cycles-1])
        ax[ist]['ax0'].set_xticks( lead_t[::4] )
        ax[ist]['ax1'].set_xticks( lead_t[::4] )
        if Storms.index(ist) < len(Storms)-1:
            ax[ist]['ax0'].set_xticklabels([]) 
            ax[ist]['ax1'].set_xticklabels([])
        # y ticks
        if ist == 'IRMA':
            ax[ist]['ax0'].set_ylim([920,980])
            ax[ist]['ax1'].set_ylim([0,70])
        else:
            ax[ist]['ax0'].set_ylim([990,1020])
            ax[ist]['ax1'].set_ylim([0,40])
        # grids
        ax[ist]['ax0'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')
        ax[ist]['ax1'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    # Set y label
    for ist in Storms:
        if Storms.index(ist) == 0:
            fig.text(0.02,0.85,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 1:
            fig.text(0.02,0.62,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 2:
            fig.text(0.02,0.39,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 3:
            fig.text(0.02,0.16,ist, fontsize=12, ha='center', va='center',rotation='vertical')
    
    # Set x label       
    ax[Storms[-1]]['ax0'].set_xlabel('EnKF Cycle',fontsize=12)      
    ax[Storms[-1]]['ax1'].set_xlabel('EnKF Cycle',fontsize=12)

    # Set panel labels
    for ist in Storms:
        if Storms.index(ist) == 0:
            fig.text(0.13,0.94,'(a1)', fontsize=12, ha='center', va='center')
            fig.text(0.59,0.94,'(a2)', fontsize=12, ha='center', va='center')
        elif Storms.index(ist) == 1:
            fig.text(0.13,0.71,'(b1)', fontsize=12, ha='center', va='center')
            fig.text(0.59,0.71,'(b2)', fontsize=12, ha='center', va='center')
        elif Storms.index(ist) == 2:
            fig.text(0.13,0.48,'(c1)', fontsize=12, ha='center', va='center')
            fig.text(0.59,0.48,'(c2)', fontsize=12, ha='center', va='center')
        elif Storms.index(ist) == 3:
            fig.text(0.13,0.25,'(d1)', fontsize=12, ha='center', va='center')
            fig.text(0.59,0.25,'(d2)', fontsize=12, ha='center', va='center')


    # Set titles
    ax[Storms[0]]['ax0'].set_title( 'MSLP (hPa)',fontsize = 12 )
    ax[Storms[0]]['ax1'].set_title( 'Vmax (m $\mathregular{s^{-1}}$)',fontsize = 12 )

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/EnKF_HPI.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

# --------------------------------------------------------
#   Precipitation 
# --------------------------------------------------------

# Obtain accumulated grid scale precipitation from WRF output
def model_ave_Precip( big_dir,istorm,iExper,DAtimes):

    dict_pp = {}
    x_pp = []
    for DAtime in DAtimes:
        xb_dir = big_dir+istorm+'/'+iExper+'/fc/'+DAtime+'/wrf_enkf_input_d03_mean'
        print('Reading the 1-h ensemble mean forecast from ', xb_dir)
        ncdir = nc.Dataset( xb_dir )
        ppb = ncdir.variables['RAINNC'][0,:,:] 
        # Make precipitation related to convection stand out
        mask = ppb <= 0.5
        ppb_masked = ma.masked_array(ppb, mask=mask)
        x_pp.append( np.ma.mean(ppb) )

    dict_pp['x_precip'] = x_pp

    return dict_pp


def plot_Precip():

    # Customize color
    colors = {'WSM6': '#FF3333','THO':'#3333FF'}

    # Customize linestyle
    lines = {'CONV':'-','IR':'--','IR+MW':(0, (1, 1))}

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200)
    grids = fig.add_gridspec(ncols=1,nrows=4,bottom=0.05,top=0.96,left=0.05,right=0.92,wspace=0.08,hspace=0.08)
    ax = {}
    for ist in Storms:
        ax[ist] = {}
        ir = Storms.index( ist )
        # gridspec inside gridspec
        #inner_grid = outer_grid[ir].subgridspec(1, 2, )#width_ratios=widths)
        ax[ist] = fig.add_subplot( grids[ir] )

    # Plot analyses
    for ist in Storms:
        for imp in MP:
            iCtg = MP.index(imp) #category
            for ida in DA:
                if not os.path.exists( big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida] ):
                    print('Experiment does not exist!')
                    continue
                # plot
                ax[ist].plot(lead_t, np.log10(d_model[ist][imp][ida]['Precip']),colors[imp],linestyle=lines[ida],linewidth=2,label=imp+':'+ida)

    # Legend
    lgd = []
    for imp in MP:
        for ida in DA:
            lgd.append( imp+':'+ida ) 
    lines_DA = ax[Storms[1]].get_lines()
    legend1 = ax[Storms[0]].legend([lines_DA[i] for i in range(len(MP)*len(DA))],lgd,ncol=2,fontsize='8',loc='upper center')

    # tick labels
    ylabel_like = [-0.75,-0.50,-0.25,0.,0.25,0.50,0.75,1.,1.25]
    yticks_loc = ylabel_like
    # Set ticks/other attributes
    for ist in Storms:
        ax[ist].set_xlim([0,cycles-1])
        ax[ist].set_xticks( lead_t[::4] )
        ax[ist].set_ylim([-0.5,1.25])
        if Storms.index(ist) < len(Storms)-1:
            ax[ist].set_xticklabels([])
        # grids
        ax[ist].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')
        ax[ist].yaxis.set_ticks_position('right')
        ax[ist].set_yticks( yticks_loc )
        ytick_label = [fr'$10^{{{i:.2f}}}$' for i in ylabel_like]
        ax[ist].set_yticklabels( ytick_label,fontsize=10 )


    # Set y label
    for ist in Storms:
        if Storms.index(ist) == 0:
            fig.text(0.02,0.85,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 1:
            fig.text(0.02,0.62,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 2:
            fig.text(0.02,0.39,ist, fontsize=12, ha='center', va='center',rotation='vertical')
        elif Storms.index(ist) == 3:
            fig.text(0.02,0.16,ist, fontsize=12, ha='center', va='center',rotation='vertical')

    # Set x label       
    ax[Storms[-1]].set_xlabel('WRF-EnKF Cycle',fontsize=12)

    # Set titles
    ax[Storms[0]].set_title( 'Averaged accumulated grid scale precipitation (mm)',fontsize = 12 )

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/Cycle_ave_Precip.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # Configuration
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    MP = ['WSM6','THO']
    DA = ['CONV','IR','IR+MW']

    var = ['track','minSLP','Vmax']

    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708231200','IRMA':'201709040000','JOSE':'201709060000','MARIA':'201709170000'}
    cycles = 25
    lead_t = list(range(0, cycles, 1))

    # observation
    if_tcvital = False
    if_btk = True

    # EnKF HPI
    if_plot_HPI = True
    plot_histogram = True
    plot_evo_wrt_time = False
    slp_xa = True
    Vmax_xa = True
    slp_xb = False
    Vmax_xb = False

    # EnKF precip
    if_plot_precip = False

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = {}
        for imp in MP:
            Exper_names[istorm][imp] = {}
            for ida in DA:
                Exper_names[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # Identify DA times in the period of interest
    d_hrs = generate_times( Storms, start_time_str, end_time_str, 1 )
    d_6hrs = generate_times( Storms, start_time_str, end_time_str, 6 )

    # Read simulations
    d_model = {}
    for ist in Storms:
        d_model[ist] = {}
        for imp in MP:
            d_model[ist][imp] = {}
            for ida in DA:
                d_model[ist][imp][ida] = {}
                if not os.path.exists( big_dir+'/'+ist+'/'+Exper_names[ist][imp][ida] ):
                    for iv in var:
                        d_model[ist][imp][ida][iv] = None
                else:
                    for iv in var:
                        if iv == 'Vmax':
                            d_vmax = Vx.model_Vmax(big_dir,ist,Exper_names[ist][imp][ida],d_hrs[ist],Vmax_xa,Vmax_xb )
                            d_model[ist][imp][ida][iv] = d_vmax['xa_vmax']
                        elif iv == 'Precip':
                            d_precip = model_ave_Precip( big_dir,ist,Exper_names[ist][imp][ida],d_hrs[ist])
                            d_model[ist][imp][ida][iv] = d_precip['x_precip']
                        else:
                            d_minslp = SC.model_minSLP(big_dir,ist,Exper_names[ist][imp][ida],d_hrs[ist],slp_xa,slp_xb )
                            if iv == 'minSLP':
                                d_model[ist][imp][ida][iv] = d_minslp['xa_slp']
                            elif iv == 'track':
                                d_model[ist][imp][ida][iv] = [d_minslp['xa_lon'],d_minslp['xa_lat']] # lon,lat
    
    # plot HPI
    if if_plot_HPI:
        # Read TCvital data
        d_tcvHr = {}
        # Read best-track data
        d_btkHr = {}
        for istorm in Storms:
            d_tcvHr[istorm] = {}
            d_btkHr[istorm] = {}
            for iv in var:
                if iv == 'Vmax':
                    if if_tcvital:
                        d_tcv6Hrs = Vx.Vmax_TCvitals( small_dir,istorm,d_6hrs[istorm] )
                        d_tcvHr[istorm][iv] = Vx.Interpolate_hourly( istorm,d_6hrs,d_tcv6Hrs,d_hrs)
                    if if_btk:
                        d_btk6Hrs = Vx.Vmax_btk( small_dir,istorm,d_6hrs[istorm] )
                        d_btkHr[istorm][iv] = Vx.Interpolate_hourly( istorm,d_6hrs,d_btk6Hrs,d_hrs)
                else:
                    if if_tcvital:
                        d_tcvital= SC.assimilated_obs( big_dir,istorm,Exper_names[istorm]['WSM6']['CONV'], d_hrs[istorm],slp_xa,slp_xb=False )
                        if iv == 'minSLP':
                            d_tcvHr[istorm][iv] = d_tcvital['minslp_obs']
                        elif iv == 'track':
                            d_tcvHr[istorm][iv] = [d_tcvital['lon_obs'],d_tcvital['lat_obs']] # lon,lat
                    if if_btk:
                        d_btk6Hrs = SC.MSLP_btk( small_dir,istorm,d_6hrs[istorm] )
                        if iv == 'minSLP':
                            d_btkHr[istorm][iv] = Vx.Interpolate_hourly( istorm,d_6hrs,d_btk6Hrs['mslp'],d_hrs)
                        elif iv == 'track':
                            interp_lat = Vx.Interpolate_hourly( istorm,d_6hrs,d_btk6Hrs['lat'],d_hrs)
                            interp_lon = Vx.Interpolate_hourly( istorm,d_6hrs,d_btk6Hrs['lon'],d_hrs)
                            d_btkHr[istorm][iv] = [interp_lon,interp_lat]   # lon, lat

        # plot evolution with respect to time
        if plot_evo_wrt_time:
            plot_EnKF_HPI()

        # plot histogram over all storms 
        if plot_histogram:
            
            # stack obs over storms
            if if_tcvital:
                d_tcvHr_stack = {}
                for iv in var:
                    d_tcvHr_stack[iv] = stack_obs( d_tcvHr,iv )
            if if_btk:
                d_btkHr_stack = {}
                for iv in var:
                    d_btkHr_stack[iv] = stack_obs( d_btkHr,iv )

            # stack simulations over storms
            d_model_stack = stack_model( d_model )

            # calculate the bias: exp - obs
            if if_tcvital:
                hist_var = bin_bias( d_tcvHr_stack,d_model_stack )
            if if_btk:
                hist_var = bin_bias( d_btkHr_stack,d_model_stack )

            # plot
            plot_hist_1by2()

    # Note: EnKF does not update precip!!! 
    if if_plot_precip:
        plot_Precip() 

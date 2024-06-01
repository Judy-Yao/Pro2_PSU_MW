
import os,fnmatch # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
from netCDF4 import Dataset
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
from fast_histogram import histogram2d as hist2d

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

# Best track: start and end of the period
def t_range_btk( Storm ):
    if Storm == 'HARVEY':
        Btk_start = '201708221800' # '201709161800' #'201709030600'
        Btk_end = '201708270000' # '201709210000' #'201709090000'
    elif Storm == 'IRMA':
        Btk_start = '201709030600'
        Btk_end = '201709080000' #201709080000
    elif Storm == 'MARIA':
        Btk_start = '201709160600'#'201709160000'
        Btk_end = '201709210000'
    elif Storm == 'JOSE':
        Btk_start = '201709050600'
        Btk_end = '201709100000'
    else:
        raise ValueError('Storm does not exist!')
    return Btk_start, Btk_end

def calculate_signed_error():

    # Calculate the error with lead times of various length
    Exper_error = {}
    for istorm in Storms:
        Exper_error[istorm] = {}
        # Read best-track data
        Btk_start, Btk_end = t_range_btk( istorm )
        best_track = UD.btk_in_duration(istorm, Btk_start, Btk_end, hour_step=6)
        # Time range for model data
        Exper_initTimes = fc_iniT( istorm )
        DF_model_end = t_range_model( istorm )
        for imp in MP:
            Exper_error[istorm][imp] = {}
            for ida in DA:
                iExper = Exper_names[istorm][imp][ida]
                if os.path.exists( wrf_dir+'/'+istorm+'/'+iExper ):
                    print('Calculating HPI error for '+istorm+' :'+Exper_names[istorm][imp][ida])
                    if sameNum_sample:
                        Exper_error[istorm][imp][ida] = Track.SameL_error_eachInit(istorm,wrf_dir,iExper,best_track,Exper_initTimes,DF_model_end,fc_run_hrs)
                    else:
                        Exper_error[istorm][imp][ida] = Track.DiffL_error_eachInit(istorm,wrf_dir,iExper,best_track,Exper_initTimes,DF_model_end)
                else:
                    Exper_error[istorm][imp][ida] = None
    return Exper_error

# Stack arrays in nested dictionary
def stack_error( Exper_error ):

    allFC_error = {}
    for istorm in Storms:
        allFC_error[istorm] = {}
        for imp in MP:
            allFC_error[istorm][imp] = {}
            for ida in DA:
                stack_errs = np.zeros((4,1))
                if os.path.exists( wrf_dir+'/'+istorm+'/'+Exper_names[istorm][imp][ida] ):
                    for key, array in Exper_error[istorm][imp][ida].items():
                        stack_errs = np.concatenate((stack_errs,array),axis=1)
                    allFC_error[istorm][imp][ida] = stack_errs[:,1:] #remove the first empty column
                else:
                    continue
    return allFC_error

# Bin sample: per angle per radius
def bin_Polar( allFC_error ):

    # Parameter to bin data by distance and angle
    num_r_bins = 5
    num_theta_bins = 10
    r_bins = np.linspace(0, 250, num_r_bins + 1)
    theta_bins = np.linspace(0,360, num_theta_bins + 1)
    
    # Bin data
    hist_polar = {}
    for imp in MP:
        hist_polar[imp] = {}
        for ida in DA:
            stack_errs = np.zeros((2,1))
            for istorm in Storms:
                if os.path.exists( wrf_dir+'/'+istorm+'/'+Exper_names[istorm][imp][ida] ):
                    stack_errs = np.concatenate((stack_errs,allFC_error[istorm][imp][ida][:2,:]),axis=1)
            # assemble all samples for one kind
            all_storms = stack_errs[:,1:] #remove the first empty column                                
            # bin data
            hist, r_edges, theta_edges = np.histogram2d(all_storms[0,:],np.mod(all_storms[1,:], 360), bins=[r_bins, theta_bins])   
            #hist = hist2d(all_storms[0,:],all_storms[1,:],range=[[0,250],[0,2*np.pi]],bins=50)
            hist_polar[imp][ida] = hist.T
   
    # Step 4: Plot contourf of the sample distribution
    r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])
    theta_mid = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    R, Theta = np.meshgrid(r_mid, np.deg2rad(theta_mid))
    
    return [R,Theta,hist_polar]

# ------------------------------------------------------------------------------------------------------
#            Operation: Plot PDF
# ------------------------------------------------------------------------------------------------------

def plot_polarCoor(ax, num_circles, radius_increment):

    radii = np.linspace(0.1, 250.0, 5)  # 5 concentric circles
    angles = np.linspace(0, 2 * np.pi, 360)  # 360 points for full circle

    # Draw concentric circles
    for radius in radii:
        circle = plt.Circle((0, 0), radius, color='grey', fill=False)
        ax.add_artist(circle)

    # Draw directional lines
    directions = ['E', 'N', 'W', 'S']
    angle_step = 90
    for i, direction in enumerate(directions):
        angle = np.deg2rad(angle_step * i)
        x = np.cos(angle) * radius_increment * num_circles
        y = np.sin(angle) * radius_increment * num_circles
        ax.plot([0, x], [0, y], 'k--', linewidth=1)
        ax.text(x, y, direction, fontsize=12, ha='center', va='center')

    # Set equal aspect ratio
    ax.set_aspect('equal', 'box')

    # Set limits
    max_radius = 250
    ax.set_xlim(-max_radius, max_radius)
    ax.set_ylim(-max_radius, max_radius)

    # Add grid
    ax.grid(True)

def plot_polarCoor2(ax):

    radii = np.linspace(0.1, 250.0, 5)  # 5 concentric circles
    angles = np.linspace(0, 2 * np.pi, 100)  # 360 points for full circle

    # Draw concentric circles
    for radius in radii:
        ax.plot(angles, np.full(100, radius),linestyle='')


# plot
def plot_sys_errs():

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200)
    grids = fig.add_gridspec(ncols=3,nrows=3,bottom=0.05,top=0.96,left=0.10,right=0.98,wspace=0.08,hspace=0.06)

    ax = {}
    for ida in DA:
        ax[ida] = {}
        ir = DA.index( ida )
        # Track: simulation relative to best track in a polar coordinate
        ax[ida]['ax0'] = fig.add_subplot( grids[ir,0] )
        plot_polarCoor(ax[ida]['ax0'], num_circles, radius_increment)    
        # Intensity
        ax[ida]['ax1'] = fig.add_subplot( grids[ir,1] )
        ax[ida]['ax2'] = fig.add_subplot( grids[ir,2] )

    # Customize color
    colors = {}
    colorset = {'HARVEY':'#FF13EC','JOSE':'#0D13F6','MARIA':'#FFA500','IRMA':'#DB2824'}
    alphas = np.linspace(0.2,1,8)

    # Customize marker type
    marker = {'WSM6': 'P','THO':'o'}

    # Plot simulations
    for imp in MP:
        for ida in DA:
            for istorm in Storms:
                # Plot each storm
                if not os.path.exists( wrf_dir+'/'+istorm+'/'+Exper_names[istorm][imp][ida] ):
                    continue
                distance = errs[istorm][imp][ida][0,:]
                direction = errs[istorm][imp][ida][1,:]
                # Convert polar coordinates (distance, direction) to Cartesian coordinates (x, y)
                point_angle = np.deg2rad(direction)
                px = distance * np.cos(point_angle)
                py = distance * np.sin(point_angle)
                if imp == 'WSM6':
                    ax[ida]['ax0'].plot(px,py,linestyle='',markersize=4,marker='P',color=colorset[istorm],alpha=1)
                else:
                    ax[ida]['ax0'].plot(px,py,linestyle='',markersize=4,marker='o',color=colorset[istorm],alpha=0.5)

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


# plot PDF
def plot_sys_errs_pdf( allFC_error ):

    # Bin
    returned = bin_Polar( allFC_error )
    R = returned[0]
    Theta = returned[1]
    hist_polar = returned[2]

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200)
    grids = fig.add_gridspec(ncols=3,nrows=3,bottom=0.05,top=0.96,left=0.10,right=0.98,wspace=0.08,hspace=0.06)

    ax = {}
    for ida in DA:
        ax[ida] = {}
        ir = DA.index( ida )
        # Track: simulation relative to best track in a polar coordinate
        ax[ida]['ax0'] = fig.add_subplot( grids[ir,0],projection='polar' )
        plot_polarCoor2(ax[ida]['ax0'])
        # Intensity
        ax[ida]['ax1'] = fig.add_subplot( grids[ir,1] )
        ax[ida]['ax2'] = fig.add_subplot( grids[ir,2] )

    # Customize color
    colors = {}
    colorset = {'HARVEY':'#FF13EC','JOSE':'#0D13F6','MARIA':'#FFA500','IRMA':'#DB2824'}
    alphas = np.linspace(0.2,1,8)

    # Customize marker type
    marker = {'WSM6': 'P','THO':'o'}

    # Plot simulations
    for ida in DA:
        for imp in MP:
            # Plot track distance/direction
            print(np.amax(hist_polar[imp][ida]))
            if imp == 'WSM6':
                c = ax[ida]['ax0'].pcolormesh(Theta, R, hist_polar[imp][ida], cmap='hot_r',shading='auto')        
                fig.colorbar(c, ax=ax[ida]['ax0'])
            else:
                c = ax[ida]['ax0'].pcolormesh(Theta, R, hist_polar[imp][ida], cmap='hot_r',shading='auto')        
                fig.colorbar(c, ax=ax[ida]['ax0'])


    #fig.colorbar(c, ax=ax)

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )



if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY','JOSE','IRMA','MARIA']
    DA = ['CONV','IR','IR+MW']
    MP = ['THO',] #

    # if operate over the same number of samples for all forecasts
    sameNum_sample = False
    if sameNum_sample:
        fc_run_hrs = 60

    # Track: number of concentric circles and the radius increment
    num_circles = 5
    radius_increment = 50

    # if calculate absolute error or signed error
    signed_error = True # absolute error
    
    # if plot population
    plot_ppl = True
    plot_bin = True

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

    # Calculate error
    if signed_error:
        Exper_error = calculate_signed_error()
    else:
        Exper_error = calculate_abs_error()

    # If plot population distribution or statisitcal attributes
    if plot_ppl:
        # For one kind of experiment, assmble samples from all forecasts.
        allFC_error = stack_error( Exper_error )
        # bin data if necessary
        if plot_bin:
            plot_sys_errs_pdf( allFC_error )
            #[R,Theta,hist_polar] = bin_Polar( allFC_error ) 
        # plot distribution
        #plot_sys_errs_pdf()












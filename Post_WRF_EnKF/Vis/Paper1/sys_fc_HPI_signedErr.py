
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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import linregress

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




def storm_color():
    colorset = {'HARVEY':'#FFA500','JOSE':'#0D13F6','IRMA':'#FF13EC','MARIA':'#097969',}
    return colorset

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
        best_track = UD.btk_in_duration(small_dir, istorm, Btk_start, Btk_end, hour_step=6)
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

# Stack arrays in nested dictionaries
# From MP*DA*Storm*Times to MP*DA*Storms
def stack_FCs( Exper_error ):

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

# Assemble samples of one kind (MP*DA) from different storms
def stack_Storms( allFC_error ):
    
    stack_storms = {}
    for imp in MP:
        stack_storms[imp] = {}
        for ida in DA:
            stack_errs = np.zeros((4,1)) #distance,angle,minSLP,Vmax
            for istorm in Storms:
                if os.path.exists( wrf_dir+'/'+istorm+'/'+Exper_names[istorm][imp][ida] ):
                    stack_errs = np.concatenate((stack_errs,allFC_error[istorm][imp][ida]),axis=1)
            stack_storms[imp][ida] = stack_errs[:,1:] #remove the first empty column 
    return stack_storms

# calculate the mean and median
def calc_mean_median( allFC_error ):

    mean_median = {}
    for ist in Storms:
        mean_median[ist] = {}
        for imp in MP:
            mean_median[ist][imp] = {}
            for ida in DA:
                mean_median[ist][imp][ida] = {}
                mean_median[ist][imp][ida]['mean'] = np.nanmean( allFC_error[ist][imp][ida],axis=1 )
                mean_median[ist][imp][ida]['median'] = np.nanmedian( allFC_error[ist][imp][ida],axis=1 )

    return mean_median

# For each experiment of MP*DA:
# track: bin data per distance per angle in a polar coordinate centered at the best-track position
# intensity: for each da scheme, bin data per intensity of WSM6 per intensify of THO 
def bin_allFc( allFC_error ):

    # parameters to bin track by distance and angle
    num_r_bins = 5 # bin width: 100 km
    num_theta_bins = 12 # bin width: 30 degrees
    r_bins = np.linspace(0, max_radii, num_r_bins + 1)
    theta_bins = np.linspace(0,360, num_theta_bins + 1)
    # parameters to bin MSLP
    minMSLPerr = -75
    maxMSLPerr = 75
    num_binsMSLP = 25 #maxMSLPerr-minMSLPerr
    #MSLP_edges = np.linspace(minMSLPerr,maxMSLPerr,num_binsMSLP+1)
    #MSLP_range_mid = 0.5 * (MSLP_edges[:-1] + MSLP_edges[1:])
    # parameters to bin Vmax
    minVmax_err = -60
    maxVmax_err = 60
    num_binsVmax = 25 #maxVmax_err-minVmax_err
    #Vmax_edges_cal = np.linspace(minVmax_err,maxVmax_err,num_binsVmax+1)
    #Vmax_range_mid = 0.5 * (Vmax_edges[:-1] + Vmax_edges[1:])

    # Assemble all samples for one kind (MP*DA)
    stack_storms = stack_Storms( allFC_error )

    # Marginally-bin data 
    hist_HPI = {}
    hist_HPI['track'] = {}
    hist_HPI['MSLP'] = {}
    hist_HPI['Vmax'] = {}
    for imp in MP:
        hist_HPI['track'][imp] = {}
        hist_HPI['MSLP'][imp] = {}
        hist_HPI['Vmax'][imp] = {}
        for ida in DA:
            all_storms = stack_storms[imp][ida]
            # track
            hist, r_edges, theta_edges = np.histogram2d(all_storms[0,:],np.mod(all_storms[1,:], 360), bins=[r_bins, theta_bins])  #np.mod(X, 360): Convert theta_data to the range [0, 360) 
            percentages = 100 * hist.T / len(all_storms[0,:])
            hist_HPI['track'][imp][ida] = percentages #hist.T
            # MSLP
            hist, MSLP_edges = np.histogram(all_storms[2,:],range=[minMSLPerr,maxMSLPerr],bins=num_binsMSLP)
            hist_HPI['MSLP'][imp][ida] = hist.T
            # Vmax
            hist, Vmax_edges = np.histogram(all_storms[3,:],range=[minVmax_err,maxVmax_err],bins=num_binsVmax)
            hist_HPI['Vmax'][imp][ida] = hist.T


    # meshgrid for coordinates
    r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])
    theta_mid = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    R, Theta = np.meshgrid(r_mid, np.deg2rad(theta_mid))
    hist_HPI['track']['R'] = R
    hist_HPI['track']['Theta'] = Theta
    
    hist_HPI['MSLP']['MSLP_edges'] = MSLP_edges
    hist_HPI['Vmax']['Vmax_edges'] = Vmax_edges

    # Jointly-bin intensity data
    hist_HPI['joint_MSLP'] = {}
    hist_HPI['joint_Vmax'] = {}
    for ida in DA:
        # bin minSLP
        hist, MSLP_edges, MSLP_edges = np.histogram2d(stack_storms['WSM6'][ida][2,:],stack_storms['THO'][ida][2,:],range=[[minMSLPerr,maxMSLPerr],[minMSLPerr,maxMSLPerr]],bins=[num_binsMSLP, num_binsMSLP])
        hist_HPI['joint_MSLP'][ida] = hist.T

        hist, Vmax_edges, Vmax_edges = np.histogram2d(stack_storms['WSM6'][ida][3,:],stack_storms['THO'][ida][3,:],range=[[minVmax_err,maxVmax_err],[minVmax_err,maxVmax_err]],bins=[num_binsVmax,num_binsVmax])
        hist_HPI['joint_Vmax'][ida] = hist

    return hist_HPI

# ------------------------------------------------------------------------------------------------------
#            Operation: Plot data samples (scatter)
# ------------------------------------------------------------------------------------------------------
def plot_ConcentricCircle(ax, num_circles, radius_increment):

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

# plot
def scatter_sys_errs( errs ):

    # Set up figure
    fig = plt.figure( figsize=(6.5,7.5),dpi=200)
    grids = fig.add_gridspec(ncols=3,nrows=3,bottom=0.05,top=0.96,left=0.02,right=0.98,wspace=0.08,hspace=0.06)

    ax = {}
    for ida in DA:
        ax[ida] = {}
        ir = DA.index( ida )
        # Track: simulation relative to best track in a polar coordinate
        ax[ida]['ax0'] = fig.add_subplot( grids[ir,0] )
        plot_ConcentricCircle(ax[ida]['ax0'], num_circles, radius_increment)    
        # Intensity
        ax[ida]['ax1'] = fig.add_subplot( grids[ir,1] )
        ax[ida]['ax2'] = fig.add_subplot( grids[ir,2] )

    # Customize color
    colors = {}
    colorset = {'HARVEY':'#FF13EC','JOSE':'#0D13F6','MARIA':'#FFA500','IRMA':'#DB2824'}
    alphas = np.linspace(0.2,1,8)

    # Customize marker type
    marker = {'WSM6': 'P','THO':'o'}

    # Scatter track errors
    for imp in MP:
        for ida in DA:
            for istorm in Storms: # scatter each storm
                if not os.path.exists( wrf_dir+'/'+istorm+'/'+Exper_names[istorm][imp][ida] ):
                    continue
                # Track
                distance = errs[istorm][imp][ida][0,:]
                direction = errs[istorm][imp][ida][1,:]
                # convert polar coordinates (distance, direction) to Cartesian coordinates (x, y)
                point_angle = np.deg2rad(direction)
                px = distance * np.cos(point_angle)
                py = distance * np.sin(point_angle)
                if imp == 'WSM6':
                    ax[ida]['ax0'].plot(px,py,linestyle='',markersize=4,marker='P',color=colorset[istorm],alpha=1)
                else:
                    ax[ida]['ax0'].plot(px,py,linestyle='',markersize=4,marker='o',color=colorset[istorm],alpha=0.5)

    # Scatter intensities
    for ida in DA:
        for istorm in Storms:
            # MSLP
            ax[ida]['ax1'].set_aspect('equal', 'box')
            ax[ida]['ax1'].plot(errs[istorm]['WSM6'][ida][2,:],errs[istorm]['THO'][ida][2,:],linestyle='',markersize=4,marker='o',color=colorset[istorm])
            ax[ida]['ax1'].axvline(x=0, color='k', linestyle='-',linewidth=1)
            ax[ida]['ax1'].axhline(y=0, color='k', linestyle='-',linewidth=1)
            #Vmax
            ax[ida]['ax2'].set_aspect('equal', 'box')
            ax[ida]['ax2'].plot(errs[istorm]['WSM6'][ida][3,:],errs[istorm]['THO'][ida][3,:],linestyle='',markersize=4,marker='o',color=colorset[istorm])
            ax[ida]['ax2'].axvline(x=0, color='k', linestyle='-',linewidth=1)
            ax[ida]['ax2'].axhline(y=0, color='k', linestyle='-',linewidth=1)

    # Set ticks/other attributes for intensity subplots
    for ida in DA:
        ax[ida]['ax1'].set_xlim([-75,75])
        ax[ida]['ax1'].set_ylim([-75,75])
        ax[ida]['ax2'].set_xlim([-60,60])
        ax[ida]['ax2'].set_ylim([-60,60])

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_scatter.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


# ------------------------------------------------------------------------------------------------------
#            Operation: Plot PDF
# ------------------------------------------------------------------------------------------------------

def plot_PolarCoord(ax):

    # Set up parameters
    radii = np.linspace(0.1, max_radii, num_circles)  # X concentric circles
    num_angle = 100 # number of points for full circle
    angles = np.linspace(0, 360, num_angle)  # X points for full circle

    # Draw concentric circles
    for radius in radii:
        ax.plot(angles, np.full(num_angle, radius),linestyle='')

# plot PDF
def plot_mslpVSmslp( allFC_error,hist_HPI ):

    # collect data
    stack_storms = stack_Storms( allFC_error )

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=2,top=0.93,left=0.1,hspace=0.09)

    ax = {}
    ax['wsm6_track'] = {}
    ax['tho_track'] = {}
    ax['wsm6'] = {}
    ax['tho'] = {}
    # gridspec inside gridspec
    track_grid = outer_grids[0].subgridspec(2, 3, wspace=0.05, hspace=0.14)
    its_grid = outer_grids[1].subgridspec(2, 3, wspace=0.05, hspace=0.3)

    for ida in DA:
        ir = DA.index( ida )
        # Track: simulation relative to best track in a polar coordinate
        ax['wsm6_track'][ida] = fig.add_subplot( track_grid[0,ir],projection='polar' )
        plot_PolarCoord( ax['wsm6_track'][ida] )
        ax['tho_track'][ida] = fig.add_subplot( track_grid[1,ir],projection='polar' )
        plot_PolarCoord( ax['tho_track'][ida] ) 
        # Intensity
        ax['wsm6'][ida] = fig.add_subplot( its_grid[0,ir] )
        ax['tho'][ida] = fig.add_subplot( its_grid[1,ir] )

    # Customize color
    colors = {}
    colorset = storm_color()
    alphas = {'HARVEY':1,'IRMA':0.75,'JOSE':0.5,'MARIA':0.35}

    # Plot simulation error
    R = hist_HPI['track']['R']
    Theta = hist_HPI['track']['Theta']
    mslp_rg = hist_HPI['MSLP']['MSLP_edges']
    mslp_rg_mid =  0.5 * (mslp_rg[:-1] + mslp_rg[1:])
    mslp_Xrg, mslp_Yrg = np.meshgrid(mslp_rg_mid,mslp_rg_mid)
    vmax_rg = hist_HPI['Vmax']['Vmax_edges']
    vmax_rg_mid =  0.5 * (vmax_rg[:-1] + vmax_rg[1:])
    vmax_Xrg, vmax_Yrg = np.meshgrid(vmax_rg_mid,vmax_rg_mid)

    ax_wsm6_margMSLP = {}
    ax_wsm6_margVmax = {}
    ax_tho_margMSLP = {}
    ax_tho_margVmax = {}
    for ida in DA: # column
        ax['wsm6_track'][ida].pcolormesh(Theta, R, hist_HPI['track']['WSM6'][ida], cmap='gist_heat_r',vmin=0,vmax=20)
        t_pdf = ax['tho_track'][ida].pcolormesh(Theta, R, hist_HPI['track']['THO'][ida], cmap='gist_heat_r',vmin=0,vmax=20)
        # Marginal histogram for MSLP
        # create inset axes for the marginal histograms
        ax_wsm6_margMSLP[ida] = inset_axes(ax['wsm6'][ida], width="100%", height="20%", loc='lower center')
        ax_wsm6_margMSLP[ida].set_frame_on(False)
        ax_wsm6_margVmax[ida] = inset_axes(ax['wsm6'][ida], width="20%", height="100%", loc='center left')
        ax_wsm6_margVmax[ida].set_frame_on(False)
        # plot the marginal histogram
        ax_wsm6_margVmax[ida].barh(mslp_rg_mid,hist_HPI['MSLP']['THO'][ida],height=np.diff(mslp_rg),color='gray', alpha=0.7)
        ax_wsm6_margMSLP[ida].bar(mslp_rg_mid,hist_HPI['MSLP']['WSM6'][ida],width=np.diff(mslp_rg), color='gray', alpha=0.7)
		# scatter
        for istorm in Storms:
            # scatter
            ax['wsm6'][ida].plot(allFC_error[istorm]['WSM6'][ida][2,:],allFC_error[istorm]['THO'][ida][2,:],linestyle='',markersize=4,marker='o',color=colorset[istorm],alpha=alphas[istorm])
       #i_pdf = ax['its'][ida].pcolormesh(mslp_Xrg, mslp_Yrg, hist_HPI['joint_MSLP'][ida], cmap='gist_heat_r',shading='auto')
       # Add diagonal line
        x = np.linspace(-75, 75, num=2)
        ax['wsm6'][ida].plot(x, x, color='gray',linewidth=1,alpha=0.5)
        x = np.linspace(-60, 60, num=2) 
        ax['tho'][ida].plot(x, x, color='gray',linewidth=1,alpha=0.5)

        # Marginal histogram for Vmax
        # create inset axes for the marginal histograms
        ax_tho_margMSLP[ida] = inset_axes(ax['tho'][ida], width="100%", height="20%", loc='upper center')
        ax_tho_margMSLP[ida].set_frame_on(False)
        ax_tho_margVmax[ida] = inset_axes(ax['tho'][ida], width="20%", height="100%", loc='center right')
        ax_tho_margVmax[ida].set_frame_on(False)
        # plot the marginal histogram
        ax_tho_margVmax[ida].barh(vmax_rg_mid,hist_HPI['Vmax']['THO'][ida],height=np.diff(vmax_rg),color='gray', alpha=0.7)
        ax_tho_margMSLP[ida].bar(vmax_rg_mid,hist_HPI['Vmax']['WSM6'][ida],width=np.diff(vmax_rg), color='gray', alpha=0.7)
        ax_tho_margVmax[ida].invert_xaxis()
        ax_tho_margMSLP[ida].invert_yaxis()
        # scatter
        for istorm in Storms:
            # scatter
            ax['tho'][ida].plot(allFC_error[istorm]['WSM6'][ida][3,:],allFC_error[istorm]['THO'][ida][3,:],linestyle='',markersize=4,marker='o',color=colorset[istorm],alpha=alphas[istorm])

	# Create a colorbar above the first row of subplots
    cbar_ax = fig.add_axes([0.925, 0.52, 0.03, 0.43])
    cbar = fig.colorbar(t_pdf, cax=cbar_ax, orientation='vertical')
    cbar.set_ticks([0, 5, 10, 15, 20])
    cbar.set_ticklabels(['0%', '5%', '10%', '15%', '20%'])

    # Add experiment name
    for ida in DA:
        if DA.index(ida) == 0:
            fig.text(0.23,0.97,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 1:
            fig.text(0.50,0.97,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 2:
            fig.text(0.78,0.97,ida, fontsize=12, ha='center', va='center')
	
	# Add y label
    fig.text(0.05,0.85,'WSM6: track', fontsize=10, ha='center', va='center',rotation='vertical') 
    fig.text(0.05,0.63,'THO: track', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.05,0.42,'THO: MSLP', fontsize=10, ha='center', va='center' ,rotation='vertical')
    fig.text(0.05,0.19,'THO: Vmax', fontsize=10, ha='center', va='center',rotation='vertical')

	# Add circle (legend)
    circle = plt.Circle((0.94, 0.46), 0.005, color=colorset['HARVEY'], alpha=alphas['HARVEY'],transform=fig.transFigure, clip_on=False) 
    fig.add_artist(circle)
    fig.text(0.94,0.44,'HARVEY', fontsize=10, ha='center', va='center')

    circle = plt.Circle((0.94, 0.38), 0.005, color=colorset['IRMA'], alpha=alphas['IRMA'], transform=fig.transFigure, clip_on=False) 
    fig.add_artist(circle)
    fig.text(0.94,0.36,'IRMA', fontsize=10, ha='center', va='center')

    circle = plt.Circle((0.94, 0.23), 0.005, color=colorset['JOSE'], alpha=alphas['JOSE'], transform=fig.transFigure, clip_on=False) 
    fig.add_artist(circle)
    fig.text(0.94,0.21,'JOSE', fontsize=10, ha='center', va='center')

    circle = plt.Circle((0.94, 0.16), 0.005, color=colorset['MARIA'], alpha=alphas['MARIA'], transform=fig.transFigure, clip_on=False) 
    fig.add_artist(circle)
    fig.text(0.94,0.14,'MARIA', fontsize=10, ha='center', va='center')


    # Set axis attributes
    yticks = np.linspace(0,max_radii,num_circles+1)
    for ida in DA:
        # Track
        ax['wsm6_track'][ida].set_rticks(yticks[:-1])
        ax['wsm6_track'][ida].set_rlabel_position(90)
        ax['tho_track'][ida].set_rticks(yticks[:-1])
        ax['tho_track'][ida].set_rlabel_position(90)
        if DA.index(ida) == 0:
            ax['wsm6_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['wsm6_track'][ida].set_xticklabels(['','NE', 'N', 'NW', 'W', 'SW','S/N', 'SE'])
            ax['tho_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['tho_track'][ida].set_xticklabels(['','NE','', 'NW', 'W', 'SW', 'S', 'SE'])
        elif DA.index(ida) == 1:
            ax['wsm6_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['wsm6_track'][ida].set_xticklabels(['','NE', 'N', 'NW', '','SW', 'S/N','SE'])
            ax['tho_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['tho_track'][ida].set_xticklabels(['','NE','', 'NW','', 'SW', 'S', 'SE'])
        else:
            ax['wsm6_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['wsm6_track'][ida].set_xticklabels(['E','NE', 'N', 'NW','', 'SW','S/N', 'SE'])
            ax['tho_track'][ida].set_xticks(np.deg2rad([0,45,90,135,190,225,270,315]))
            ax['tho_track'][ida].set_xticklabels(['E','NE','', 'NW', '','SW', 'S', 'SE'])

        # MSLP
        ax['wsm6'][ida].set_xlabel('WSM6: MSLP',fontsize='10')
        mslp_ticks = np.linspace(-75, 75, num=7, dtype=int)  
        ax['wsm6'][ida].set_xticks(mslp_ticks)
        ax['wsm6'][ida].set_xticklabels([str(i) for i in mslp_ticks] )
        ax['wsm6'][ida].set_xlim([-75,75])
        ax['wsm6'][ida].set_ylim([-75,75])
        # set equal aspect ratio
        ax['wsm6'][ida].set_aspect('equal', 'box')
        # set reference lines
        ax['wsm6'][ida].axvline(x=0, color='gray', linestyle='-',linewidth=1,alpha=0.5)
        ax['wsm6'][ida].axhline(y=0, color='gray', linestyle='-',linewidth=1,alpha=0.5)
        # minimize the use of y axis label
        #if DA.index( ida ) != 0:
        #    ax['wsm6'][ida].set_yticklabels([])
        # hide the marginal axes labels
        ax_wsm6_margMSLP[ida].xaxis.set_visible(False)
        ax_wsm6_margMSLP[ida].yaxis.set_ticks([])
        ax_wsm6_margVmax[ida].yaxis.set_visible(False)
        ax_wsm6_margVmax[ida].xaxis.set_ticks([])

        # Vmax
        ax['tho'][ida].set_xlabel('WSM6: Vmax', fontsize='10')
        vmax_ticks = np.linspace(-60, 60, num=7, dtype=int)
        ax['tho'][ida].set_xticks(vmax_ticks)
        ax['tho'][ida].set_xticklabels([str(i) for i in vmax_ticks])
        ax['tho'][ida].set_xlim([-60,60])
        ax['tho'][ida].set_ylim([-60,60])
        # set equal aspect ratio
        ax['tho'][ida].set_aspect('equal', 'box')
        # set reference lines
        ax['tho'][ida].axvline(x=0, color='gray', linestyle='-',linewidth=1,alpha=0.5)
        ax['tho'][ida].axhline(y=0, color='gray', linestyle='-',linewidth=1,alpha=0.5)
        # minimize the use of y axis label
        #if DA.index( ida ) != 0:
        #    ax['tho'][ida].set_yticklabels([])
        # hide the marginal axes labels
        ax_tho_margMSLP[ida].xaxis.set_visible(False)
        ax_tho_margMSLP[ida].yaxis.set_ticks([])
        ax_tho_margVmax[ida].yaxis.set_visible(False)
        ax_tho_margVmax[ida].xaxis.set_ticks([])

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_signedError_mslpVSmslp.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


# plot PDF
def plot_mslpVSvmax( allFC_error,hist_HPI ):

    mm = calc_mean_median( allFC_error )

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=2,top=0.93,left=0.1,hspace=0.09)

    ax = {}
    ax['wsm6_track'] = {}
    ax['tho_track'] = {}
    ax['wsm6'] = {}
    ax['tho'] = {}
    # gridspec inside gridspec
    track_grid = outer_grids[0].subgridspec(2, 3, wspace=0.05, hspace=0.14)
    its_grid = outer_grids[1].subgridspec(2, 3, wspace=0.05, hspace=0.3)

    for ida in DA:
        ir = DA.index( ida )
        # Track: simulation relative to best track in a polar coordinate
        ax['wsm6_track'][ida] = fig.add_subplot( track_grid[0,ir],projection='polar' )
        plot_PolarCoord( ax['wsm6_track'][ida] )
        ax['tho_track'][ida] = fig.add_subplot( track_grid[1,ir],projection='polar' )
        plot_PolarCoord( ax['tho_track'][ida] )
        # Intensity
        ax['wsm6'][ida] = fig.add_subplot( its_grid[0,ir] )
        ax['tho'][ida] = fig.add_subplot( its_grid[1,ir] )
    # Customize color
    colors = {}
    colorset = storm_color()
    alphas = {'HARVEY':1,'IRMA':0.75,'JOSE':0.5,'MARIA':0.35}

    # Plot simulation error
    R = hist_HPI['track']['R']
    Theta = hist_HPI['track']['Theta']
    mslp_rg = hist_HPI['MSLP']['MSLP_edges']
    mslp_rg_mid =  0.5 * (mslp_rg[:-1] + mslp_rg[1:])
    mslp_Xrg, mslp_Yrg = np.meshgrid(mslp_rg_mid,mslp_rg_mid)
    vmax_rg = hist_HPI['Vmax']['Vmax_edges']
    vmax_rg_mid =  0.5 * (vmax_rg[:-1] + vmax_rg[1:])
    vmax_Xrg, vmax_Yrg = np.meshgrid(vmax_rg_mid,vmax_rg_mid)

    ax_wsm6_margMSLP = {}
    ax_wsm6_margVmax = {}
    ax_tho_margMSLP = {}
    ax_tho_margVmax = {}
    for ida in DA: # column
        ax['wsm6_track'][ida].pcolormesh(Theta, R, hist_HPI['track']['WSM6'][ida], cmap='gist_heat_r',vmin=0,vmax=20)
        t_pdf = ax['tho_track'][ida].pcolormesh(Theta, R, hist_HPI['track']['THO'][ida], cmap='gist_heat_r',vmin=0,vmax=20)
        # Marginal histogram for WSM6 intensity
        # create inset axes for the marginal histograms
        ax_wsm6_margMSLP[ida] = inset_axes(ax['wsm6'][ida], width="100%", height="20%", loc='lower center')
        ax_wsm6_margMSLP[ida].set_frame_on(False)
        ax_wsm6_margVmax[ida] = inset_axes(ax['wsm6'][ida], width="20%", height="100%", loc='center left')
        ax_wsm6_margVmax[ida].set_frame_on(False)
        # plot the marginal histogram
        ax_wsm6_margVmax[ida].barh(vmax_rg_mid,hist_HPI['Vmax']['WSM6'][ida],height=np.diff(vmax_rg),color='gray', alpha=0.7)
        ax_wsm6_margMSLP[ida].bar(mslp_rg_mid,hist_HPI['MSLP']['WSM6'][ida],width=np.diff(mslp_rg), color='gray', alpha=0.7)
        # scatter
        stack_sts = np.zeros((4,1))
        for ist in Storms:
            # scatter the mean and median
            ax['wsm6'][ida].scatter(mm[ist]['WSM6'][ida]['mean'][2],mm[ist]['WSM6'][ida]['mean'][3],linewidths=1.5,marker='o',edgecolor=colorset[ist],facecolors='none',s=25,alpha=0.6)
            ax['wsm6'][ida].scatter(mm[ist]['WSM6'][ida]['median'][2],mm[ist]['WSM6'][ida]['median'][3],linewidths=1.5,marker='+',facecolor=colorset[ist],s=25)
            #ax['wsm6'][ida].plot(allFC_error[istorm]['WSM6'][ida][2,:],allFC_error[istorm]['WSM6'][ida][3,:],linestyle='',markersize=4,marker='o',color=colorset[istorm],alpha=alphas[istorm])
            # stack
            stack_sts = np.concatenate((stack_sts,allFC_error[ist]['WSM6'][ida]),axis=1)

        stack_wsm6 = stack_sts[:,1:] #remove the first empty column
        # remove data point (column) if any element is nan
        nan_mask = np.isnan(stack_wsm6).any(axis=0) # identify nan elements using a boolean mask
        clean_wsm6 = stack_wsm6[:,~nan_mask] # remove NaN elements by inverting the mask
        # Fit the least squares line
        slope, intercept, r_value, p_value, std_err = linregress(clean_wsm6[2,:], clean_wsm6[3,:]) # linear regression
        # Generate y values for the fitted line
        x = np.linspace(-75, 75, num=2)
        y_fit = slope * x + intercept
        ax['wsm6'][ida].plot(x, y_fit, color='black',linewidth=1,alpha=0.5)
        #ax['wsm6'][ida].text(-50,50,"y="+f"{slope:.2f}"+"*x+"+f"({intercept:.2f})",fontsize=8,horizontalalignment='left')
        ax['wsm6'][ida].text(10,20,"$\mathregular{R^2}$: "+f"{r_value**2:.2f}")
        
        # calculate mean and median over all storms
        mean_all = np.mean(clean_wsm6,axis=1)
        median_all = np.median(clean_wsm6,axis=1)
        ax['wsm6'][ida].scatter(mean_all[2],mean_all[3],marker='o',edgecolor='black',facecolors='none',linewidths=1.5,s=25,alpha=0.6)
        ax['wsm6'][ida].scatter(median_all[2],median_all[3],marker='+',facecolors='black',linewidths=1.5,s=25)

        #i_pdf = ax['its'][ida].pcolormesh(mslp_Xrg, mslp_Yrg, hist_HPI['joint_MSLP'][ida], cmap='gist_heat_r',shading='auto')
        # Add diagonal line
        #x = np.linspace(-75, 75, num=2)
        #ax['wsm6'][ida].plot(x, x, color='gray',linewidth=1,alpha=0.5)
        #x = np.linspace(-75, 75, num=2)
        #ax['tho'][ida].plot(x, x, color='gray',linewidth=1,alpha=0.5)

        # Marginal histogram for Vmax
        # create inset axes for the marginal histograms
        ax_tho_margMSLP[ida] = inset_axes(ax['tho'][ida], width="100%", height="20%", loc='lower center')
        ax_tho_margMSLP[ida].set_frame_on(False)
        ax_tho_margVmax[ida] = inset_axes(ax['tho'][ida], width="20%", height="100%", loc='center left')
        ax_tho_margVmax[ida].set_frame_on(False)
        # plot the marginal histogram
        ax_tho_margVmax[ida].barh(vmax_rg_mid,hist_HPI['Vmax']['THO'][ida],height=np.diff(vmax_rg),color='gray', alpha=0.7)
        ax_tho_margMSLP[ida].bar(mslp_rg_mid,hist_HPI['MSLP']['THO'][ida],width=np.diff(mslp_rg), color='gray', alpha=0.7)
        #ax_tho_margVmax[ida].invert_xaxis()
        #ax_tho_margMSLP[ida].invert_yaxis()
        # scatter
        stack_sts = np.zeros((4,1))
        for ist in Storms:
            # scatter the mean and median
            ax['tho'][ida].scatter(mm[ist]['THO'][ida]['mean'][2],mm[ist]['THO'][ida]['mean'][3],marker='o',linewidths=1.5,edgecolor=colorset[ist],facecolors='none',s=25,alpha=0.6)
            ax['tho'][ida].scatter(mm[ist]['THO'][ida]['median'][2],mm[ist]['THO'][ida]['median'][3],linewidths=1.5,marker='+',facecolor=colorset[ist],s=25)
            #ax['tho'][ida].plot(allFC_error[istorm]['THO'][ida][2,:],allFC_error[istorm]['THO'][ida][3,:],linestyle='',markersize=4,marker='o',color=colorset[istorm],alpha=alphas[istorm])
            # stack
            stack_sts = np.concatenate((stack_sts,allFC_error[ist]['THO'][ida]),axis=1)
        
        stack_tho = stack_sts[:,1:] #remove the first empty column
        nan_mask = np.isnan(stack_tho).any(axis=0) # identify nan elements using a boolean mask
        clean_tho = stack_tho[:,~nan_mask] # remove NaN elements by inverting the mask
        # Fit the least squares line
        slope, intercept, r_value, p_value, std_err = linregress(clean_tho[2,:], clean_tho[3,:])  # 1 means fitting a line (degree 1 polynomial)
        # Generate y values for the fitted line
        x = np.linspace(-75, 75, num=2)
        y_fit = slope * x + intercept
        ax['tho'][ida].plot(x, y_fit, color='black',linewidth=1,alpha=0.5)
        #ax['tho'][ida].text(-50,50,"y="+f"{slope:.2f}"+"*x+"+f"({intercept:.2f})",fontsize=8,horizontalalignment='left')
        ax['tho'][ida].text(10,20,"$\mathregular{R^2}$: "+f"{r_value**2:.2f}")

        # calculate mean and median over all storms
        mean_all = np.mean(clean_tho,axis=1)
        median_all = np.median(clean_tho,axis=1)
        ax['tho'][ida].scatter(mean_all[2],mean_all[3],marker='o',edgecolor='black',facecolors='none',linewidths=1.5,s=25,alpha=0.6)
        ax['tho'][ida].scatter(median_all[2],median_all[3],marker='+',facecolors='black',linewidths=1.5,s=25)

    # Create a colorbar above the first row of subplots
    cbar_ax = fig.add_axes([0.925, 0.52, 0.03, 0.43])
    cbar = fig.colorbar(t_pdf, cax=cbar_ax, orientation='vertical',extend='max')
    cbar.set_ticks([0, 5, 10, 15, 20])
    cbar.set_ticklabels(['0%', '5%', '10%', '15%', '20%'])

    # Add legend
    scatter_DA = ax['wsm6'][DA[0]].collections
    legend_DA = ax['wsm6'][DA[0]].legend([scatter_DA[i] for i in [-2,-1]],['Mean','Median'],fontsize='7',loc='upper center')
    scatter_DA = ax['tho'][DA[0]].collections
    legend_DA = ax['tho'][DA[0]].legend([scatter_DA[i] for i in [-2,-1]],['Mean','Median'],fontsize='7',loc='upper center')

    # Add experiment name
    for ida in DA:
        if DA.index(ida) == 0:
            fig.text(0.23,0.97,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 1:
            fig.text(0.50,0.97,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 2:
            fig.text(0.78,0.97,ida, fontsize=12, ha='center', va='center')

    # Add y label
    fig.text(0.05,0.85,'WSM6: Track Bias (km)', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.05,0.63,'THO: Track Bias (km)', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.03,0.42,'WSM6', fontsize=10, ha='center', va='center' ,rotation='vertical')
    fig.text(0.06,0.42,'Vmax Bias (m $\mathregular{s^{-1}}$)', fontsize=8, ha='center', va='center' ,rotation='vertical')
    fig.text(0.03,0.19,'THO', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.06,0.19,'Vmax Bias (m $\mathregular{s^{-1}}$)', fontsize=8, ha='center', va='center' ,rotation='vertical')

    # Add circle (legend)
    circle = plt.Circle((0.94, 0.46), 0.005, color=colorset['HARVEY'], alpha=alphas['HARVEY'], transform=fig.transFigure, clip_on=False)
    fig.add_artist(circle)
    fig.text(0.94,0.44,'HARVEY', fontsize=10, ha='center', va='center')

    circle = plt.Circle((0.94, 0.38), 0.005, color=colorset['IRMA'], alpha=alphas['IRMA'], transform=fig.transFigure, clip_on=False)
    fig.add_artist(circle)
    fig.text(0.94,0.36,'IRMA', fontsize=10, ha='center', va='center')

    circle = plt.Circle((0.94, 0.30), 0.005, color=colorset['JOSE'], alpha=alphas['JOSE'], transform=fig.transFigure, clip_on=False)
    fig.add_artist(circle)
    fig.text(0.94,0.28,'JOSE', fontsize=10, ha='center', va='center')

    circle = plt.Circle((0.94, 0.22), 0.005, color=colorset['MARIA'], alpha=alphas['MARIA'], transform=fig.transFigure, clip_on=False)
    fig.add_artist(circle)
    fig.text(0.94,0.19,'MARIA', fontsize=10, ha='center', va='center')

    circle = plt.Circle((0.94, 0.14), 0.005, color='black', alpha=1, transform=fig.transFigure, clip_on=False)
    fig.add_artist(circle)
    fig.text(0.94,0.12,'All', fontsize=10, ha='center', va='center')


    # Set axis attributes
    yticks = np.linspace(0,max_radii,num_circles+1)
    for ida in DA:
        # Track
        ax['wsm6_track'][ida].set_rticks(yticks[:-1])
        ax['wsm6_track'][ida].set_rlabel_position(90)
        ax['tho_track'][ida].set_rticks(yticks[:-1])
        ax['tho_track'][ida].set_rlabel_position(90)
        if DA.index(ida) == 0:
            ax['wsm6_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['wsm6_track'][ida].set_xticklabels(['','NE', 'N', 'NW', 'W', 'SW','S/N', 'SE'])
            ax['tho_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['tho_track'][ida].set_xticklabels(['','NE','', 'NW', 'W', 'SW', 'S', 'SE'])
        elif DA.index(ida) == 1:
            ax['wsm6_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['wsm6_track'][ida].set_xticklabels(['','NE', 'N', 'NW', '','SW', 'S/N','SE'])
            ax['tho_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['tho_track'][ida].set_xticklabels(['','NE','', 'NW','', 'SW', 'S', 'SE'])
        else:
            ax['wsm6_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['wsm6_track'][ida].set_xticklabels(['E','NE', 'N', 'NW','', 'SW','S/N', 'SE'])
            ax['tho_track'][ida].set_xticks(np.deg2rad([0,45,90,135,190,225,270,315]))
            ax['tho_track'][ida].set_xticklabels(['E','NE','', 'NW', '','SW', 'S', 'SE'])

        # WSM6 intensity
        horizontalalignment='center',ax['wsm6'][ida].set_xlabel('MSLP Bias (hPa)',fontsize='9')
        mslp_ticks = np.linspace(-75, 75, num=7, dtype=int)
        ax['wsm6'][ida].set_xticks(mslp_ticks)
        ax['wsm6'][ida].set_xticklabels([str(i) for i in mslp_ticks] )
        ax['wsm6'][ida].set_xlim([-75,75])
        ax['wsm6'][ida].set_ylim([-75,75])
        # set equal aspect ratio
        ax['wsm6'][ida].set_aspect('equal', 'box')
        # set reference lines
        ax['wsm6'][ida].axvline(x=0, color='gray', linestyle='-',linewidth=0.5,alpha=0.5)
        ax['wsm6'][ida].axhline(y=0, color='gray', linestyle='-',linewidth=0.5,alpha=0.5)
        # minimize the use of y axis label
        #if DA.index( ida ) != 0:
        #    ax['wsm6'][ida].set_yticklabels([])
        # hide the marginal axes labels
        ax_wsm6_margMSLP[ida].xaxis.set_visible(False)
        ax_wsm6_margMSLP[ida].yaxis.set_ticks([])
        ax_wsm6_margVmax[ida].yaxis.set_visible(False)
        ax_wsm6_margVmax[ida].xaxis.set_ticks([])

        # THO intensity
        ax['tho'][ida].set_xlabel('MSLP Bias (hPa)', fontsize='9')
        vmax_ticks = np.linspace(-75, 75, num=7, dtype=int)
        ax['tho'][ida].set_xticks(vmax_ticks)
        ax['tho'][ida].set_xticklabels([str(i) for i in vmax_ticks])
        ax['tho'][ida].set_xlim([-75,75])
        ax['tho'][ida].set_ylim([-75,75])
        # set equal aspect ratio
        ax['tho'][ida].set_aspect('equal', 'box')
        # set reference lines
        ax['tho'][ida].axvline(x=0, color='gray', linestyle='-',linewidth=1,alpha=0.5)
        ax['tho'][ida].axhline(y=0, color='gray', linestyle='-',linewidth=1,alpha=0.5)
        # minimize the use of y axis label
        #if DA.index( ida ) != 0:
        #    ax['tho'][ida].set_yticklabels([])
        # hide the marginal axes labels
        ax_tho_margMSLP[ida].xaxis.set_visible(False)
        ax_tho_margMSLP[ida].yaxis.set_ticks([])
        ax_tho_margVmax[ida].yaxis.set_visible(False)
        ax_tho_margVmax[ida].xaxis.set_ticks([])

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_signedError_mslpVSvmax.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


def plot_mslpVSvmax_track( allFC_error,hist_HPI ):

    # Set up figure
    fig = plt.figure( figsize=(6.5,4.5),dpi=200) # standard: 6.5,8.5
    track_grids = fig.add_gridspec(ncols=3,nrows=2,top=0.88,left=0.08,hspace=0.14,wspace=0.05)

    ax = {}
    ax['wsm6_track'] = {}
    ax['tho_track'] = {}
    # gridspec inside gridspec
    for ida in DA:
        ir = DA.index( ida )
        # Track: simulation relative to best track in a polar coordinate
        ax['wsm6_track'][ida] = fig.add_subplot( track_grids[0,ir],projection='polar' )
        plot_PolarCoord( ax['wsm6_track'][ida] )
        ax['tho_track'][ida] = fig.add_subplot( track_grids[1,ir],projection='polar' )
        plot_PolarCoord( ax['tho_track'][ida] )
    
    # Customize color
    colors = {}
    colorset = storm_color()
    alphas = {'HARVEY':1,'IRMA':0.75,'JOSE':0.5,'MARIA':0.35}

    # Plot simulation error
    R = hist_HPI['track']['R']
    Theta = hist_HPI['track']['Theta']

    for ida in DA: # column
        ax['wsm6_track'][ida].pcolormesh(Theta, R, hist_HPI['track']['WSM6'][ida], cmap='gist_heat_r',vmin=0,vmax=20)
        t_pdf = ax['tho_track'][ida].pcolormesh(Theta, R, hist_HPI['track']['THO'][ida], cmap='gist_heat_r',vmin=0,vmax=20)

    # Create a colorbar above the first row of subplots
    cbar_ax = fig.add_axes([0.925, 0.08, 0.03, 0.83])
    cbar = fig.colorbar(t_pdf, cax=cbar_ax, orientation='vertical',extend='max')
    cbar.set_ticks([0, 5, 10, 15, 20])
    cbar.set_ticklabels(['0%', '5%', '10%', '15%', '20%'])

    # Add experiment name
    for ida in DA:
        if DA.index(ida) == 0:
            fig.text(0.21,0.97,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 1:
            fig.text(0.49,0.97,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 2:
            fig.text(0.78,0.97,ida, fontsize=12, ha='center', va='center')

    # Add y label
    fig.text(0.03,0.75,'WSM6: Track Bias (km)', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.03,0.30,'THO: Track Bias (km)', fontsize=10, ha='center', va='center',rotation='vertical')

    # Set axis attributes
    yticks = np.linspace(0,max_radii,num_circles+1)
    for ida in DA:
        # Track
        ax['wsm6_track'][ida].set_rticks(yticks[:-1])
        ax['wsm6_track'][ida].set_yticklabels("")
        #ax['wsm6_track'][ida].set_rlabel_position(90)
        ax['tho_track'][ida].set_rticks(yticks[:-1])
        ax['tho_track'][ida].set_yticklabels("")
        #ax['tho_track'][ida].set_rlabel_position(90)
        if DA.index(ida) == 0:
            ax['wsm6_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['wsm6_track'][ida].set_xticklabels(['','NE', 'N', 'NW', 'W', 'SW','S/N', 'SE'])
            ax['tho_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['tho_track'][ida].set_xticklabels(['','NE','', 'NW', 'W', 'SW', 'S', 'SE'])
        elif DA.index(ida) == 1:
            ax['wsm6_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['wsm6_track'][ida].set_xticklabels(['','NE', 'N', 'NW', '','SW', 'S/N','SE'])
            ax['tho_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['tho_track'][ida].set_xticklabels(['','NE','', 'NW','', 'SW', 'S', 'SE'])
        else:
            ax['wsm6_track'][ida].set_xticks(np.deg2rad([0,45,90,135,180,225,270,315]))
            ax['wsm6_track'][ida].set_xticklabels(['E','NE', 'N', 'NW','', 'SW','S/N', 'SE'])
            ax['tho_track'][ida].set_xticks(np.deg2rad([0,45,90,135,190,225,270,315]))
            ax['tho_track'][ida].set_xticklabels(['E','NE','', 'NW', '','SW', 'S', 'SE'])

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_signedError_track.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

def plot_mslpVSvmax_intensity( allFC_error,hist_HPI ):

    mm = calc_mean_median( allFC_error )

    # Set up figure
    fig = plt.figure( figsize=(6.5,4.5),dpi=200) # standard: 6.5,8.5
    its_grids = fig.add_gridspec(ncols=3,nrows=2,top=0.93,left=0.1,hspace=0.05)

    ax = {}
    ax['wsm6'] = {}
    ax['tho'] = {}
    # gridspec inside gridspec
    for ida in DA:
        ir = DA.index( ida )
        # Intensity
        ax['wsm6'][ida] = fig.add_subplot( its_grids[0,ir] )
        ax['tho'][ida] = fig.add_subplot( its_grids[1,ir] )
    # Customize color
    colors = {}
    colorset = storm_color()
    alphas = {'HARVEY':1,'IRMA':0.75,'JOSE':0.5,'MARIA':0.35}

    # Plot simulation error
    mslp_rg = hist_HPI['MSLP']['MSLP_edges']
    mslp_rg_mid =  0.5 * (mslp_rg[:-1] + mslp_rg[1:])
    mslp_Xrg, mslp_Yrg = np.meshgrid(mslp_rg_mid,mslp_rg_mid)
    vmax_rg = hist_HPI['Vmax']['Vmax_edges']
    vmax_rg_mid =  0.5 * (vmax_rg[:-1] + vmax_rg[1:])
    vmax_Xrg, vmax_Yrg = np.meshgrid(vmax_rg_mid,vmax_rg_mid)

    ax_wsm6_margMSLP = {}
    ax_wsm6_margVmax = {}
    ax_tho_margMSLP = {}
    ax_tho_margVmax = {}
    for ida in DA: # column
        # Marginal histogram for WSM6 intensity
        # create inset axes for the marginal histograms
        ax_wsm6_margMSLP[ida] = inset_axes(ax['wsm6'][ida], width="100%", height="20%", loc='lower center')
        ax_wsm6_margMSLP[ida].set_frame_on(False)
        ax_wsm6_margVmax[ida] = inset_axes(ax['wsm6'][ida], width="20%", height="100%", loc='center left')
        ax_wsm6_margVmax[ida].set_frame_on(False)
        # plot the marginal histogram
        ax_wsm6_margVmax[ida].barh(vmax_rg_mid,hist_HPI['Vmax']['WSM6'][ida],height=np.diff(vmax_rg),color='gray', alpha=0.7)
        ax_wsm6_margMSLP[ida].bar(mslp_rg_mid,hist_HPI['MSLP']['WSM6'][ida],width=np.diff(mslp_rg), color='gray', alpha=0.7)
        # scatter
        stack_sts = np.zeros((4,1))
        for ist in Storms:
            # scatter the mean and median
            ax['wsm6'][ida].scatter(mm[ist]['WSM6'][ida]['mean'][2],mm[ist]['WSM6'][ida]['mean'][3],linewidths=1.5,marker='o',edgecolor=colorset[ist],facecolors='none',s=25,alpha=0.6)
            ax['wsm6'][ida].scatter(mm[ist]['WSM6'][ida]['median'][2],mm[ist]['WSM6'][ida]['median'][3],linewidths=1.5,marker='+',facecolor=colorset[ist],s=25)
            #ax['wsm6'][ida].plot(allFC_error[istorm]['WSM6'][ida][2,:],allFC_error[istorm]['WSM6'][ida][3,:],linestyle='',markersize=4,marker='o',color=colorset[istorm],alpha=alphas[istorm])
            # stack
            stack_sts = np.concatenate((stack_sts,allFC_error[ist]['WSM6'][ida]),axis=1)

        stack_wsm6 = stack_sts[:,1:] #remove the first empty column
        # remove data point (column) if any element is nan
        nan_mask = np.isnan(stack_wsm6).any(axis=0) # identify nan elements using a boolean mask
        clean_wsm6 = stack_wsm6[:,~nan_mask] # remove NaN elements by inverting the mask
        # Fit the least squares line
        slope, intercept, r_value, p_value, std_err = linregress(clean_wsm6[2,:], clean_wsm6[3,:]) # linear regression
        # Generate y values for the fitted line
        x = np.linspace(-75, 75, num=2)
        y_fit = slope * x + intercept
        ax['wsm6'][ida].plot(x, y_fit, color='black',linewidth=1,alpha=0.5)
        #ax['wsm6'][ida].text(-50,50,"y="+f"{slope:.2f}"+"*x+"+f"({intercept:.2f})",fontsize=8,horizontalalignment='left')
        ax['wsm6'][ida].text(10,20,"$\mathregular{R^2}$: "+f"{r_value**2:.2f}")

        # calculate mean and median over all storms
        mean_all = np.mean(clean_wsm6,axis=1)
        median_all = np.median(clean_wsm6,axis=1)
        ax['wsm6'][ida].scatter(mean_all[2],mean_all[3],marker='o',edgecolor='black',facecolors='none',linewidths=1.5,s=25,alpha=0.6)
        ax['wsm6'][ida].scatter(median_all[2],median_all[3],marker='+',facecolors='black',linewidths=1.5,s=25)
        # Marginal histogram for Vmax
        # create inset axes for the marginal histograms
        ax_tho_margMSLP[ida] = inset_axes(ax['tho'][ida], width="100%", height="20%", loc='lower center')
        ax_tho_margMSLP[ida].set_frame_on(False)
        ax_tho_margVmax[ida] = inset_axes(ax['tho'][ida], width="20%", height="100%", loc='center left')
        ax_tho_margVmax[ida].set_frame_on(False)
        # plot the marginal histogram
        ax_tho_margVmax[ida].barh(vmax_rg_mid,hist_HPI['Vmax']['THO'][ida],height=np.diff(vmax_rg),color='gray', alpha=0.7)
        ax_tho_margMSLP[ida].bar(mslp_rg_mid,hist_HPI['MSLP']['THO'][ida],width=np.diff(mslp_rg), color='gray', alpha=0.7)
        #ax_tho_margVmax[ida].invert_xaxis()
        #ax_tho_margMSLP[ida].invert_yaxis()
        # scatter
        stack_sts = np.zeros((4,1))
        for ist in Storms:
            # scatter the mean and median
            ax['tho'][ida].scatter(mm[ist]['THO'][ida]['mean'][2],mm[ist]['THO'][ida]['mean'][3],marker='o',linewidths=1.5,edgecolor=colorset[ist],facecolors='none',s=25,alpha=0.6)
            ax['tho'][ida].scatter(mm[ist]['THO'][ida]['median'][2],mm[ist]['THO'][ida]['median'][3],linewidths=1.5,marker='+',facecolor=colorset[ist],s=25)
            #ax['tho'][ida].plot(allFC_error[istorm]['THO'][ida][2,:],allFC_error[istorm]['THO'][ida][3,:],linestyle='',markersize=4,marker='o',color=colorset[istorm],alpha=alphas[istorm])
            # stack
            stack_sts = np.concatenate((stack_sts,allFC_error[ist]['THO'][ida]),axis=1)

        stack_tho = stack_sts[:,1:] #remove the first empty column
        nan_mask = np.isnan(stack_tho).any(axis=0) # identify nan elements using a boolean mask
        clean_tho = stack_tho[:,~nan_mask] # remove NaN elements by inverting the mask
        # Fit the least squares line
        slope, intercept, r_value, p_value, std_err = linregress(clean_tho[2,:], clean_tho[3,:])  # 1 means fitting a line (degree 1 polynomial)
        # Generate y values for the fitted line
        x = np.linspace(-75, 75, num=2)
        y_fit = slope * x + intercept
        ax['tho'][ida].plot(x, y_fit, color='black',linewidth=1,alpha=0.5)
        #ax['tho'][ida].text(-50,50,"y="+f"{slope:.2f}"+"*x+"+f"({intercept:.2f})",fontsize=8,horizontalalignment='left')
        ax['tho'][ida].text(10,20,"$\mathregular{R^2}$: "+f"{r_value**2:.2f}")

        # calculate mean and median over all storms
        mean_all = np.mean(clean_tho,axis=1)
        median_all = np.median(clean_tho,axis=1)
        ax['tho'][ida].scatter(mean_all[2],mean_all[3],marker='o',edgecolor='black',facecolors='none',linewidths=1.5,s=25,alpha=0.6)
        ax['tho'][ida].scatter(median_all[2],median_all[3],marker='+',facecolors='black',linewidths=1.5,s=25)

    # Add legend
    scatter_DA = ax['wsm6'][DA[0]].collections
    legend_DA = ax['wsm6'][DA[0]].legend([scatter_DA[i] for i in [-2,-1]],['Mean','Median'],fontsize='7',loc='upper center')
    scatter_DA = ax['tho'][DA[0]].collections
    legend_DA = ax['tho'][DA[0]].legend([scatter_DA[i] for i in [-2,-1]],['Mean','Median'],fontsize='7',loc='upper center')


    # Add experiment name
    for ida in DA:
        if DA.index(ida) == 0:
            fig.text(0.21,0.93,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 1:
            fig.text(0.50,0.93,ida, fontsize=12, ha='center', va='center')
        elif DA.index(ida) == 2:
            fig.text(0.78,0.93,ida, fontsize=12, ha='center', va='center')

    # Add y label
    fig.text(0.02,0.75,'WSM6', fontsize=10, ha='center', va='center' ,rotation='vertical')
    fig.text(0.05,0.75,'Vmax Bias (m $\mathregular{s^{-1}}$)', fontsize=8, ha='center', va='center' ,rotation='vertical')
    fig.text(0.02,0.30,'THO', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.05,0.30,'Vmax Bias (m $\mathregular{s^{-1}}$)', fontsize=8, ha='center', va='center' ,rotation='vertical')

    # Add circle (legend)
    circle = plt.Circle((0.94, 0.85), 0.006, color=colorset['HARVEY'], alpha=0.5, transform=fig.transFigure, clip_on=False)
    fig.add_artist(circle)
    fig.text(0.945,0.82,'HARVEY', fontsize=10, ha='center', va='center')

    circle = plt.Circle((0.94, 0.70), 0.006, color=colorset['IRMA'], alpha=0.5, transform=fig.transFigure, clip_on=False)
    fig.add_artist(circle)
    fig.text(0.94,0.67,'IRMA', fontsize=10, ha='center', va='center')

    circle = plt.Circle((0.94, 0.55), 0.006, color=colorset['JOSE'], alpha=0.5, transform=fig.transFigure, clip_on=False)
    fig.add_artist(circle)
    fig.text(0.94,0.52,'JOSE', fontsize=10, ha='center', va='center')

    circle = plt.Circle((0.94, 0.40), 0.006, color=colorset['MARIA'], alpha=0.5, transform=fig.transFigure, clip_on=False)
    fig.add_artist(circle)
    fig.text(0.94,0.37,'MARIA', fontsize=10, ha='center', va='center')

    circle = plt.Circle((0.94, 0.25), 0.006, color='black', alpha=1, transform=fig.transFigure, clip_on=False)
    fig.add_artist(circle)
    fig.text(0.94,0.22,'All', fontsize=10, ha='center', va='center')

    # Set axis attributes
    yticks = np.linspace(0,max_radii,num_circles+1)
    for ida in DA:
        # WSM6 intensity
        horizontalalignment='center',ax['wsm6'][ida].set_xlabel('MSLP Bias (hPa)',fontsize='9')
        mslp_ticks = np.linspace(-75, 75, num=7, dtype=int)
        ax['wsm6'][ida].set_xticks(mslp_ticks)
        ax['wsm6'][ida].set_xticklabels([str(i) for i in mslp_ticks] )
        ax['wsm6'][ida].set_xlim([-75,75])
        ax['wsm6'][ida].set_ylim([-75,75])
        # set equal aspect ratio
        ax['wsm6'][ida].set_aspect('equal', 'box')
        # set reference lines
        ax['wsm6'][ida].axvline(x=0, color='gray', linestyle='-',linewidth=0.5,alpha=0.5)
        ax['wsm6'][ida].axhline(y=0, color='gray', linestyle='-',linewidth=0.5,alpha=0.5)
        # minimize the use of y axis label
        #if DA.index( ida ) != 0:
        #    ax['wsm6'][ida].set_yticklabels([])
        # hide the marginal axes labels
        ax_wsm6_margMSLP[ida].xaxis.set_visible(False)
        ax_wsm6_margMSLP[ida].yaxis.set_ticks([])
        ax_wsm6_margVmax[ida].yaxis.set_visible(False)
        ax_wsm6_margVmax[ida].xaxis.set_ticks([])
        # THO intensity
        ax['tho'][ida].set_xlabel('MSLP Bias (hPa)', fontsize='9')
        vmax_ticks = np.linspace(-75, 75, num=7, dtype=int)
        ax['tho'][ida].set_xticks(vmax_ticks)
        ax['tho'][ida].set_xticklabels([str(i) for i in vmax_ticks])
        ax['tho'][ida].set_xlim([-75,75])
        ax['tho'][ida].set_ylim([-75,75])
        # set equal aspect ratio
        ax['tho'][ida].set_aspect('equal', 'box')
        # set reference lines
        ax['tho'][ida].axvline(x=0, color='gray', linestyle='-',linewidth=1,alpha=0.5)
        ax['tho'][ida].axhline(y=0, color='gray', linestyle='-',linewidth=1,alpha=0.5)
        # minimize the use of y axis label
        #if DA.index( ida ) != 0:
        #    ax['tho'][ida].set_yticklabels([])
        # hide the marginal axes labels
        ax_tho_margMSLP[ida].xaxis.set_visible(False)
        ax_tho_margMSLP[ida].yaxis.set_ticks([])
        ax_tho_margVmax[ida].yaxis.set_visible(False)
        ax_tho_margVmax[ida].xaxis.set_ticks([])

    # Save figure
    des_name = small_dir+'SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_signedError_intensity.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/Clean_results/'

    #--------Configuration------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    DA = ['CONV','IR','MW']
    MP = ['WSM6','THO'] #

    # if operate over the same number of samples for all forecasts
    sameNum_sample = False
    if sameNum_sample:
        fc_run_hrs = 60

    # Track: number of concentric circles and the radius increment
    max_radii = 500.0
    num_circles = 5
    radius_increment = max_radii/num_circles

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
        allFC_error = stack_FCs( Exper_error )
        # bin data if necessary
        if plot_bin:
            hist_HPI = bin_allFc( allFC_error ) 
            #plot_mslpVSmslp( allFC_error,hist_HPI )
            #plot_mslpVSvmax( allFC_error,hist_HPI )
            plot_mslpVSvmax_track( allFC_error,hist_HPI )
            plot_mslpVSvmax_intensity( allFC_error,hist_HPI )
        else:
            pass
            #scatter_sys_errs( allFC_error ) 












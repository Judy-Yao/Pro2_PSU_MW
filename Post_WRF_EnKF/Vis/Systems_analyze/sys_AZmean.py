#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import pyproj
import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import netCDF4 as nc
#from wrf import getvar, interplevel
from scipy import interpolate
import scipy as sp
import scipy.ndimage
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import math
import warnings
import pickle
from matplotlib.colors import LinearSegmentedColormap

import Util_data as UD

def load_AZmean_oneExper(istorm,imp,ida,Exper_names,DAtimes):

    d_az = {}
    for DAtime in DAtimes[istorm]:
        # Read coordinate-transformed data
        saved_des = small_dir+istorm+'/'+Exper_names[istorm][imp][ida]+'/Data_analyze/Model/Azmean_'+var_name+'_'+DAtime+'.pickle'
        with open(saved_des, 'rb') as handle:
            d_az[DAtime] = pickle.load( handle )
    return d_az

def plot_0thCycle_expers( storm_mean ):

    fig, ax=plt.subplots(1, 3, sharex='all', sharey='all', linewidth=0.5, figsize=(22,7), dpi=400)

    # x axis: radius 
    xv = np.linspace(0,400,100)
    x_axis_rg = range(len(xv))
    f_xinterp = interpolate.interp1d( xv, x_axis_rg)
    # y axis: model vertical coordinate
    if use_pressure:
        y_range = np.arange(900,50,50)
    else:
        y_range = np.arange(0,31,1)
    y_axis_rg = range(len(y_range))
    f_yinterp = interpolate.interp1d( y_range, y_axis_rg)
    yv = f_yinterp( storm_mean[MP[0]][DA[0]]['cycle0']['ver_coor'] )
    # Make a mesh grid
    xcoor, ycoor = np.meshgrid( x_axis_rg, yv )

    # Set colorbar
    if 'radial_wind' in var_name:
        bounds = np.linspace(-12,12,9)
        cmap = 'PiYG_r'
    elif 'tangential_wind' in var_name:
        color_intervals = list(np.linspace(0,50,6))
        color_intervals.insert(0,-10.0)
        exist_cmap = plt.cm.jet
        colors = exist_cmap(np.linspace(0,1,len(color_intervals)))
        cmap = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals))
        bounds = color_intervals
    elif 'hor_wind' in var_name:
        # tangential_wind
        color_intervals = list(np.linspace(0,21,8))
        color_intervals.insert(0,-5.0)
        exist_cmap = plt.cm.rainbow
        colors = exist_cmap(np.linspace(0,1,len(color_intervals)))
        cmap = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals))
        tan_bounds = color_intervals
        # radial wind
        wind_level = [-2,-1,1,2]
        colors = [(0, 0, 1), (1, 1, 1), (1, 0, 0)]  # Blue, White, Red
        cmap_name = 'blue_black_red'
        wind_level_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=len(wind_level))

    elif 'W' in var_name:
        bounds = [-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
        cmap = 'jet'
    elif 'Q' in var_name:
        bounds = 10.**np.arange(-6,0,1)
        cmap = 'ocean_r'

    # Plot 
    ax[0].contourf( xcoor, ycoor, storm_mean[MP[0]][DA[0]]['cycle0']['input_tangential_wind'], cmap=cmap,levels=tan_bounds,extend='both')
    wind_contour = ax[0].contour(xcoor,ycoor,storm_mean[MP[0]][DA[0]]['cycle0']['input_radial_wind'],levels=wind_level,cmap=wind_level_cmap,linewidths=2)
    plt.clabel(wind_contour,wind_level,inline=True, fmt="%i", use_clabeltext=True, fontsize=18)
    for ida in DA:
        j = DA.index(ida)+1
        data = storm_mean[MP[0]][ida]['cycle0']
        wind_contourf = ax[j].contourf( xcoor, ycoor, data['output_tangential_wind'], cmap=cmap,levels=tan_bounds,extend='both')
        wind_contour = ax[j].contour(xcoor,ycoor,data['output_radial_wind'],levels=wind_level,cmap=wind_level_cmap,linewidths=2)
        plt.clabel(wind_contour,wind_level,inline=True, fmt="%i", use_clabeltext=True, fontsize=18)

    # Colorbar
    cbaxes = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    wind_bar = fig.colorbar(wind_contourf,cax=cbaxes,fraction=0.046, pad=0.04) #Make a colorbar for the ContourSet returned by the contourf call.
    wind_bar.ax.set_ylabel('Wind Speed (m/s)',fontsize=23)
    wind_bar.ax.tick_params(labelsize='22')

   # subplot title and labels
    if use_pressure:
        pass
    else:
        ylabel_like = [0.0,5.0,10.0,15.0,20.0]
        yticks = []
        list_y_range = list(y_range)
        for it in ylabel_like:
            yticks.append( list_y_range.index(it) )
        fig.text(0.06,0.5,'Height (km)',ha='center',va='center',rotation='vertical',fontsize=28)

    for idx in range(3):
        ax.flat[idx].set_ylim(ymin=f_yinterp(0),ymax=f_yinterp(20.5)) # cut off data above 25km 
        ax.flat[idx].set_yticks( yticks )
        ax.flat[idx].set_yticklabels( [str(it) for it in ylabel_like],fontsize=23 )

    if 'wind' in var_name:
        xb_title = 'Xb:'+var_name+' (m/s)'
        xa_title = 'Xa:'+var_name+' (m/s)'
    elif 'W' in var_name:
        xb_title = 'Xb: vertical velocity(m/s)'
        xa_title = 'Xa: vertical velocity (m/s)'
    elif 'Q' in var_name:
        xb_title = 'Xb:'+var_name+' (kg/kg)'
        xa_title = 'Xa:'+var_name+' (kg/kg)'

    ax[0].set_title( 'EnKF Prior', fontsize = 25, fontweight='bold')
    ax[1].set_title( 'conv: Posterior', fontsize = 25,  fontweight='bold')
    ax[2].set_title( 'conv+IR: Posterior', fontsize = 25, fontweight='bold')

    # set X label
    xlabel_like = np.linspace(0,300,7)
    xticks = []
    for it in xlabel_like:
        xticks.append(  f_xinterp( it )) 

    for idx in range(3):
        ax[idx].set_xticks( xticks )
        ax[idx].set_xticklabels(  [str(it) for it in xlabel_like],fontsize=18 )
        ax[idx].set_xlabel('Radius (km)',fontsize=20)
        ax[idx].set_xlim(xmin=f_xinterp(0),xmax=f_xinterp(300))

    # set X label
    xlabel_like = np.linspace(0,300,7)
    xticks = []
    for it in xlabel_like:
        xticks.append(  f_xinterp( it ) )

    for idx in range(3):
        ax[idx].set_xticks( xticks )
        ax[idx].set_xticklabels(  [str(it) for it in xlabel_like],fontsize=18 )
        ax[idx].set_xlabel('Radius (km)',fontsize=20)
        ax[idx].set_xlim(xmin=f_xinterp(0),xmax=f_xinterp(300))

    # title
    title_name = 'Azimuthally Averaged Wind (m/s)'
    fig.suptitle(title_name, fontsize=20, fontweight='bold')

    # Save figures
    figure_des=plot_dir+'AZmean'+var_name+'_cycle0.png'
    plt.savefig(figure_des, dpi=400)
    print('Saving the figure: ', figure_des)


def plot_twoTimes_azmean( storm_mean ):

    # ------------------ Plot -----------------------
    fig, ax=plt.subplots(2, 2, sharex='all', sharey='all', figsize=(15,15), dpi=400)

    # Set colorbar
    if 'radial_wind' in var_name:
        bounds = np.linspace(-12,12,9)
        cmap = 'PiYG_r'
    elif 'tangential_wind' in var_name:
        color_intervals = list(np.linspace(0,50,6))
        color_intervals.insert(0,-10.0)
        exist_cmap = plt.cm.jet
        colors = exist_cmap(np.linspace(0,1,len(color_intervals)))
        cmap = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals))
        bounds = color_intervals
    elif 'hor_wind' in var_name:
        # tangential_wind
        color_intervals = list(np.linspace(0,21,8))
        color_intervals.insert(0,-5.0)
        exist_cmap = plt.cm.rainbow
        colors = exist_cmap(np.linspace(0,1,len(color_intervals)))
        cmap = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals))
        tan_bounds = color_intervals
        # radial wind
        wind_level = [-2,-1,1,2]
        colors = [(0, 0, 1), (1, 1, 1), (1, 0, 0)]  # Blue, White, Red
        cmap_name = 'blue_black_red'
        wind_level_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=len(wind_level))

    elif 'W' in var_name:
        bounds = [-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
        cmap = 'jet'
    elif 'Q' in var_name:
        bounds = 10.**np.arange(-6,0,1)
        cmap = 'ocean_r'

    # x axis: radius 
    xv = np.linspace(0,400,100)
    x_axis_rg = range(len(xv))
    f_xinterp = interpolate.interp1d( xv, x_axis_rg)
    # y axis: model vertical coordinate
    if use_pressure:
        y_range = np.arange(900,50,50)
    else:
        y_range = np.arange(0,31,1)
    y_axis_rg = range(len(y_range))
    f_yinterp = interpolate.interp1d( y_range, y_axis_rg)

    for imp in MP:
        for i in range(2):
            # Make a mesh grid
            yv = f_yinterp( storm_mean[imp][DA[0]]['cycle0']['ver_coor'] )
            xcoor, ycoor = np.meshgrid( x_axis_rg, yv )
            # Plot 
            for ida in DA:
                j = DA.index(ida)
                data = storm_mean[imp][ida]['cycle0']
                if i == 0:
                    ax[i,j].contourf( xcoor, ycoor, data['input_tangential_wind'], cmap=cmap,levels=tan_bounds,extend='both')
                    wind_contour = ax[i,j].contour(xcoor,ycoor,data['input_radial_wind'],levels=wind_level,cmap=wind_level_cmap,linewidths=2)
                    plt.clabel(wind_contour,wind_level,inline=True, fmt="%i", use_clabeltext=True, fontsize=18)
                else:
                    wind_contourf = ax[i,j].contourf( xcoor, ycoor, data['output_tangential_wind'], cmap=cmap,levels=tan_bounds,extend='both')
                    wind_contour = ax[i,j].contour(xcoor,ycoor,data['output_radial_wind'],levels=wind_level,cmap=wind_level_cmap,linewidths=2)
                    plt.clabel(wind_contour,wind_level,inline=True, fmt="%i", use_clabeltext=True, fontsize=18)

    # Colorbar
    caxes = fig.add_axes([0.12, 0.03, 0.8, 0.02])
    wind_bar = fig.colorbar(wind_contourf,ax=ax[0:2],orientation="horizontal", cax=caxes)
    wind_bar.ax.tick_params(labelsize=18)

    # subplot title and labels
    if use_pressure:
        pass
    else:
        ylabel_like = [0.0,5.0,10.0,15.0,20.0]
        yticks = []
        list_y_range = list(y_range)
        for it in ylabel_like:
            yticks.append( list_y_range.index(it) )
        fig.text(0.06,0.5,'Height (km)',ha='center',va='center',rotation='vertical',fontsize=25)

    for idx in range(4):
        ax.flat[idx].set_ylim(ymin=f_yinterp(0),ymax=f_yinterp(20.5)) # cut off data above 25km 
        ax.flat[idx].set_yticks( yticks )
        ax.flat[idx].set_yticklabels( [str(it) for it in ylabel_like],fontsize=20 )

    if 'wind' in var_name:
        xb_title = 'Xb:'+var_name+' (m/s)'
        xa_title = 'Xa:'+var_name+' (m/s)'
    elif 'W' in var_name:
        xb_title = 'Xb: vertical velocity(m/s)'
        xa_title = 'Xa: vertical velocity (m/s)'
    elif 'Q' in var_name:
        xb_title = 'Xb:'+var_name+' (kg/kg)'
        xa_title = 'Xa:'+var_name+' (kg/kg)'

    ax[0,0].set_title( DA[0], fontsize = 20, fontweight='bold')
    ax[0,1].set_title( DA[1], fontsize = 20, fontweight='bold')

    # set X label
    xlabel_like = np.linspace(0,300,7)
    xticks = []
    for it in xlabel_like:
        xticks.append(  f_xinterp( it ) )

    for idx in range(2):
        ax[1,idx].set_xticks( xticks )
        ax[1,idx].set_xticklabels(  [str(it) for it in xlabel_like],fontsize=18 )
        ax[1,idx].set_xlabel('Radius (KM)',fontsize=20)
        ax[1,idx].set_xlim(xmin=f_xinterp(0),xmax=f_xinterp(300))

    # title
    title_name = 'Azimuthally Averaged Wind (m/s)'
    fig.suptitle(title_name, fontsize=20, fontweight='bold')

    # Save figures
    figure_des=plot_dir+'AZmean'+var_name+'_cycle0.png'
    plt.savefig(figure_des, dpi=400)
    print('Saving the figure: ', figure_des)


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['JOSE','MARIA','IRMA']#['HARVEY','IRMA','JOSE','MARIA']
    DA = ['conv','IR']
    MP = ['WSM6',]
    var_name = 'hor_wind' #hor_wind

    DAtimes = {'IRMA':['201709030000',],'JOSE':['201709050000'],'MARIA':['201709160000']}
    # model dimension
    xmax = 297
    ymax = 297
    nLevel = 42

    use_pressure = False
    azimuths = np.linspace(0,360,73) #units: degree
    radius = np.linspace(0,300,100) # units: km
    Plot_sys_azmean = True

    #-----------------------------------

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = {}
        for imp in MP:
            Exper_names[istorm][imp] = {}
            for ida in DA:
                Exper_names[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # Obtain azimuthally averaged variable for each experiment
    Exper_AZ = {}
    for istorm in Storms:
        Exper_AZ[istorm] = {}
        for imp in MP:
            Exper_AZ[istorm][imp] = {}
            for ida in DA:
                iExper = Exper_names[istorm][imp][ida]
                if iExper is not None:
                    Exper_AZ[istorm][imp][ida] = load_AZmean_oneExper(istorm,imp,ida,Exper_names,DAtimes)
                else:
                    Exper_AZ[istorm][imp][ida] = None  

    # Average over storms for each kind of experiment
    storm_mean = {}
    for imp in MP:
        storm_mean[imp] = {}
        for ida in DA:
            storm_mean[imp][ida] = {}
            tmp_mean_vcoor = np.zeros((nLevel,))
            tmp_tan_input = np.zeros((nLevel,len(radius))) 
            tmp_tan_output = np.zeros((nLevel,len(radius)))
            tmp_rad_input = np.zeros((nLevel,len(radius))) 
            tmp_rad_output = np.zeros((nLevel,len(radius)))
            for cycle in range(1):
                storm_mean[imp][ida]['cycle'+str(cycle)] = {}
                for istorm in Storms:
                    tmp_mean_vcoor = tmp_mean_vcoor + Exper_AZ[istorm][imp][ida][DAtimes[istorm][cycle]]['d_coor']['ver_coor']
                    tmp_tan_input = tmp_tan_input + Exper_AZ[istorm][imp][ida][DAtimes[istorm][cycle]]['d_var']['input']['tangential_wind']
                    tmp_tan_output = tmp_tan_output + Exper_AZ[istorm][imp][ida][DAtimes[istorm][cycle]]['d_var']['output']['tangential_wind']
                    tmp_rad_input = tmp_rad_input + Exper_AZ[istorm][imp][ida][DAtimes[istorm][cycle]]['d_var']['input']['radial_wind']
                    tmp_rad_output = tmp_rad_output + Exper_AZ[istorm][imp][ida][DAtimes[istorm][cycle]]['d_var']['output']['radial_wind']
                storm_mean[imp][ida]['cycle'+str(cycle)]['ver_coor'] = tmp_mean_vcoor/len(Storms)
                storm_mean[imp][ida]['cycle'+str(cycle)]['input_tangential_wind'] = tmp_tan_input/len(Storms)
                storm_mean[imp][ida]['cycle'+str(cycle)]['output_tangential_wind'] = tmp_tan_output/len(Storms)
                storm_mean[imp][ida]['cycle'+str(cycle)]['input_radial_wind'] = tmp_rad_input/len(Storms)
                storm_mean[imp][ida]['cycle'+str(cycle)]['output_radial_wind'] = tmp_rad_output/len(Storms)


    # Plot storm-averaged distribution
    if Plot_sys_azmean:
        # Create plot dir
        plot_dir = small_dir+'/SYSTEMS/Vis_analyze/Model/'
        plotdir_exists = os.path.exists( plot_dir )
        if plotdir_exists == False:
            os.mkdir(plot_dir)
        # plot
        plot_0thCycle_expers( storm_mean  )
        #plot_twoTimes_azmean( storm_mean )







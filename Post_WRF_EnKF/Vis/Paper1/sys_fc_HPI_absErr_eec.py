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
import matplotlib.patches as mpatches

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

# Exper_error: storms * microphysics * DA * fc_init 
def calculate_abs_error():

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
                Exper_error[istorm][imp][ida] = {}
                iExper = Exper_names[istorm][imp][ida]
                if os.path.exists( wrf_dir+'/'+istorm+'/'+iExper ):
                    print('Calculating HPI error for '+istorm+' :'+Exper_names[istorm][imp][ida])
                    if sameNum_sample:
                        tmp_err = Track.SameL_error_eachInit(istorm,wrf_dir,iExper,best_track,Exper_initTimes,DF_model_end,fc_run_hrs)
                        for key in tmp_err.keys():
                            tmp_err[key][2:,:] = np.abs( tmp_err[key][2:,:] ) # MSLP and Vamx error
                            Exper_error[istorm][imp][ida][key] = tmp_err[key]
                    else:
                        tmp_err = Track.DiffL_error_eachInit(istorm,wrf_dir,iExper,best_track,Exper_initTimes,DF_model_end)
                        for key in tmp_err.keys():
                            tmp_err[key][2:,:] = np.abs( tmp_err[key][2:,:] ) # MSLP and Vamx error
                            Exper_error[istorm][imp][ida][key] = tmp_err[key]
                else:
                    Exper_error[istorm][imp][ida] = None
    return Exper_error


# For each storm, each MP scheme, at an initialization time, there are num(DA) forecasts.
# Calculate the mean absolute error over sysnoptic times, for each forecast.
# Find the max(MAE) and min(MAE) of these forecasts, normalize MAEs of these forecasts.
def normalize_MAE():

    # Calculate MAEs
    MAE = {}
    for istorm in Storms:
        fc_inits = fc_iniT( istorm )
        MAE[istorm] = {}
        for imp in MP:
            MAE[istorm][imp] = {}
            for ida in DA:
                MAE[istorm][imp][ida] = {}
                for fc_init in fc_inits:
                    # only calculate MAE for absolute error of track, MSLP, and Vmax
                    MAE[istorm][imp][ida][fc_init] = np.nanmean(Exper_error[istorm][imp][ida][fc_init][[0,2,3],:],axis=1) 

    # normalize MAEs
    normMAE = {}
    for istorm in Storms:
        fc_inits = fc_iniT( istorm )
        normMAE[istorm] = {}
        for imp in MP:
            normMAE[istorm][imp] = {}
            # for each storm,MP: normalize MAE
            across_error = np.empty((3, 1))
            # assemble MAEs of forecasts of DA kinds * fc_inits
            for fc_init in fc_inits:
                # assemble MAEs of forecasts of DA kinds
                for ida in DA:
                    across_error = np.concatenate( (across_error,MAE[istorm][imp][ida][fc_init].reshape((3, 1))),axis=1)
            across_error = across_error[:,1:]
            min_error = np.amin(across_error,axis=1)
            max_error = np.amax(across_error,axis=1)
                
            # normalize
            for fc_init in fc_inits:
                normMAE[istorm][imp][fc_init] = {}
                for ida in DA:
                    normMAE[istorm][imp][fc_init][ida] = (MAE[istorm][imp][ida][fc_init]-min_error)/(max_error-min_error)

    return normMAE

# Calculate FSP (frequency of superior performance)
# For a forecast, fraction of times that the experiment is superior to the control
# 1>= FSP > 0.5: Experiment is usually superior
# 0<= FSP < 0.5: Experiment is usually inferior
def calculate_FSP():

    fsp = {}

    for istorm in Storms:
        fc_inits = fc_iniT( istorm )
        fsp[istorm] = {}
        for imp in MP:
            fsp[istorm][imp] = {}
            for ida in DA[1:]:
                fsp[istorm][imp][ida] = {}
                for fc_init in fc_inits:
                    fc_err = Exper_error[istorm][imp][ida][fc_init][[0,2,3],:]
                    conv_err = Exper_error[istorm][imp]['CONV'][fc_init][[0,2,3],:]
                    judge = fc_err <= conv_err
                    fsp[istorm][imp][ida][fc_init] = np.sum(judge,axis=1)/np.shape(fc_err)[1]

    return fsp


# ------------------------------------------------------------------------------------------------------
#            Operation: Plot 
# ------------------------------------------------------------------------------------------------------


# Custom handler for multiple markers in a row
class HandlerMultipleMarkers:
    def __init__(self, color, alphas):
        self.alphas = alphas
        self.color = color

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        for i, alpha in enumerate(self.alphas):
            handlebox.add_artist(plt.Line2D([x0 + (i + 0.5) * width / len(self.alphas)], [y0 + height / 2], linestyle='',
                                            marker='o', markersize=6, markerfacecolor=self.color, alpha=alpha))
        return handlebox


# layout:
# WSM6_normMAE, WSM6_FSP
# THO_normMAE, THO_FSP
def plot_2by2():

    # Set up figure
    fig = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=2,top=0.93,right=0.95,hspace=0.09)

    ax = {}
    for imp in MP:
        ax[imp] = {}
        ir = MP.index(imp)
        ist_grids = outer_grids[ir].subgridspec(1, 2, wspace=0.0)
        ax[imp]['mae'] = fig.add_subplot( ist_grids[0] )
        ax[imp]['fsp'] = fig.add_subplot( ist_grids[1] )

    # Customization
    colors = {}
    colorset = {'HARVEY':'#FF13EC','JOSE':'#0D13F6','MARIA':'#FFA500','IRMA':'#DB2824'}
    alphas = np.linspace(0.2,1,4)
    #marker_type= {'CONV':'o','IR':'P','IR+MW':'X'}
    marker_type= {'CONV':'P','IR':'o','IR+MW':'s'}


    # Plot simulations
    for imp in MP:
        for ist in Storms:
            fc_inits = fc_iniT( ist )
            for ida in DA:
                for fc_init in fc_inits:
                    # Plot MAE
                    ax[imp]['mae'].set_aspect('equal', 'box')
                    x = np.linspace(-0.02,1.02)
                    ax[imp]['mae'].plot(x, x, color='grey', linestyle='-', linewidth=1)
                    mae_track = normMAE[ist][imp][fc_init][ida][0]
                    mae_mslp = normMAE[ist][imp][fc_init][ida][1] # Vmax: [2]
                    ax[imp]['mae'].scatter(mae_mslp,mae_track,s=30,marker=marker_type[ida],facecolor=colorset[ist],alpha=alphas[fc_inits.index(fc_init)])
                    # Plot FSP
                    ax[imp]['fsp'].set_aspect('equal', 'box')
                    if ida == 'CONV':
                        continue
                    # add mpatches to distinguish inferior and superior area
                    if_corners = [(-0.02, -0.02), (0.5, -0.02), (0.5, 0.5),(-0.02, 0.5)]
                    sp_corners = [(0.5, 0.5), (1.02, 0.5), (1.02, 1.02),(0.5, 1.02)]
                    # Add a patch (polygon) to the figure based on the corners
                    if_polygon = mpatches.Polygon(if_corners, closed=True, edgecolor='r', facecolor='none')
                    sp_polygon = mpatches.Polygon(sp_corners, closed=True, edgecolor='b', facecolor='none')
                    ax[imp]['fsp'].add_patch(if_polygon)
                    ax[imp]['fsp'].add_patch(sp_polygon)

                    fsp_track = fsp[ist][imp][ida][fc_init][0]
                    fsp_mslp = fsp[ist][imp][ida][fc_init][1] # Vmax: [2]
                    ax[imp]['fsp'].scatter(fsp_mslp,fsp_track,s=30,marker=marker_type[ida],facecolor=colorset[ist],alpha=alphas[fc_inits.index(fc_init)],label='_nolegend_')

    # Add legends
    # DAs
    scatter_DA = ax[MP[0]]['mae'].collections
    legend_DA = ax[MP[0]]['mae'].legend([scatter_DA[i] for i in [0,4,8]],DA,fontsize='7',loc='upper left')
    # Add the first legend manually to the current Axes
    ax[MP[0]]['mae'].add_artist(legend_DA)
    # Storms
    # function to create rows for the legend
    def create_legend_rows(colors, alphas, labels):
        rows = [{'color': color, 'alphas': alphas, 'label': label} for color, label in zip(colors, labels)]
        return rows
    # define colors, alphas, and labels for the legend rows
    legend_colors = ['#FF13EC','#0D13F6','#FFA500','#DB2824']
    legend_alphas = alphas
    legend_labels = Storms
    # create legend rows
    rows = create_legend_rows(legend_colors, legend_alphas, legend_labels)

    handles = []
    for row in rows:
        handle = mpatches.FancyBboxPatch((0, 0), 1, 1, color='none')  # Dummy handle
        handles.append(handle) 

    ax[MP[1]]['fsp'].legend(handles, [row['label'] for row in rows],
          handler_map={handle: HandlerMultipleMarkers(row['color'], row['alphas']) for handle, row in zip(handles, rows)},
          loc='upper right')

    #scatter_storm = ax[MP[1]]['fsp'].collections
    #legend_storm = ax[MP[1]]['fsp'].legend([scatter_storm[i] for i in [0,8,16,24]],Storms,fontsize='7',loc='upper right')
    #ax[MP[1]]['fsp'].add_artist(legend_storm) 

    # Add texts
    for imp in MP:
        ax[imp]['fsp'].text(0.38,0.05,'inferior', horizontalalignment='center', fontsize=10, color='red')
        ax[imp]['fsp'].text(0.88,0.55,'superior', horizontalalignment='center', fontsize=10, color='blue')

    fig.text(0.05,0.73,'WSM6', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.05,0.32,'THO', fontsize=10, ha='center', va='center',rotation='vertical')

    # Set attributes
    for imp in MP:
        # limitation
        ax[imp]['mae'].set_xlim([-0.02,1.02])
        ax[imp]['mae'].set_ylim([-0.02,1.02])
        ax[imp]['fsp'].set_xlim([-0.02,1.02])
        ax[imp]['fsp'].set_ylim([-0.02,1.02])
        # hide yticklabel
        ax[imp]['fsp'].set_yticklabels([])
        
        # X/Y labels
        ax[imp]['mae'].set_ylabel('Track: Normalized MAE')
        ax[imp]['fsp'].set_ylabel('Track: FSP')
        if MP.index( imp ) == 1:
            ax[imp]['mae'].set_xlabel('MSLP: Normalized MAE')
            ax[imp]['fsp'].set_xlabel('MSLP: FSP')
    

    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_absError.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )






if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY','JOSE','IRMA','MARIA']
    DA = ['CONV','IR','IR+MW']
    MP = ['WSM6','THO'] #

    # if operate over the same number of samples for all forecasts
    sameNum_sample = False
    if sameNum_sample:
        fc_run_hrs = 60
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
    Exper_error = calculate_abs_error()
    
    # Calculate normalized MAE
    normMAE = normalize_MAE()

    # Calculate FSP
    fsp = calculate_FSP()

    plot_2by2()

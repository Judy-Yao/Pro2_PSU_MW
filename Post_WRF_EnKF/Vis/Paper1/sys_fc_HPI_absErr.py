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
import matplotlib.lines as mlines

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

def storm_color():
    colorset = {'HARVEY':'#FFA500','IRMA':'#FF13EC','JOSE':'#0D13F6','MARIA':'#097969',}
    return colorset

def alpha_fc_init():
    alphas = np.linspace(0.2,1,4)
    return alphas

def DA_marker():
    marker_type= {'CONV':'P','IR':'o','IR+MW':'s'}
    return marker_type


# Exper_error: storms * microphysics * DA * fc_init 
def calculate_abs_error():

    # Calculate the error with lead times of various length
    Exper_error = {}
    for istorm in Storms:
        Exper_error[istorm] = {}
        # Read best-track data
        Btk_start, Btk_end = t_range_btk( istorm )
        best_track = UD.btk_in_duration(small_dir,istorm, Btk_start, Btk_end, hour_step=6)
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


# Calculate the mean absolute error over sysnoptic times, for each forecast.
def calc_MAE_overTime():

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
    return MAE


# For each storm, each MP scheme, at an initialization time, there are num(DA) forecasts.
# Calculate the mean absolute error over sysnoptic times, for each forecast.
# Find the max(MAE) and min(MAE) of these forecasts, normalize MAEs of these forecasts.
def normalize_MAE():

    # Calculate MAEs
    MAE = calc_MAE_overTime()

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
def FSP_overCONV():

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

def FSP_overIR():

    fsp = {}

    for istorm in Storms:
        fc_inits = fc_iniT( istorm )
        fsp[istorm] = {}
        for imp in MP:
            fsp[istorm][imp] = {}
            fsp[istorm][imp]['IR+MW'] = {}
            for fc_init in fc_inits:
                fc_err = Exper_error[istorm][imp]['IR+MW'][fc_init][[0,2,3],:]
                ir_err = Exper_error[istorm][imp]['IR'][fc_init][[0,2,3],:]
                judge = fc_err <= ir_err
                fsp[istorm][imp]['IR+MW'][fc_init] = np.sum(judge,axis=1)/np.shape(fc_err)[1]

    return fsp


# For each storm, mean of absolute errors is averaged over all forecasts 
# with respect to the shortest forecast period
def calculate_MAEs_wrt_time():

    # find the shortest forecast period 
    fc_srt_len = {}
    for ist in Storms:
        fc_inits = fc_iniT( ist )
        fc_lens = [ np.shape( Exper_error[ist][MP[0]][DA[0]][fc_init] )[1] for fc_init in fc_inits] 
        fc_srt_len[ist] = np.min( fc_lens )

    # calculate the mean absolute error
    MAEs = {}
    for istorm in Storms:
        fc_inits = fc_iniT( istorm )
        MAEs[istorm] = {}
        for imp in MP:
            MAEs[istorm][imp] = {}
            for ida in DA:
                MAEs[istorm][imp][ida] = {}
                # assemble errors of forecasted attributes of len(fc_inits)
                across_track = np.empty((1,fc_srt_len[istorm]))
                across_MSLP = np.empty((1,fc_srt_len[istorm]))
                across_Vmax = np.empty((1,fc_srt_len[istorm]))
                for fc_init in fc_inits:
                    # only assemble errors to the end of the shortest forecast period
                    across_track = np.concatenate( (across_track,Exper_error[istorm][imp][ida][fc_init][0,:fc_srt_len[istorm]].reshape(1,fc_srt_len[istorm])),axis=0) 
                    across_MSLP = np.concatenate( (across_MSLP,Exper_error[istorm][imp][ida][fc_init][2,:fc_srt_len[istorm]].reshape(1,fc_srt_len[istorm])),axis=0) 
                    across_Vmax = np.concatenate( (across_MSLP,Exper_error[istorm][imp][ida][fc_init][3,:fc_srt_len[istorm]].reshape(1,fc_srt_len[istorm])),axis=0) 
                across_track = across_track[1:,:] # remove the first empty row
                across_MSLP = across_MSLP[1:,:]
                across_Vmax = across_Vmax[1:,:]
                # average over all forecasts for that forecast kind
                MAEs[istorm][imp][ida]['track'] = np.nanmean(across_track,axis=0)
                MAEs[istorm][imp][ida]['mslp'] = np.nanmean(across_MSLP,axis=0)
                MAEs[istorm][imp][ida]['vmax'] = np.nanmean(across_Vmax,axis=0)
    # special treatment to abnormal experiments
    for ida in ['IR','IR+MW']:
        MAEs['JOSE']['THO'][ida]['track_all'] = MAEs['JOSE']['THO'][ida]['track'] 
        MAEs['JOSE']['THO'][ida]['mslp_all'] = MAEs['JOSE']['THO'][ida]['mslp'] 
        MAEs['JOSE']['THO'][ida]['vmax_all'] = MAEs['JOSE']['THO'][ida]['vmax'] 

    # abnormal points: JOSE,THO,20151200
    for ida in ['IR','IR+MW']:
        across_track = np.empty((1,fc_srt_len[istorm]))
        across_MSLP = np.empty((1,fc_srt_len[istorm]))
        across_Vmax = np.empty((1,fc_srt_len[istorm]))
        fc_inits = fc_iniT( 'JOSE' )
        for fc_init in fc_inits:
            if fc_init == '201709051200':
                continue
            if ida == 'IR+MW' and fc_init == '201709050600':
                continue
            across_track = np.concatenate( (across_track,Exper_error['JOSE']['THO'][ida][fc_init][0,:fc_srt_len[istorm]].reshape(1,fc_srt_len['JOSE'])),axis=0)
            across_MSLP = np.concatenate( (across_MSLP,Exper_error['JOSE']['THO'][ida][fc_init][2,:fc_srt_len[istorm]].reshape(1,fc_srt_len['JOSE'])),axis=0)
            across_Vmax = np.concatenate( (across_MSLP,Exper_error['JOSE']['THO'][ida][fc_init][3,:fc_srt_len[istorm]].reshape(1,fc_srt_len['JOSE'])),axis=0)
        across_track = across_track[1:,:] # remove the first empty row
        across_MSLP = across_MSLP[1:,:]
        across_Vmax = across_Vmax[1:,:]
        # average over all forecasts for that forecast kind
        MAEs['JOSE']['THO'][ida]['track'] = np.nanmean(across_track,axis=0)
        MAEs['JOSE']['THO'][ida]['mslp'] = np.nanmean(across_MSLP,axis=0)
        MAEs['JOSE']['THO'][ida]['vmax'] = np.nanmean(across_Vmax,axis=0)

    return fc_srt_len,MAEs

def MAE_overStorms():

    mae_sts = {}
    for imp in MP:
        mae_sts[imp] = {}
        for ida in DA:
            mae_sts[imp][ida] = {}
            # initialize containers
            across_track = np.full( (len(Storms),max(fc_srt_len.values())), np.nan)
            across_mslp = np.full( (len(Storms),max(fc_srt_len.values())), np.nan)
            across_vmax = np.full( (len(Storms),max(fc_srt_len.values())), np.nan)
            # loop thru storms
            for ist in Storms:
                across_track[Storms.index(ist),:fc_srt_len[ist]] = MAEs[ist][imp][ida]['track']
                across_mslp[Storms.index(ist),:fc_srt_len[ist]] = MAEs[ist][imp][ida]['mslp'] 
                across_vmax[Storms.index(ist),:fc_srt_len[ist]] = MAEs[ist][imp][ida]['vmax'] 
            # average
            mae_sts[imp][ida]['track'] = np.nanmean(across_track,axis=0)
            mae_sts[imp][ida]['mslp'] = np.nanmean(across_mslp,axis=0)
            mae_sts[imp][ida]['vmax'] = np.nanmean(across_vmax,axis=0)

    return mae_sts


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
    colorset = storm_color()
    alphas = alpha_fc_init()
    #marker_type= {'CONV':'o','IR':'P','IR+MW':'X'}
    marker_type= DA_marker()

    # Plot simulations
    for imp in MP:
        ax[imp]['mae'].set_aspect('equal', 'box')
        x = np.linspace(-0.02,1.02)
        ax[imp]['mae'].plot(x, x, color='grey', linestyle='-', linewidth=1)

        ax[imp]['fsp'].set_aspect('equal', 'box')
        # add mpatches to distinguish inferior and superior area
        if_corners = [(-0.02, -0.02), (0.5, -0.02), (0.5, 0.5),(-0.02, 0.5)]
        sp_corners = [(0.5, 0.5), (1.02, 0.5), (1.02, 1.02),(0.5, 1.02)]
        # Add a patch (polygon) to the figure based on the corners
        if_polygon = mpatches.Polygon(if_corners, closed=True, edgecolor='r', facecolor='none')
        sp_polygon = mpatches.Polygon(sp_corners, closed=True, edgecolor='b', facecolor='none')
        ax[imp]['fsp'].add_patch(if_polygon)
        ax[imp]['fsp'].add_patch(sp_polygon)

        for ist in Storms:
            fc_inits = fc_iniT( ist )
            for ida in DA:
                for fc_init in fc_inits:
                    # Plot MAE
                    mae_track = normMAE[ist][imp][fc_init][ida][0]
                    mae_mslp = normMAE[ist][imp][fc_init][ida][1] # Vmax: [2]
                    ax[imp]['mae'].scatter(mae_mslp,mae_track,s=30,marker=marker_type[ida],facecolor=colorset[ist],alpha=alphas[fc_inits.index(fc_init)])
                    # Plot FSP
                    ax[imp]['fsp'].set_aspect('equal', 'box')
                    if ida == 'CONV':
                        continue
                    fsp_track = fsp[ist][imp][ida][fc_init][0]
                    fsp_mslp = fsp[ist][imp][ida][fc_init][1] # Vmax: [2]
                    ax[imp]['fsp'].scatter(fsp_mslp,fsp_track,s=30,marker=marker_type[ida],facecolor=colorset[ist],alpha=alphas[fc_inits.index(fc_init)],label='_nolegend_')

    # Add legends
    # DAs
    scatter_DA = ax[MP[0]]['mae'].collections
    legend_DA = ax[MP[0]]['mae'].legend([scatter_DA[i] for i in [3,7,11]],DA,fontsize='7',loc='upper left')
    # Add the first legend manually to the current Axes
    ax[MP[0]]['mae'].add_artist(legend_DA)
    # Storms
    # function to create rows for the legend
    def create_legend_rows(colors, alphas, labels):
        rows = [{'color': color, 'alphas': alphas, 'label': label} for color, label in zip(colors, labels)]
        return rows
    # define colors, alphas, and labels for the legend rows
    legend_colors = list(colorset.values()) #['#FFA500','#0D13F6','#FF13EC','#097969'] ##097969
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


# layout:
# WSM6_normMAE, WSM6_FSP
# THO_normMAE, THO_FSP
# Different DA is distinguished thru colors
def plot_2by2_byColor():

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
    colorset = storm_color()
    alphas = alpha_fc_init()
    colorset= {'CONV':'black','IR':'red','IR+MW':'blue'}   #DA_marker()

    # Plot simulations
    for imp in MP:
        ax[imp]['mae'].set_aspect('equal', 'box')
        x = np.linspace(-0.02,1.02)
        ax[imp]['mae'].plot(x, x, color='grey', linestyle='-', linewidth=1)

        ax[imp]['fsp'].set_aspect('equal', 'box')
        # add mpatches to distinguish inferior and superior area
        if_corners = [(-0.02, -0.02), (0.5, -0.02), (0.5, 0.5),(-0.02, 0.5)]
        sp_corners = [(0.5, 0.5), (1.02, 0.5), (1.02, 1.02),(0.5, 1.02)]
        # Add a patch (polygon) to the figure based on the corners
        if_polygon = mpatches.Polygon(if_corners, closed=True, edgecolor='gray', facecolor='none')
        sp_polygon = mpatches.Polygon(sp_corners, closed=True, edgecolor='black', facecolor='none')
        ax[imp]['fsp'].add_patch(if_polygon)
        ax[imp]['fsp'].add_patch(sp_polygon)

        for ist in Storms:
            fc_inits = fc_iniT( ist )
            for ida in DA:
                for fc_init in fc_inits:
                    # Plot MAE
                    mae_track = normMAE[ist][imp][fc_init][ida][0]
                    mae_mslp = normMAE[ist][imp][fc_init][ida][1] # Vmax: [2]
                    ax[imp]['mae'].scatter(mae_mslp,mae_track,s=30,facecolor=colorset[ida],alpha=0.6)
                    # Plot FSP
                    ax[imp]['fsp'].set_aspect('equal', 'box')
                    if ida == 'CONV':
                        continue
                    fsp_track = fsp[ist][imp][ida][fc_init][0]
                    fsp_mslp = fsp[ist][imp][ida][fc_init][1] # Vmax: [2]
                    ax[imp]['fsp'].scatter(fsp_mslp,fsp_track,s=30,facecolor=colorset[ida],alpha=0.6)

    # Add legends
    # DAs
    scatter_DA = ax[MP[0]]['mae'].collections
    legend_DA = ax[MP[0]]['mae'].legend([scatter_DA[i] for i in [3,7,11]],DA,fontsize='7',loc='upper left')
    # Add the first legend manually to the current Axes
    ax[MP[0]]['mae'].add_artist(legend_DA)

    # Add texts
    for imp in MP:
        ax[imp]['fsp'].text(0.38,0.05,'inferior', horizontalalignment='center', fontsize=10, color='gray')
        ax[imp]['fsp'].text(0.88,0.55,'superior', horizontalalignment='center', fontsize=10, color='black')

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
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_absError_byColor.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

    return None


# layout:
# WSM6_normMAE, WSM6_FSP
# THO_normMAE, THO_FSP
# Above: different DA is distinguished thru colors
# WSM6_FSP_overIR, THO_FSP_overIR
# Different colors+transparencies to indicate FSP_overIR
def plot_3by2():

    matplotlib.rcParams['font.size'] = 10

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=3,top=0.95,right=0.95,hspace=0.12)

    row = ['mae','fsp','fsp_ir']

    ax = {}
    for i in range(3):
        ax[row[i]] = {}
        for imp in MP:
            ir = MP.index(imp)
            ist_grids = outer_grids[i].subgridspec(1, 2, wspace=0.05)
            ax[row[i]][imp] = fig.add_subplot( ist_grids[ir] )

    # Customization
    storm_colors = storm_color()
    alphas = alpha_fc_init()
    DA_color= {'CONV':'black','IR':'red','IR+MW':'blue'} 

    # Plot simulations
    for imp in MP:
        #ax['mae'][imp].set_aspect('equal', 'box')
        x = np.linspace(-0.02,1.02)
        ax['mae'][imp].plot(x, x, color='grey', linestyle='-', linewidth=1)

        #ax['fsp'][imp].set_aspect('equal', 'box')
        # add mpatches to distinguish inferior and superior area
        if_corners = [(-0.02, -0.02), (0.5, -0.02), (0.5, 0.5),(-0.02, 0.5)]
        sp_corners = [(0.5, 0.5), (1.02, 0.5), (1.02, 1.02),(0.5, 1.02)]
        # Add a patch (polygon) to the figure based on the corners
        if_polygon = mpatches.Polygon(if_corners, closed=True, edgecolor='gray', facecolor='none')
        sp_polygon = mpatches.Polygon(sp_corners, closed=True, edgecolor='black', facecolor='none')
        ax['fsp'][imp].add_patch(if_polygon)
        ax['fsp'][imp].add_patch(sp_polygon)

		# Plot MAE and FSP relative to CONV
        for ist in Storms:
            fc_inits = fc_iniT( ist )
            for ida in DA:
                for fc_init in fc_inits:
                    # Plot MAE
                    mae_track = normMAE[ist][imp][fc_init][ida][0]
                    mae_mslp = normMAE[ist][imp][fc_init][ida][1] # Vmax: [2]
                    ax['mae'][imp].scatter(mae_mslp,mae_track,s=30,facecolor=DA_color[ida],alpha=0.6)
                    # Plot FSP relative to CONV
                    if ida == 'CONV':
                        continue
                    fsp_track = fsp[ist][imp][ida][fc_init][0]
                    fsp_mslp = fsp[ist][imp][ida][fc_init][1] # Vmax: [2]
                    ax['fsp'][imp].scatter(fsp_mslp,fsp_track,s=30,facecolor=DA_color[ida],alpha=0.6)

		# Plot FSP relative to IR
        for ist in Storms:
            fc_inits = fc_iniT( ist )
            #ax['fsp_ir'][imp].set_aspect('equal', 'box')
            # add mpatches to distinguish inferior and superior area
            if_corners = [(-0.02, -0.02), (0.5, -0.02), (0.5, 0.5),(-0.02, 0.5)]
            sp_corners = [(0.5, 0.5), (1.02, 0.5), (1.02, 1.02),(0.5, 1.02)]
            # Add a patch (polygon) to the figure based on the corners
            if_polygon = mpatches.Polygon(if_corners, closed=True, edgecolor='gray', facecolor='none')
            sp_polygon = mpatches.Polygon(sp_corners, closed=True, edgecolor='black', facecolor='none')
            ax['fsp_ir'][imp].add_patch(if_polygon)
            ax['fsp_ir'][imp].add_patch(sp_polygon)

            for fc_init in fc_inits:
                fsp_track = fsp_ir[ist][imp]['IR+MW'][fc_init][0]
                fsp_mslp = fsp_ir[ist][imp]['IR+MW'][fc_init][1] # Vmax: [2]
                ax['fsp_ir'][imp].scatter(fsp_mslp,fsp_track,s=30,marker='o',facecolor=storm_colors[ist],alpha=alphas[fc_inits.index(fc_init)],label='_nolegend_')			
	
    # Add legends
    # DAs
    scatter_DA = ax['mae'][MP[0]].collections
    legend_DA = ax['mae'][MP[0]].legend([scatter_DA[i] for i in [3,7,11]],DA,fontsize='10',loc='upper left')
    # Add the first legend manually to the current Axes
    ax['mae'][MP[0]].add_artist(legend_DA)

    # Storms
    # function to create rows for the legend
    def create_legend_rows(colors, alphas, labels):
        rows = [{'color': color, 'alphas': alphas, 'label': label} for color, label in zip(colors, labels)]
        return rows
        
    # define colors, alphas, and labels for the legend rows
    legend_colors = list(storm_colors.values()) #['#FFA500','#0D13F6','#FF13EC','#097969'] ##097969
    legend_alphas = alphas
    legend_labels = Storms
    # create legend rows
    rows = create_legend_rows(legend_colors, legend_alphas, legend_labels)

    handles = []
    for row in rows:
        handle = mpatches.FancyBboxPatch((0, 0), 1, 1, color='none')  # Dummy handle
        handles.append(handle)

    ax['fsp_ir'][MP[1]].legend(handles, [row['label'] for row in rows],
          handler_map={handle: HandlerMultipleMarkers(row['color'], row['alphas']) for handle, row in zip(handles, rows)},
          loc='lower right',fontsize=8)

    # Add texts
    for imp in MP:
        ax['fsp'][imp].text(0.38,0.05,'inferior', horizontalalignment='center', fontsize=10, color='gray')
        ax['fsp_ir'][imp].text(0.38,0.05,'inferior', horizontalalignment='center', fontsize=10, color='gray')
        ax['fsp'][imp].text(0.64,0.55,'superior', horizontalalignment='center', fontsize=10, color='black')
        ax['fsp_ir'][imp].text(0.64,0.55,'superior', horizontalalignment='center', fontsize=10, color='black')
    
    # Add title
    for imp in MP:
        ax['mae'][imp].set_title( imp,fontsize=12 )

    # Set attributes
    for imp in MP:
        # limitation
        ax['mae'][imp].set_xlim([-0.02,1.02])
        ax['mae'][imp].set_ylim([-0.02,1.02])
        ax['fsp'][imp].set_xlim([-0.02,1.02])
        ax['fsp'][imp].set_ylim([-0.02,1.02])
        ax['fsp_ir'][imp].set_xlim([-0.02,1.02])
        ax['fsp_ir'][imp].set_ylim([-0.02,1.02])
        # hide yticklabel
        ax['fsp'][imp].set_xticklabels([])
        ax['mae'][imp].set_xticklabels([])
        if MP.index( imp ) == 1:
            ax['fsp_ir'][imp].set_yticklabels([])
            ax['fsp'][imp].set_yticklabels([])
            ax['mae'][imp].set_yticklabels([])
        # X/Y labels
        if MP.index( imp ) == 0:
            ax['mae'][imp].set_ylabel('Track: Normalized MAE',fontsize=10)
            ax['fsp'][imp].set_ylabel('Track: FSP (rel. CONV)',fontsize=10)
            ax['fsp_ir'][imp].set_ylabel('Track: FSP (rel. IR)',fontsize=10)
        ax['mae'][imp].set_xlabel('MSLP: Normalized MAE',fontsize=10)
        ax['fsp'][imp].set_xlabel('MSLP: FSP (rel. CONV)',fontsize=10)
        ax['fsp_ir'][imp].set_xlabel('MSLP: FSP (rel. IR)',fontsize=10)


    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_absError_3by2.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

# layout:
# WSM6_FSP, WSM6_FSP
# Above: different DA is distinguished thru colors
# WSM6_FSP_overIR, THO_FSP_overIR
# Different colors+transparencies to indicate FSP_overIR
def plot_2by2_FSP():

    matplotlib.rcParams['font.size'] = 10

    # Set up figure
    fig = plt.figure( figsize=(6.5,6.5),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=2,top=0.95,right=0.95,hspace=0.12)

    row = ['fsp','fsp_ir']

    ax = {}
    for i in range(2):
        ax[row[i]] = {}
        for imp in MP:
            ir = MP.index(imp)
            ist_grids = outer_grids[i].subgridspec(1, 2, wspace=0.05)
            ax[row[i]][imp] = fig.add_subplot( ist_grids[ir] )

    # Customization
    storm_colors = storm_color()
    alphas = alpha_fc_init()
    DA_color= {'IR':'red','IR+MW':'blue'}

    # Plot simulations
    for imp in MP:
        # add mpatches to distinguish inferior and superior area
        if_corners = [(-0.02, -0.02), (0.5, -0.02), (0.5, 0.5),(-0.02, 0.5)]
        sp_corners = [(0.5, 0.5), (1.02, 0.5), (1.02, 1.02),(0.5, 1.02)]
        # Add a patch (polygon) to the figure based on the corners
        if_polygon = mpatches.Polygon(if_corners, closed=True, edgecolor='gray', facecolor='none')
        sp_polygon = mpatches.Polygon(sp_corners, closed=True, edgecolor='black', facecolor='none')
        ax['fsp'][imp].add_patch(if_polygon)
        ax['fsp'][imp].add_patch(sp_polygon)

        # Plot FSP relative to CONV
        for ist in Storms:
            fc_inits = fc_iniT( ist )
            for ida in DA:
                for fc_init in fc_inits:
                    # Plot FSP relative to CONV
                    if ida == 'CONV':
                        continue
                    fsp_track = fsp[ist][imp][ida][fc_init][0]
                    fsp_mslp = fsp[ist][imp][ida][fc_init][1] # Vmax: [2]
                    ax['fsp'][imp].scatter(fsp_mslp,fsp_track,s=30,facecolor=DA_color[ida],alpha=0.6)

        # Plot FSP relative to IR
        for ist in Storms:
            fc_inits = fc_iniT( ist )
            #ax['fsp_ir'][imp].set_aspect('equal', 'box')
            # add mpatches to distinguish inferior and superior area
            if_corners = [(-0.02, -0.02), (0.5, -0.02), (0.5, 0.5),(-0.02, 0.5)]
            sp_corners = [(0.5, 0.5), (1.02, 0.5), (1.02, 1.02),(0.5, 1.02)]
            # Add a patch (polygon) to the figure based on the corners
            if_polygon = mpatches.Polygon(if_corners, closed=True, edgecolor='gray', facecolor='none')
            sp_polygon = mpatches.Polygon(sp_corners, closed=True, edgecolor='black', facecolor='none')
            ax['fsp_ir'][imp].add_patch(if_polygon)
            ax['fsp_ir'][imp].add_patch(sp_polygon)

            for fc_init in fc_inits:
                fsp_track = fsp_ir[ist][imp]['IR+MW'][fc_init][0]
                fsp_mslp = fsp_ir[ist][imp]['IR+MW'][fc_init][1] # Vmax: [2]
                ax['fsp_ir'][imp].scatter(fsp_mslp,fsp_track,s=30,marker='o',facecolor=storm_colors[ist],alpha=alphas[fc_inits.index(fc_init)],label='_nolegend_')

    # Add legends
    # DAs
    scatter_DA = ax['fsp'][MP[1]].collections
    legend_DA = ax['fsp'][MP[1]].legend([scatter_DA[i] for i in [3,7]],['IR','IR+MW'],fontsize='10',loc='upper right')
    # Add the first legend manually to the current Axes
    ax['fsp'][MP[1]].add_artist(legend_DA)

    # Storms
    # function to create rows for the legend
    def create_legend_rows(colors, alphas, labels):
        rows = [{'color': color, 'alphas': alphas, 'label': label} for color, label in zip(colors, labels)]
        return rows

    # define colors, alphas, and labels for the legend rows
    legend_colors = list(storm_colors.values()) #['#FFA500','#0D13F6','#FF13EC','#097969'] ##097969
    legend_alphas = alphas
    legend_labels = Storms
    # create legend rows
    rows = create_legend_rows(legend_colors, legend_alphas, legend_labels)

    handles = []
    for row in rows:
        handle = mpatches.FancyBboxPatch((0, 0), 1, 1, color='none')  # Dummy handle
        handles.append(handle)

    ax['fsp_ir'][MP[1]].legend(handles, [row['label'] for row in rows],
          handler_map={handle: HandlerMultipleMarkers(row['color'], row['alphas']) for handle, row in zip(handles, rows)},
          loc='lower right',fontsize=8)

    # Add texts
    for imp in MP:
        ax['fsp'][imp].text(0.38,0.05,'inferior', horizontalalignment='center', fontsize=10, color='gray')
        ax['fsp_ir'][imp].text(0.38,0.05,'inferior', horizontalalignment='center', fontsize=10, color='gray')
        ax['fsp'][imp].text(0.64,0.55,'superior', horizontalalignment='center', fontsize=10, color='black')
        ax['fsp_ir'][imp].text(0.64,0.55,'superior', horizontalalignment='center', fontsize=10, color='black')

    # Add title
    for imp in MP:
        ax['fsp'][imp].set_title( imp,fontsize=12 )

    # Set attributes
    for imp in MP:
        # limitation
        ax['fsp'][imp].set_xlim([-0.02,1.02])
        ax['fsp'][imp].set_ylim([-0.02,1.02])
        ax['fsp_ir'][imp].set_xlim([-0.02,1.02])
        ax['fsp_ir'][imp].set_ylim([-0.02,1.02])
        # hide yticklabel
        ax['fsp'][imp].set_xticklabels([])
        if MP.index( imp ) == 1:
            ax['fsp_ir'][imp].set_yticklabels([])
            ax['fsp'][imp].set_yticklabels([])
        # X/Y labels
        if MP.index( imp ) == 0:
            ax['fsp'][imp].set_ylabel('Track: FSP (rel. CONV)',fontsize=10)
            ax['fsp_ir'][imp].set_ylabel('Track: FSP (rel. IR)',fontsize=10)
        ax['fsp'][imp].set_xlabel('MSLP: FSP (rel. CONV)',fontsize=10)
        ax['fsp_ir'][imp].set_xlabel('MSLP: FSP (rel. IR)',fontsize=10)


    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_absError_2by2_FSP.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )



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

# layout:
# WSM6_MAE_track, WSM6_MAE_MSLP
# THO_MAE_track, THO_MAE_MSLP
def plot_2by2_MAEs():

    # Set up figure
    fig = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=2,top=0.93,right=0.96,hspace=0.03)

    ax = {}
    for imp in MP:
        ax[imp] = {}
        ir = MP.index(imp)
        ist_grids = outer_grids[ir].subgridspec(1, 2, wspace=0.15)
        ax[imp]['track'] = fig.add_subplot( ist_grids[0] )
        ax[imp]['mslp'] = fig.add_subplot( ist_grids[1] )
 
    # Customization
    #colorset = storm_color()
    #marker_type= DA_marker()
   
    # Option2 : DA by color; Storm by marker
    colorset = {'CONV': 'black','IR':'red','IR+MW':'blue'}
    marker_type = {'HARVEY':'|','IRMA':'x','JOSE':'_','MARIA':'.'}

    line_width = {'HARVEY':1,'IRMA':1,'JOSE':1,'MARIA':1,'Mean':2.5}

    # Plot individual simulations
    #subax_track = add_subplot_axes(ax['THO']['track'],[0.1,0.66,0.45,0.3]) # small subplot
    #subax_mslp = add_subplot_axes(ax['THO']['mslp'],[0.1,0.66,0.45,0.3]) # small subplot

    for imp in MP:
        for ist in Storms:
            lead_t = list(range(0,fc_srt_len[ist]))
            for ida in DA:
                ax[imp]['track'].plot(lead_t,MAEs[ist][imp][ida]['track'],color=colorset[ida],marker=marker_type[ist],markersize=5,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4) 
                ax[imp]['mslp'].plot(lead_t,MAEs[ist][imp][ida]['mslp'],color=colorset[ida],marker=marker_type[ist],markersize=5,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4)
                if ist == 'JOSE' and imp == "THO":
                    if ida == 'CONV':
                        pass
                    #else:
                        #subax_track.plot(lead_t,MAEs[ist][imp][ida]['track_all'],color=colorset[ida],marker=marker_type[ist],markersize=3,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4)
                        #subax_track.set_xlim( [-0.1,max(fc_srt_len.values())-1] )
                        #subax_track.set_yticks( [0,450,900,1350] )
                        #subax_track.set_ylim([-0.1,1350.1])
                        #subax_track.set_yticklabels(['0','450','900','1350'],fontsize=6)
                        #subax_track.set_xticks([0,4,8,12,16] )
                        #subax_track.set_xticklabels(['D0','D1','D2','D3','D4'])
                        #subax_track.grid(True,linewidth=0.8, color='gray', alpha=0.3, linestyle='-')

                        #subax_mslp.plot(lead_t,MAEs[ist][imp][ida]['mslp_all'],color=colorset[ida],marker=marker_type[ist],markersize=3,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4)
                        #subax_mslp.set_xlim( [-0.1,max(fc_srt_len.values())-1] )
                        #subax_mslp.set_yticks( [0,45,90] )
                        #subax_mslp.set_ylim([-0.1,90])
                        #subax_mslp.set_yticklabels(['0','45','90'],fontsize=6)
                        #subax_mslp.set_xticks([0,4,8,12,16] )
                        #subax_mslp.set_xticklabels(['D0','D1','D2','D3','D4'])
                        #subax_mslp.grid(True,linewidth=0.8, color='gray', alpha=0.3, linestyle='-')

    # Plot mean over storms
    mae_sts = MAE_overStorms()
    lead_t = list(range(0,max(fc_srt_len.values())))
    alphas = {'CONV':0.5,'IR':0.7,'IR+MW':0.9}
    for imp in MP:
        for ida in DA:
            ax[imp]['track'].plot(lead_t,mae_sts[imp][ida]['track'],color=colorset[ida],linewidth=3)
            ax[imp]['mslp'].plot(lead_t,mae_sts[imp][ida]['mslp'],color=colorset[ida],linewidth=3)


    # Manully add legends
    # create proxy artists for the legend with different line widths and colors

    ## Option 1: DA by marker; Storm by color
    #lgd_1 = Storms + ['Mean']
    #legend_colors = list(colorset.values())
    #legend_colors.append( 'black' )
    #list_widths = list(line_width.values())
    #proxy_artists = [plt.Line2D([0], [0], color=color, lw=lw) for color,lw in zip( legend_colors,list_widths )]
    #legend0.set_alpha( 0.5 )
    # Add the first legend manually to the current Axes
    #ax['WSM6']['track'].legend(proxy_artists,lgd_1,fontsize='8',loc='upper center',ncol=2)

    ## Option 2: DA by color; Storm by marker
    lines = ax['WSM6']['track'].get_lines()
    lgd_0 = Storms + ['Mean']
    legend = ax['WSM6']['track'].legend([lines[i] for i in [0,3,6,9,-3]], lgd_0,fontsize='8',loc='upper center',ncol=2)
    # Add the first legend manually to the current Axes
    ax['WSM6']['track'].add_artist(legend)

    lines = ax['WSM6']['mslp'].get_lines()
    lgd_0 = DA
    legend0 = ax['WSM6']['mslp'].legend([lines[i] for i in [-3,-2,-1]], lgd_0,fontsize='8',loc='upper left')
    #legend0.set_alpha( 0.5 )
    # Add the first legend manually to the current Axes
    ax['WSM6']['mslp'].add_artist(legend0)

    fig.text(0.32,0.95,'Track', fontsize=12, ha='center', va='center',rotation='horizontal')
    fig.text(0.78,0.95,'MSLP', fontsize=12, ha='center', va='center',rotation='horizontal')
    fig.text(0.03,0.73,'WSM6', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.03,0.32,'THO', fontsize=10, ha='center', va='center',rotation='vertical')

    # Set axis attributes
    ax['WSM6']['track'].set_ylim( [-0.1,400] )
    ax['THO']['track'].set_ylim( [-0.1,800] )
    ax['WSM6']['mslp'].set_ylim( [-0.1,45] )
    ax['THO']['mslp'].set_ylim( [-0.1,90] )
    # y ticks
    wsm6_yticks = list(range(0,45+1,5))
    ax['WSM6']['mslp'].set_yticks( wsm6_yticks )
    ax['WSM6']['mslp'].set_yticklabels( [str(it) for it in wsm6_yticks] )
    tho_yticks = list(range(0,90+1,15))
    ax['THO']['mslp'].set_yticks( tho_yticks)
    ax['THO']['mslp'].set_yticklabels( [str(it) for it in tho_yticks] )
    # x axis
    for imp in MP:
        ax[imp]['track'].set_xlim( [-0.1,max(fc_srt_len.values())-1] )
        ax[imp]['mslp'].set_xlim( [-0.1,max(fc_srt_len.values())-1] )
        ax[imp]['track'].set_xticks([0,4,8,12,16] )
        ax[imp]['mslp'].set_xticks([0,4,8,12,16] )
        if MP.index(imp) != 0:
            ax[imp]['track'].set_xticklabels(['D0','D1','D2','D3','D4'])
            ax[imp]['mslp'].set_xticklabels(['D0','D1','D2','D3','D4'])
            ax[imp]['track'].set_xlabel('Forecast Time (days)')
            ax[imp]['mslp'].set_xlabel('Forecast Time (days)')
        else:
            ax[imp]['track'].set_xticklabels(['','','','',''])
            ax[imp]['mslp'].set_xticklabels(['','','','',''])
        # y label
        ax[imp]['track'].set_ylabel('Track: MAE (km)')
        ax[imp]['mslp'].set_ylabel('MSLP: MAE (hPa)')
        # grid lines
        ax[imp]['track'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')
        ax[imp]['mslp'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_absError_MAEs_withFCtime.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

    return None



# layout
# one MSLP MAE
def plot_MSLP_MAE():

    # Set up figure
    fig,ax = plt.subplots( 1,1,figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    
    # DA by color
    colorset = {'CONV': 'black','IR':'red','IR+MW':'blue'}

    # 
    for imp in MP:
        for ist in Storms:
            lead_t = list(range(0,fc_srt_len[ist]))
            for ida in DA:
                if imp == 'WSM6':
                    ax.plot(lead_t,MAEs[ist][imp][ida]['mslp'],color=colorset[ida],linewidth=3)

    # labels
    ax.set_title('Mean Absolute Error of MSLP (hPa)',fontsize=15)
    #fig.text(0.78,0.95,'MSLP', fontsize=12, ha='center', va='center',rotation='horizontal')

    lines = ax.get_lines()
    legend = ax.legend([lines[i] for i in [0,1,2]], list(colorset.keys()),fontsize='15',loc='upper center',ncol=3)
    # Add the first legend manually to the current Axes
    ax.add_artist(legend)

    # Set axis attributes
    ax.set_ylim( [-0.1,25] )
    # y ticks
    wsm6_yticks = list(range(0,25+1,5))
    ax.set_yticks( wsm6_yticks )
    ax.set_yticklabels( [str(it) for it in wsm6_yticks],fontsize='15' )
    # x axis
    for imp in MP:
        ax.set_xlim( [-0.1,max(fc_srt_len.values())-1] )
        ax.set_xticks([0,4,8,12,16] )
        ax.set_xticklabels(['D0','D1','D2','D3','D4'],fontsize='15')
        ax.set_xlabel('Forecast Time (days)',fontsize='15')
        # y label
        ax.set_ylabel('MSLP: MAE (hPa)',fontsize='15')
        # grid lines
        ax.grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    # Save figure
    des_name = 'WSM6_MSLP_JOSE.png'#small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_absError_MAEs_withFCtime.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


# layout:
# WSM6_FSP_overIR, THO_FSP_overIR
def plot_1by2_FSP():

    # Set up figure
    fig = plt.figure( figsize=(6.5,3.25),dpi=200) # standard: 6.5,8.5
    grids = fig.add_gridspec(ncols=2,nrows=1,top=0.93,hspace=0.0,wspace=0.05)
    ax = {}
    for imp in MP:
        ax[imp] = fig.add_subplot( grids[ MP.index(imp) ] )

    # Customization
    colorset = storm_color()
    alphas = alpha_fc_init()
    marker_type= DA_marker()

    # Plot simulations
    for imp in MP:
        for ist in Storms:
            fc_inits = fc_iniT( ist )
            ax[imp].set_aspect('equal', 'box')
            # add mpatches to distinguish inferior and superior area
            if_corners = [(-0.02, -0.02), (0.5, -0.02), (0.5, 0.5),(-0.02, 0.5)]
            sp_corners = [(0.5, 0.5), (1.02, 0.5), (1.02, 1.02),(0.5, 1.02)]
            # Add a patch (polygon) to the figure based on the corners
            if_polygon = mpatches.Polygon(if_corners, closed=True, edgecolor='r', facecolor='none')
            sp_polygon = mpatches.Polygon(sp_corners, closed=True, edgecolor='b', facecolor='none')
            ax[imp].add_patch(if_polygon)
            ax[imp].add_patch(sp_polygon)

            for fc_init in fc_inits:
                fsp_track = fsp[ist][imp]['IR+MW'][fc_init][0]
                fsp_mslp = fsp[ist][imp]['IR+MW'][fc_init][1] # Vmax: [2]
                ax[imp].scatter(fsp_mslp,fsp_track,s=30,marker='o',facecolor=colorset[ist],alpha=alphas[fc_inits.index(fc_init)],label='_nolegend_')

    # Add legends
    # Storms
    # function to create rows for the legend
    def create_legend_rows(colors, alphas, labels):
        rows = [{'color': color, 'alphas': alphas, 'label': label} for color, label in zip(colors, labels)]
        return rows
    # define colors, alphas, and labels for the legend rows
    legend_colors = list(colorset.values()) #['#FFA500','#0D13F6','#FF13EC','#097969'] ##097969
    legend_alphas = alphas
    legend_labels = Storms
    # create legend rows
    rows = create_legend_rows(legend_colors, legend_alphas, legend_labels)

    handles = []
    for row in rows:
        handle = mpatches.FancyBboxPatch((0, 0), 1, 1, color='none')  # Dummy handle
        handles.append(handle)

    ax[MP[1]].legend(handles, [row['label'] for row in rows],
          handler_map={handle: HandlerMultipleMarkers(row['color'], row['alphas']) for handle, row in zip(handles, rows)},
          loc='lower right')

    # Add texts
    for imp in MP:
        ax[imp].text(0.38,0.05,'inferior', horizontalalignment='center', fontsize=10, color='red')
        ax[imp].text(0.88,0.55,'superior', horizontalalignment='center', fontsize=10, color='blue')

    # Add title
    for imp in MP:
        ax[imp].set_title( imp,fontsize=12 )

    # Set x/y attributes
    for imp in MP:
        # limitation
        ax[imp].set_xlim([-0.02,1.02])
        ax[imp].set_ylim([-0.02,1.02])
        # hide yticklabel
        ax[imp].set_yticklabels([])
        # X/Y labels
        if MP.index( imp ) == 0:
            ax[imp].set_ylabel('Track: FSP')
        ax[imp].set_xlabel('MSLP: FSP')

    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_absError_FSPoverIR.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

def test_plot():

    matplotlib.rcParams['font.size'] = 8

    # Calculate MAE
    MAE = calc_MAE_overTime()

    # Set up figure
    fig,ax = plt.subplots( figsize=(6.5,3.25),dpi=200) # standard: 6.5,8.5
    
    # Customization
    colorset = {'CONV':'black','IR':'red','IR+MW':'blue'}
    alphas = alpha_fc_init()

    # Plot
    for ist in Storms:
        fc_inits = fc_iniT( ist )
        for imp in MP:
            for ida in DA:
                for fc_init in fc_inits:
                    mae_track = MAE[ist][imp][ida][fc_init][0]
                    mae_mslp = MAE[ist][imp][ida][fc_init][2]
                    ax.scatter(mae_mslp,mae_track,s=30,marker=marker_type[ida],facecolor=colorset[imp],alpha=alphas[fc_inits.index(fc_init)],label='_nolegend_')

    # Add legends
    # DAs
    scatter_DA = ax.collections
    legend_DA = ax.legend([scatter_DA[i] for i in [3,7,11]],DA,fontsize='7',loc='upper left')
    # Add the first legend manually to the current Axes
    ax.add_artist(legend_DA)
    # Storms
    # function to create rows for the legend
    def create_legend_rows(colors, alphas, labels):
        rows = [{'color': color, 'alphas': alphas, 'label': label} for color, label in zip(colors, labels)]
        return rows
    # define colors, alphas, and labels for the legend rows
    legend_colors = list(colorset.values()) #['#FFA500','#0D13F6','#FF13EC','#097969'] ##097969
    legend_alphas = alphas
    legend_labels = MP
    # create legend rows
    rows = create_legend_rows(legend_colors, legend_alphas, legend_labels)

    handles = []
    for row in rows:
        handle = mpatches.FancyBboxPatch((0, 0), 1, 1, color='none')  # Dummy handle
        handles.append(handle)

    ax.legend(handles, [row['label'] for row in rows],
          handler_map={handle: HandlerMultipleMarkers(row['color'], row['alphas']) for handle, row in zip(handles, rows)},
          loc='lower right')

    ax.set_ylim([0,300]) 
    ax.set_xlim([0,30])
    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_MAE_track_mslp.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )




if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/Clean_results/'

    #--------Configuration------------
    Storms = ['JOSE']  #['HARVEY','IRMA','JOSE','MARIA']
    DA = ['CONV','IR','IR+MW']
    MP = ['WSM6','THO'] #['WSM6','THO'] #

    # if operate over the same number of samples for all forecasts
    sameNum_sample = False
    if sameNum_sample:
        fc_run_hrs = 60

    MAE_only = False
    normMAE_FSP = False
    MAE_wrt_time = True
    FSP_IR = False
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

    # Plot

    if MAE_only:
        test_plot()


    if normMAE_FSP and not FSP_IR:
        # Calculate normalized MAE
        normMAE = normalize_MAE()

        # Calculate FSP
        fsp = FSP_overCONV()

        plot_2by2_byColor()
        #plot_2by2()

    if MAE_wrt_time:
        # Calculate MAE with same samples
        fc_srt_len, MAEs = calculate_MAEs_wrt_time()

        plot_MSLP_MAE()
        #plot_2by2_MAEs()

    if FSP_IR and not normMAE_FSP:
        fsp = FSP_overIR()
        plot_1by2_FSP()

    if normMAE_FSP and FSP_IR:
        # Calculate normalized MAE
        normMAE = normalize_MAE()

        # Calculate FSP relative to CONV
        fsp = FSP_overCONV()

        # Calculate FSP relative to IR
        fsp_ir = FSP_overIR()
    
        #plot_3by2()
        plot_2by2_FSP()

    #plot_4by2():unable to plot it with 6.5 inch wide and 8.5 inch tall
    # layout:
    # WSM6_normMAE, WSM6_FSP
    # THO_normMAE, THO_FSP
    # WSM6_MAE_track, WSM6_MAE_MSLP
    # THO_MAE_track, THO_MAE_MSLP











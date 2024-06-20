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

# For each storm, mean of absolute errors is averaged over all forecasts 
# with respect to the shortest forecast period
def calculate_MAEs():

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
    colorset = storm_color()
    #marker_type= {'CONV':'o','IR':'P','IR+MW':'X'}
    marker_type= DA_marker()
    line_width = {'HARVEY':1,'IRMA':1,'JOSE':1,'MARIA':1,'Mean':1.5}

    # Plot individual simulations
    subax = add_subplot_axes(ax['THO']['track'],[0.1,0.66,0.45,0.3]) # small subplot

    for imp in MP:
        for ist in Storms:
            lead_t = list(range(0,fc_srt_len[ist]))
            for ida in DA:
                ax[imp]['track'].plot(lead_t,MAEs[ist][imp][ida]['track'],color=colorset[ist],marker=marker_type[ida],markersize=5,linewidth=line_width[ist],alpha=0.4) 
                ax[imp]['mslp'].plot(lead_t,MAEs[ist][imp][ida]['mslp'],color=colorset[ist],marker=marker_type[ida],markersize=5,linewidth=line_width[ist],alpha=0.4)
                if ist == 'JOSE' and imp == "THO":
                    if ida == 'CONV':
                        pass
                    else:
                        subax.plot(lead_t,MAEs[ist][imp][ida]['track'],color=colorset[ist],marker=marker_type[ida],markersize=3,linewidth=line_width[ist],alpha=0.4)
                        subax.set_xlim( [-0.1,max(fc_srt_len.values())-1] )
                        subax.set_yticks( [0,450,900,1350] )
                        subax.set_ylim([-0.1,1350.1])
                        subax.set_yticklabels(['0','450','900','1350'],fontsize=6)
                        subax.set_xticks([0,4,8,12,16] )
                        subax.set_xticklabels(['D0','D1','D2','D3','D4'])


    # Plot mean over storms
    mae_sts = MAE_overStorms()
    lead_t = list(range(0,max(fc_srt_len.values())))
    alphas = {'CONV':0.5,'IR':0.7,'IR+MW':0.9}
    for imp in MP:
        for ida in DA:
            ax[imp]['track'].plot(lead_t,mae_sts[imp][ida]['track'],color='black',marker=marker_type[ida],markersize=6,linewidth=1.5,alpha=alphas[ida])
            ax[imp]['mslp'].plot(lead_t,mae_sts[imp][ida]['mslp'],color='black',marker=marker_type[ida],markersize=6,linewidth=1.5,alpha=alphas[ida])


    # Manully add legends
    # create proxy artists for the legend with different line widths and colors
    lgd_1 = Storms + ['Mean']
    legend_colors = list(colorset.values())
    legend_colors.append( 'black' )
    list_widths = list(line_width.values())
    proxy_artists = [plt.Line2D([0], [0], color=color, lw=lw) for color,lw in zip( legend_colors,list_widths )]
    #legend0.set_alpha( 0.5 )
    # Add the first legend manually to the current Axes
    ax['WSM6']['track'].legend(proxy_artists,lgd_1,fontsize='8',loc='upper center',ncol=2)

    lines = ax['WSM6']['mslp'].get_lines()
    lgd_0 = DA
    legend0 = ax['WSM6']['mslp'].legend([lines[i] for i in [-3,-2,-1]], lgd_0,fontsize='8',loc='upper left')
    #legend0.set_alpha( 0.5 )
    # Add the first legend manually to the current Axes
    ax['WSM6']['mslp'].add_artist(legend0)

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

# layout:
# WSM6_normMAE, WSM6_FSP
# THO_normMAE, THO_FSP
# WSM6_MAE_track, WSM6_MAE_MSLP
# THO_MAE_track, THO_MAE_MSLP
def plot_4by2():

    # Set up figure
    fig = plt.figure( figsize=(6.5,8.5),dpi=200) # standard: 6.5,8.5
    outer_grids = fig.add_gridspec(ncols=1,nrows=2,top=0.93,right=0.95,hspace=0.09)

    ax = {}
    for imp in MP:
        ax[imp] = {}

    # gridspec inside gridspec
    kind1_grid = outer_grids[0].subgridspec(2, 2, wspace=0.0)
    for imp in MP:
        ir = MP.index(imp)
        ax[imp]['mae'] = fig.add_subplot( kind1_grid[ir,0])
        ax[imp]['fsp'] = fig.add_subplot( kind1_grid[ir,1])

    kind2_grid = outer_grids[1].subgridspec(2, 2, wspace=0.15)
    for imp in MP: 
        ir = MP.index(imp)
        ax[imp]['track'] = fig.add_subplot( kind2_grid[ir,0])
        ax[imp]['mslp'] = fig.add_subplot( kind2_grid[ir,1])

    # Customization
    colorset = storm_color()
    alphas = alpha_fc_init()
    #marker_type= {'CONV':'o','IR':'P','IR+MW':'X'}
    marker_type= DA_marker()
    line_width = {'HARVEY':1,'IRMA':1,'JOSE':1,'MARIA':1,'Mean':1.5}

    # Plot individual simulations
    subax = add_subplot_axes(ax['THO']['track'],[0.1,0.66,0.45,0.3]) # small subplot

    for imp in MP:
        for ist in Storms:
            fc_inits = fc_iniT( ist )
            for ida in DA:
                # Plot metrics of each forecast
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

                # Plot mean absolute error with respect to forecast time
                lead_t = list(range(0,fc_srt_len[ist]))
                ax[imp]['track'].plot(lead_t,MAEs[ist][imp][ida]['track'],color=colorset[ist],marker=marker_type[ida],markersize=5,linewidth=line_width[ist],alpha=0.4)
                ax[imp]['mslp'].plot(lead_t,MAEs[ist][imp][ida]['mslp'],color=colorset[ist],marker=marker_type[ida],markersize=5,linewidth=line_width[ist],alpha=0.4)
                if ist == 'JOSE' and imp == "THO":
                    if ida == 'CONV':
                        pass
                    else:
                        subax.plot(lead_t,MAEs[ist][imp][ida]['track'],color=colorset[ist],marker=marker_type[ida],markersize=3,linewidth=line_width[ist],alpha=0.4)
                        subax.set_yticks( [0,450,900,1350] )
                        subax.set_ylim([-0.1,1350.1])
                        subax.set_yticklabels(['0','450','900','1350'],fontsize=6)
                        subax.set_xticks([0,4,8,12,16] )
                        subax.set_xticklabels(['D0','D1','D2','D3','D4'])

    # Plot mean over storms
    mae_sts = MAE_overStorms()
    lead_t = list(range(0,max(fc_srt_len.values())))
    alpha_bt = {'CONV':0.5,'IR':0.7,'IR+MW':0.9}
    for imp in MP:
        for ida in DA:
            ax[imp]['track'].plot(lead_t,mae_sts[imp][ida]['track'],color='black',marker=marker_type[ida],markersize=6,linewidth=1.5,alpha=alpha_bt[ida])
            ax[imp]['mslp'].plot(lead_t,mae_sts[imp][ida]['mslp'],color='black',marker=marker_type[ida],markersize=6,linewidth=1.5,alpha=alpha_bt[ida])

    # Manually add legends for the top figures
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
    legend_colors = list(colorset.values()) ##097969
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

    # Manully add legends for the bottom figures
    # create proxy artists for the legend with different line widths and colors
    lgd_1 = Storms + ['Mean']
    legend_colors = list(colorset.values())
    legend_colors.append( 'black' )
    list_widths = list(line_width.values())
    proxy_artists = [plt.Line2D([0], [0], color=color, lw=lw) for color,lw in zip( legend_colors,list_widths )]
    #legend0.set_alpha( 0.5 )
    # Add the first legend manually to the current Axes
    ax['WSM6']['track'].legend(proxy_artists,lgd_1,fontsize='8',loc='upper center',ncol=2)

    lines = ax['WSM6']['mslp'].get_lines()
    lgd_0 = DA
    legend0 = ax['WSM6']['mslp'].legend([lines[i] for i in [-3,-2,-1]], lgd_0,fontsize='8',loc='upper left')
    #legend0.set_alpha( 0.5 )
    # Add the first legend manually to the current Axes
    ax['WSM6']['mslp'].add_artist(legend0)

    # Set attributes for X/Y axes
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
        if MP.index(imp) != 0:
            ax[imp]['track'].set_xticks([0,4,8,12,16] )
            ax[imp]['track'].set_xticklabels(['D0','D1','D2','D3','D4'])
            ax[imp]['mslp'].set_xticks([0,4,8,12,16] )
            ax[imp]['mslp'].set_xticklabels(['D0','D1','D2','D3','D4'])
            ax[imp]['track'].set_xlabel('Forecast Time (days)')
            ax[imp]['mslp'].set_xlabel('Forecast Time (days)')
        # y label
        ax[imp]['track'].set_ylabel('Track: MAE (km)')
        ax[imp]['mslp'].set_ylabel('MSLP: MAE (hPa)')
        # grid lines
        ax[imp]['track'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')
        ax[imp]['mslp'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    # Add y label
    fig.text(0.05,0.85,'WSM6', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.05,0.63,'THO', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.05,0.42,'WSM6', fontsize=10, ha='center', va='center' ,rotation='vertical')
    fig.text(0.05,0.19,'THO', fontsize=10, ha='center', va='center',rotation='vertical')

    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_absError.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    #--------Configuration------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
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

    # Calculate MAE with same samples
    fc_srt_len, MAEs = calculate_MAEs()

    plot_2by2_MAEs()

    #plot_4by2()












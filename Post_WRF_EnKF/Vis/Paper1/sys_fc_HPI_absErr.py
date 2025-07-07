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
from scipy.stats import wilcoxon, ttest_ind, anderson, ttest_rel

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
    marker_type= {'CONV':'P','IR':'o','MW':'s'}
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
                #print(istorm,imp,ida)
                if os.path.exists( wrf_dir+'/'+istorm+'/'+iExper ):
                    print('Calculating HPI error for '+istorm+' :'+Exper_names[istorm][imp][ida])
                    if sameNum_sample:
                        tmp_err = Track.SameL_error_eachInit(istorm,wrf_dir,iExper,best_track,Exper_initTimes,DF_model_end,fc_run_hrs)
                        for key in tmp_err.keys():
                            tmp_err[key][2:,:] = np.abs( tmp_err[key][2:,:] ) # MSLP and Vmax error
                            Exper_error[istorm][imp][ida][key] = tmp_err[key]
                    else:
                        tmp_err = Track.DiffL_error_eachInit(istorm,wrf_dir,iExper,best_track,Exper_initTimes,DF_model_end)
                        for key in tmp_err.keys(): # key: fc_init_t
                            #print('Initialization time:',key)
                            # calculate the absolute value for MSLP and Vamx error
                            tmp_err[key][2:,:] = np.abs( tmp_err[key][2:,:] )
                            # for each dictionary element: track error (km), azimuth error, mslp error, and Vmax error
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


# Modify the FPS: if two errors are less than a threshold, then this pair of errors is not significantly different from each other 
def threshold_for_FSP( threshold, error_control, error ):

    # Pairwise comparison
    superior = 0
    superior_control = 0
    not_significant = 0

    for ec, er in zip(error_control, error):
        #print('ec,er:', ec, er)
        if abs(ec - er) < threshold:
            not_significant += 1
        elif ec < er:
            superior_control += 1
        else:
            superior += 1
    print( 'not_significant',not_significant )

    # Total valid comparisons
    total_comparisons = superior_control + superior

    # Frequency of Superior Performance
    fsp = superior / total_comparisons if total_comparisons > 0 else 0

    return fsp


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
            for ida in ['IR','MW']:
                fsp[istorm][imp][ida] = {}
                for fc_init in fc_inits:
                    #print(istorm,imp,ida,fc_init)
                    fsp[istorm][imp][ida][fc_init] = []
                    # Track error
                    track_thd = 1 
                    fc_err = Exper_error[istorm][imp][ida][fc_init][0,:] 
                    conv_err = Exper_error[istorm][imp]['CONV'][fc_init][0,:]
                    #print('track: fc_err',fc_err)
                    #print('track: conv_err',conv_err)
                    fsp[istorm][imp][ida][fc_init].append( threshold_for_FSP( track_thd, conv_err, fc_err) )
                    # MSLP error
                    mslp_thd = 0.5
                    fc_err = Exper_error[istorm][imp][ida][fc_init][2,:] 
                    conv_err = Exper_error[istorm][imp]['CONV'][fc_init][2,:]
                    #print('mslp: fc_err',fc_err)
                    #print('mslp: conv_err',conv_err)
                    fsp[istorm][imp][ida][fc_init].append( threshold_for_FSP( mslp_thd, conv_err, fc_err) )
                    #judge = fc_err <= conv_err
                    #fsp[istorm][imp][ida][fc_init] = np.sum(judge,axis=1)/np.shape(fc_err)[1]

    return fsp

def FSP_overIR():

    fsp = {}

    for istorm in Storms:
        fc_inits = fc_iniT( istorm )
        fsp[istorm] = {}
        for imp in MP:
            fsp[istorm][imp] = {}
            fsp[istorm][imp]['MW'] = {}
            for fc_init in fc_inits:
                #print('FSP_overIR:',istorm,imp,ida,fc_init)
                fsp[istorm][imp][ida][fc_init] = []
                # Track error
                track_thd = 1
                fc_err = Exper_error[istorm][imp][ida][fc_init][0,:]
                ir_err = Exper_error[istorm][imp]['IR'][fc_init][0,:]
                fsp[istorm][imp]['MW'][fc_init].append( threshold_for_FSP( track_thd, ir_err, fc_err) )
                # MSLP error
                mslp_thd = 0.5
                fc_err = Exper_error[istorm][imp][ida][fc_init][2,:]
                ir_err = Exper_error[istorm][imp]['IR'][fc_init][2,:]
                fsp[istorm][imp]['MW'][fc_init].append( threshold_for_FSP( mslp_thd, ir_err, fc_err) )
                #judge = fc_err <= ir_err
                #fsp[istorm][imp]['MW'][fc_init] = np.sum(judge,axis=1)/np.shape(fc_err)[1]

    return fsp


# For each storm:
# 1. find the shortest forecast period among the four deterministic forecasts
# 2. Stack the forecasted track/mslp/vmax for one set of experiment (imp,ida)
def Assemble_err_wrt_time( Exper_error ):
    # find the shortest forecast period 
    fc_srt_len = {}
    for ist in Storms:
        fc_inits = fc_iniT( ist )
        fc_lens = [ np.shape( Exper_error[ist][MP[0]][DA[0]][fc_init] )[1] for fc_init in fc_inits] 
        fc_srt_len[ist] = np.min( fc_lens )

    # Assemble the errors into an array
    aerr = {} # assembled errors
    for ist in Storms:
        fc_inits = fc_iniT( ist )
        aerr[ist] = {}
        for imp in MP:
            aerr[ist][imp] = {}
            for ida in DA:
                aerr[ist][imp][ida] = {}
                # assemble errors of forecasted attributes of len(fc_inits)
                across_track = np.empty((1,fc_srt_len[ist]))
                across_MSLP = np.empty((1,fc_srt_len[ist]))
                across_Vmax = np.empty((1,fc_srt_len[ist]))
                for fc_init in fc_inits:
                    # only assemble errors to the end of the shortest forecast period
                    across_track = np.concatenate( (across_track,Exper_error[ist][imp][ida][fc_init][0,:fc_srt_len[ist]].reshape(1,fc_srt_len[ist])),axis=0)
                    across_MSLP = np.concatenate( (across_MSLP,Exper_error[ist][imp][ida][fc_init][2,:fc_srt_len[ist]].reshape(1,fc_srt_len[ist])),axis=0)
                    across_Vmax = np.concatenate( (across_MSLP,Exper_error[ist][imp][ida][fc_init][3,:fc_srt_len[ist]].reshape(1,fc_srt_len[ist])),axis=0)
                across_track = across_track[1:,:] # remove the first empty row
                across_MSLP = across_MSLP[1:,:]
                across_Vmax = across_Vmax[1:,:]
                # average over all forecasts for that forecast kind
                aerr[ist][imp][ida]['track'] = across_track
                aerr[ist][imp][ida]['mslp'] = across_MSLP
                aerr[ist][imp][ida]['vmax'] = across_Vmax
                    
    # special treatment to abnormal experiments: JOSE,THO_IR,20151200
    #for ida in ['IR','MW']:
    if 'JOSE' in Storms:
        for ida in ['IR',]:
            aerr['JOSE']['THO'][ida]['track_all'] = aerr['JOSE']['THO'][ida]['track']
            aerr['JOSE']['THO'][ida]['mslp_all'] = aerr['JOSE']['THO'][ida]['mslp']
            aerr['JOSE']['THO'][ida]['vmax_all'] = aerr['JOSE']['THO'][ida]['vmax']

            across_track = np.empty((1,fc_srt_len['JOSE']))
            across_MSLP = np.empty((1,fc_srt_len['JOSE']))
            across_Vmax = np.empty((1,fc_srt_len['JOSE']))
            fc_inits = fc_iniT( 'JOSE' )
            for fc_init in fc_inits:
                if fc_init == '201709051200':
                    continue
                #if ida == 'MW' and fc_init == '201709050600':
                #    continue
                across_track = np.concatenate( (across_track,Exper_error['JOSE']['THO'][ida][fc_init][0,:fc_srt_len[istorm]].reshape(1,fc_srt_len['JOSE'])),axis=0)
                across_MSLP = np.concatenate( (across_MSLP,Exper_error['JOSE']['THO'][ida][fc_init][2,:fc_srt_len[istorm]].reshape(1,fc_srt_len['JOSE'])),axis=0)
                across_Vmax = np.concatenate( (across_MSLP,Exper_error['JOSE']['THO'][ida][fc_init][3,:fc_srt_len[istorm]].reshape(1,fc_srt_len['JOSE'])),axis=0)
            across_track = across_track[1:,:] # remove the first empty row
            across_MSLP = across_MSLP[1:,:]
            across_Vmax = across_Vmax[1:,:]
            # average over all forecasts for that forecast kind
            aerr['JOSE']['THO'][ida]['track'] = across_track
            aerr['JOSE']['THO'][ida]['mslp'] = across_MSLP
            aerr['JOSE']['THO'][ida]['vmax'] = across_Vmax

    return aerr


# For each storm, mean of absolute errors is averaged over all forecasts 
# with respect to the shortest forecast period
def calculate_MAEs_wrt_time( aerr ):

    # find the shortest forecast period 
    fc_srt_len = {}
    for ist in Storms:
        fc_inits = fc_iniT( ist )
        fc_lens = [ np.shape( Exper_error[ist][MP[0]][DA[0]][fc_init] )[1] for fc_init in fc_inits]
        fc_srt_len[ist] = np.min( fc_lens )

    # calculate the mean absolute error
    MAEs = {}
    for istorm in Storms:
        MAEs[istorm] = {}
        for imp in MP:
            MAEs[istorm][imp] = {}
            for ida in DA:
                MAEs[istorm][imp][ida] = {}
                across_track = aerr[istorm][imp][ida]['track'] # 2D (num_fc,fc_len)
                across_MSLP = aerr[istorm][imp][ida]['mslp']
                across_Vmax = aerr[istorm][imp][ida]['vmax']
                # average over all forecasts for that forecast kind
                MAEs[istorm][imp][ida]['track'] = np.nanmean(across_track,axis=0)
                MAEs[istorm][imp][ida]['mslp'] = np.nanmean(across_MSLP,axis=0)
                MAEs[istorm][imp][ida]['vmax'] = np.nanmean(across_Vmax,axis=0)

    # special treatment to abnormal experiments: JOSE,THO_IR,20151200
    #for ida in ['IR','MW']:
    if 'JOSE' in Storms:
        for ida in ['IR',]:
            across_track = aerr['JOSE'][imp][ida]['track_all'] # 2D (num_fc,fc_len)
            across_MSLP = aerr['JOSE'][imp][ida]['mslp_all']
            across_Vmax = aerr['JOSE'][imp][ida]['vmax_all']
            # average over all forecasts for that forecast kind
            MAEs['JOSE'][imp][ida]['track_all'] = np.nanmean(across_track,axis=0)
            MAEs['JOSE'][imp][ida]['mslp_all'] = np.nanmean(across_MSLP,axis=0)
            MAEs['JOSE'][imp][ida]['vmax_all'] = np.nanmean(across_Vmax,axis=0)

    return fc_srt_len,MAEs


# Min-Max scale formula
def cal_min_max_scale( arr,min_v,max_v ):

    if (max_v - min_v) == 0:
        return (arr - min_v)/((max_v - min_v)+1e-10 )
    else:
        return (arr - min_v)/(max_v - min_v)

def cal_z_score( arr,miu,sigma ):

    return (arr - miu)/sigma

# For each storm, scaled mean of absolute errors is averaged over all forecasts 
# with respect to the shortest forecast period
def calculate_scaled_MAEs_wrt_time( aerr ):

    # find the shortest forecast period 
    fc_srt_len = {}
    for ist in Storms:
        fc_inits = fc_iniT( ist )
        fc_lens = [ np.shape( Exper_error[ist][MP[0]][DA[0]][fc_init] )[1] for fc_init in fc_inits]
        fc_srt_len[ist] = np.min( fc_lens )

    # Find statistics for scaling methods
    if min_max_scale:
        var = ['track','mslp','vmax']
        min_st = {}
        max_st = {}
        for ist in Storms:
            min_st[ist] = {}
            max_st[ist] = {}
            for iv in var:
                min_v = 0
                max_v = 0
                for imp in MP:
                    for ida in DA:
                        min_v = min( min_v, np.amin(aerr[ist][imp][ida][iv]) )
                        max_v = max( max_v, np.amax(aerr[ist][imp][ida][iv]) )
                min_st[ist][iv] = min_v
                max_st[ist][iv] = max_v

        #print(min_st, max_st)

    if z_score:
        var = ['track','mslp','vmax']
        miu_st = {}
        sigma_st = {}
        # Stack experiments for each storm
        for ist in Storms:
            miu_st[ist] = {}
            sigma_st[ist] = {}
            for iv in var: 
                tmp = np.empty((1,fc_srt_len[ist]))
                # stack MP and DA experiments
                for imp in MP:
                    for ida in DA:
                        tmp = np.vstack( (tmp,aerr[ist][imp][ida][iv]) )
                # remove the first empty row 
                tmp = tmp[1:,:]
                # calculate the mean and sigma
                miu_st[ist][iv] = np.mean( tmp )
                sigma_st[ist][iv] = np.std( tmp )

        #print( miu_st,sigma_st )


    # Scale the data
    # calculate the mean absolute error
    var = ['track','mslp','vmax']
    MAEs = {}
    for istorm in Storms:
        MAEs[istorm] = {}
        for imp in MP:
            MAEs[istorm][imp] = {}
            for ida in DA:
                MAEs[istorm][imp][ida] = {}
                # Loop thru vars
                for iv in var:
                    tmp = aerr[istorm][imp][ida][iv] # 2D (num_fc,fc_len)
                    # scale the data
                    if min_max_scale:
                        tmp = cal_min_max_scale( tmp,min_st[ist][iv],max_st[ist][iv] )
                        #print( tmp )
                    if z_score:
                        tmp = cal_z_score( tmp,miu_st[ist][iv],sigma_st[ist][iv] )
                    # Average over all forecasts for that forecast kind
                    MAEs[istorm][imp][ida][iv] = np.nanmean(tmp,axis=0)


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
#            Operation: Plot scaled MAEs
# ------------------------------------------------------------------------------------------------------

# layout:
# WSM6_MAE_track, WSM6_MAE_MSLP
# THO_MAE_track, THO_MAE_MSLP
def plot_2by2_scaled_MAEs():

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
    colorset = {'CONV': 'black','IR':'red','MW':'blue'}
    marker_type = {'HARVEY':'|','IRMA':'x','JOSE':'_','MARIA':'.'}
    line_width = {'HARVEY':1,'IRMA':1,'JOSE':1,'MARIA':1,'Mean':2.5}

    #for imp in MP:
    #    for ist in Storms:
    #        lead_t = list(range(0,fc_srt_len[ist]))
    #        for ida in DA:
    #            ax[imp]['track'].plot(lead_t,MAEs[ist][imp][ida]['track'],color=colorset[ida],marker=marker_type[ist],markersize=5,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4)
    #            ax[imp]['mslp'].plot(lead_t,MAEs[ist][imp][ida]['mslp'],color=colorset[ida],marker=marker_type[ist],markersize=5,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4)

    # Plot mean over storms
    mae_sts = MAE_overStorms()
    lead_t = list(range(0,max(fc_srt_len.values())))
    alphas = {'CONV':0.5,'IR':0.7,'MW':0.9}
    for imp in MP:
        for ida in DA:
            ax[imp]['track'].plot(lead_t,mae_sts[imp][ida]['track'],color=colorset[ida],linewidth=3)
            ax[imp]['mslp'].plot(lead_t,mae_sts[imp][ida]['mslp'],color=colorset[ida],linewidth=3)

    ## Option 2: DA by color; Storm by marker
    #lines = ax['WSM6']['track'].get_lines()
    #lgd_0 = Storms + ['Mean']
    #legend = ax['WSM6']['track'].legend([lines[i] for i in [0,3,6,9,-3]], lgd_0,fontsize='8',loc='upper center',ncol=2)
    # Add the first legend manually to the current Axes
    #ax['WSM6']['track'].add_artist(legend)

    lines = ax['WSM6']['track'].get_lines()
    lgd_0 = DA
    legend0 = ax['WSM6']['track'].legend([lines[i] for i in [-3,-2,-1]], lgd_0,fontsize='8',loc='upper center',ncol=3)
    #legend0.set_alpha( 0.5 )
    # Add the first legend manually to the current Axes
    ax['WSM6']['track'].add_artist(legend0)

    fig.text(0.32,0.95,'Mean Track Error', fontsize=12, ha='center', va='center',rotation='horizontal')
    fig.text(0.78,0.95,'Mean MSLP Error', fontsize=12, ha='center', va='center',rotation='horizontal')
    fig.text(0.03,0.73,'WSM6', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.03,0.32,'THO', fontsize=10, ha='center', va='center',rotation='vertical')

    # Set axis attributes
    ax['WSM6']['track'].set_ylim( [0,0.5] )
    ax['THO']['track'].set_ylim( [0,0.5] )
    ax['WSM6']['mslp'].set_ylim( [0,0.5] )
    ax['THO']['mslp'].set_ylim( [0,0.5] )
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
        ax[imp]['track'].set_ylabel('Track: Scaled MAEs', fontsize=8)
        if imp == 'WSM6':
            fig.text(0.53,0.72,'MSLP: Scaled MAEs', fontsize=8, ha='center', va='center',rotation='vertical')
            #ax[imp]['mslp'].set_ylabel('MSLP: Scaled MAEs')
        else:
            fig.text(0.53,0.30,'MSLP: Scaled MAEs', fontsize=8, ha='center', va='center',rotation='vertical')

        # grid lines
        ax[imp]['track'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')
        ax[imp]['mslp'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_absError_scaled_MAEs_withFCtime.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

    return None



# ------------------------------------------------------------------------------------------------------
#            Operation: Plot MAEs
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


# Plot stat for individual storm
# layout:
# WSM6_MAE_track, WSM6_MAE_MSLP
# THO_MAE_track, THO_MAE_MSLP
def plot_2by2_MAEs_each():

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
    colorset = {'CONV': 'black','IR':'red','MW':'blue'}
    marker_type = {'HARVEY':'|','IRMA':'x','JOSE':'_','MARIA':'.'}
    line_width = {'HARVEY':1,'IRMA':1,'JOSE':1,'MARIA':1,'Mean':2.5}

    for imp in MP:
        for ist in Storms:
            lead_t = list(range(0,fc_srt_len[ist]))
            for ida in DA:
                ax[imp]['track'].plot(lead_t,MAEs[ist][imp][ida]['track'],color=colorset[ida],marker=marker_type[ist],markersize=5,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4)
                ax[imp]['mslp'].plot(lead_t,MAEs[ist][imp][ida]['mslp'],color=colorset[ida],marker=marker_type[ist],markersize=5,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4)


    lines = ax['WSM6']['mslp'].get_lines()
    lgd_0 = DA
    legend0 = ax['WSM6']['mslp'].legend([lines[i] for i in [0,1,2]], lgd_0,fontsize='8',loc='upper left')
    #legend0.set_alpha( 0.5 )
    # Add the first legend manually to the current Axes
    ax['WSM6']['mslp'].add_artist(legend0)

    fig.text(0.32,0.95,'Track', fontsize=12, ha='center', va='center',rotation='horizontal')
    fig.text(0.78,0.95,'MSLP', fontsize=12, ha='center', va='center',rotation='horizontal')
    fig.text(0.03,0.73,'WSM6', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.03,0.32,'THO', fontsize=10, ha='center', va='center',rotation='vertical')

    # Set axis attributes
    ax['WSM6']['track'].set_ylim( [-0.1,500] )
    ax['THO']['track'].set_ylim( [-0.1,500] )
    ax['WSM6']['mslp'].set_ylim( [-0.1,45] )
    ax['THO']['mslp'].set_ylim( [-0.1,90] )
    #ax['THO']['mslp'].set_ylim( [-0.1,45] )

    # y ticks
    wsm6_yticks = list(range(0,500+1,100))
    ax['WSM6']['track'].set_yticks( wsm6_yticks )
    ax['WSM6']['track'].set_yticklabels( [str(it) for it in wsm6_yticks] )
    tho_yticks = list(range(0,500+1,100))
    ax['THO']['track'].set_yticks( tho_yticks)
    ax['THO']['track'].set_yticklabels( [str(it) for it in tho_yticks] )

    wsm6_yticks = list(range(0,45+1,5))
    ax['WSM6']['mslp'].set_yticks( wsm6_yticks )
    ax['WSM6']['mslp'].set_yticklabels( [str(it) for it in wsm6_yticks] )
    tho_yticks = list(range(0,90+1,15))
    #tho_yticks = list(range(0,45+1,5))
    ax['THO']['mslp'].set_yticks( tho_yticks)
    ax['THO']['mslp'].set_yticklabels( [str(it) for it in tho_yticks] )
    # Bold ticks for WSM6 range
    for tick in ax['THO']['mslp'].get_yticklabels():
        if int(tick.get_text()) >= 45:
            tick.set_fontweight('bold')
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
        ax[imp]['track'].set_ylabel('Track: MAEs (km)')
        ax[imp]['mslp'].set_ylabel('MSLP: MAEs (hPa)')
        # grid lines
        ax[imp]['track'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')
        ax[imp]['mslp'].grid(True,linewidth=1, color='gray', alpha=0.3, linestyle='-')

    # Save figure
    des_name = small_dir+'/SYSTEMS/Vis_analyze/Paper1/sys_fc_HPI_absError_MAEs_withFCtime.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )

    return None

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
    colorset = {'CONV': 'black','IR':'red','MW':'blue'}
    marker_type = {'HARVEY':'|','IRMA':'x','JOSE':'_','MARIA':'.'}
    line_width = {'HARVEY':1,'IRMA':1,'JOSE':1,'MARIA':1,'Mean':2.5}

    # Plot individual simulations
    subax_track = add_subplot_axes(ax['THO']['track'],[0.1,0.66,0.45,0.3]) # small subplot
    subax_mslp = add_subplot_axes(ax['THO']['mslp'],[0.1,0.66,0.45,0.3]) # small subplot

    for imp in MP:
        for ist in Storms:
            lead_t = list(range(0,fc_srt_len[ist]))
            for ida in DA:
                ax[imp]['track'].plot(lead_t,MAEs[ist][imp][ida]['track'],color=colorset[ida],marker=marker_type[ist],markersize=5,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4) 
                ax[imp]['mslp'].plot(lead_t,MAEs[ist][imp][ida]['mslp'],color=colorset[ida],marker=marker_type[ist],markersize=5,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4)
                if ist == 'JOSE' and imp == "THO":
                    if ida == 'CONV':
                        pass
                    elif ida == 'MW':
                        pass
                    else:
                        subax_track.plot(lead_t,MAEs[ist][imp][ida]['track_all'],color=colorset[ida],marker=marker_type[ist],markersize=3,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4)
                        subax_track.set_xlim( [-0.1,max(fc_srt_len.values())-1] )
                        subax_track.set_yticks( [0,450,900,1350] )
                        subax_track.set_ylim([-0.1,1350.1])
                        subax_track.set_yticklabels(['0','450','900','1350'],fontsize=6)
                        subax_track.set_xticks([0,4,8,12,16] )
                        subax_track.set_xticklabels(['D0','D1','D2','D3','D4'])
                        subax_track.grid(True,linewidth=0.8, color='gray', alpha=0.3, linestyle='-')

                        subax_mslp.plot(lead_t,MAEs[ist][imp][ida]['mslp_all'],color=colorset[ida],marker=marker_type[ist],markersize=3,linewidth=line_width[ist],markeredgewidth=2,alpha=0.4)
                        subax_mslp.set_xlim( [-0.1,max(fc_srt_len.values())-1] )
                        subax_mslp.set_yticks( [0,45,90] )
                        subax_mslp.set_ylim([-0.1,90])
                        subax_mslp.set_yticklabels(['0','45','90'],fontsize=6)
                        subax_mslp.set_xticks([0,4,8,12,16] )
                        subax_mslp.set_xticklabels(['D0','D1','D2','D3','D4'])
                        subax_mslp.grid(True,linewidth=0.8, color='gray', alpha=0.3, linestyle='-')

    # Plot mean over storms
    #mae_sts = MAE_overStorms()
    #lead_t = list(range(0,max(fc_srt_len.values())))
    #alphas = {'CONV':0.5,'IR':0.7,'MW':0.9}
    #for imp in MP:
    #    for ida in DA:
    #        ax[imp]['track'].plot(lead_t,mae_sts[imp][ida]['track'],color=colorset[ida],linewidth=3)
    #        ax[imp]['mslp'].plot(lead_t,mae_sts[imp][ida]['mslp'],color=colorset[ida],linewidth=3)


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
    #lgd_0 = Storms + ['Mean']
    lgd_0 = Storms
    #legend = ax['WSM6']['track'].legend([lines[i] for i in [0,3,6,9,-3]], lgd_0,fontsize='8',loc='upper center',ncol=2)
    legend = ax['WSM6']['track'].legend([lines[i] for i in [0,3,6,9,]], lgd_0,fontsize='8',loc='upper center',ncol=2)
    # Add the first legend manually to the current Axes
    ax['WSM6']['track'].add_artist(legend)

    #lines = ax['WSM6']['mslp'].get_lines()
    legend_colors = list(colorset.values())
    lgd_1 = DA
    proxy_artists = [plt.Line2D([0], [0], color=color, lw=2) for color in legend_colors ]
    legend1 = ax['WSM6']['mslp'].legend(proxy_artists,lgd_1,fontsize='8',loc='upper left')
    #legend0 = ax['WSM6']['mslp'].legend([lines[i] for i in [-3,-2,-1]], lgd_0,fontsize='8',loc='upper left')
    #legend0.set_alpha( 0.5 )
    # Add the first legend manually to the current Axes
    ax['WSM6']['mslp'].add_artist(legend1)

    fig.text(0.32,0.95,'Track Error', fontsize=12, ha='center', va='center',rotation='horizontal')
    fig.text(0.78,0.95,'MSLP Error', fontsize=12, ha='center', va='center',rotation='horizontal')
    fig.text(0.03,0.73,'WSM6', fontsize=10, ha='center', va='center',rotation='vertical')
    fig.text(0.03,0.32,'THO', fontsize=10, ha='center', va='center',rotation='vertical')

    # Set axis attributes
    ax['WSM6']['track'].set_ylim( [-0.1,500] )
    ax['THO']['track'].set_ylim( [-0.1,500] )
    ax['WSM6']['mslp'].set_ylim( [-0.1,50] )
    ax['THO']['mslp'].set_ylim( [-0.1,100] )
    # y ticks
    wsm6_yticks = list(range(0,500+1,100))
    ax['WSM6']['track'].set_yticks( wsm6_yticks )
    ax['WSM6']['track'].set_yticklabels( [str(it) for it in wsm6_yticks] )
    tho_yticks = list(range(0,500+1,100))
    ax['THO']['track'].set_yticks( tho_yticks)
    ax['THO']['track'].set_yticklabels( [str(it) for it in tho_yticks] )

    wsm6_yticks = list(range(0,50+1,5))
    ax['WSM6']['mslp'].set_yticks( wsm6_yticks )
    ax['WSM6']['mslp'].set_yticklabels( [str(it) for it in wsm6_yticks] )
    tho_yticks = list(range(0,100+1,10))
    ax['THO']['mslp'].set_yticks( tho_yticks)
    ax['THO']['mslp'].set_yticklabels( [str(it) for it in tho_yticks] )
    # Bold ticks for WSM6 range
    for tick in ax['THO']['mslp'].get_yticklabels():
        if int(tick.get_text()) >= 50:
            tick.set_fontweight('bold')
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
        ax[imp]['track'].set_ylabel('Track: MAEs (km)')
        if imp == 'WSM6':
            ax[imp]['mslp'].set_ylabel('MSLP: MAEs (hPa)')
        else:
            fig.text(0.53,0.30,'MSLP: MAEs (hPa)', fontsize=8, ha='center', va='center',rotation='vertical')

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
    colorset = {'CONV': 'black','IR':'red','MW':'blue'}

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


def test_plot():

    matplotlib.rcParams['font.size'] = 8

    # Calculate MAE
    MAE = calc_MAE_overTime()

    # Set up figure
    fig,ax = plt.subplots( figsize=(6.5,3.25),dpi=200) # standard: 6.5,8.5
    
    # Customization
    colorset = {'CONV':'black','IR':'red','MW':'blue'}
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

# Instead average MAEs of storms into a MAEs, this function concatenate these MAEs
def MAE_addStorms( MAEs ):

    mae_add = {}
    for imp in MP:
        mae_add[imp] = {}
        for ida in DA:
            mae_add[imp][ida] = {}
            # initialize
            add_track = MAEs[Storms[0]][imp][ida]['track']
            add_mslp = MAEs[Storms[0]][imp][ida]['mslp']
            add_vmax = MAEs[Storms[0]][imp][ida]['vmax']
            # loop thru storms:
            for ist in Storms[1:]:
                add_track = np.concatenate( (add_track,MAEs[ist][imp][ida]['track'] ),axis=0)
                add_mslp = np.concatenate( (add_mslp,MAEs[ist][imp][ida]['mslp'] ),axis=0)
                add_vmax = np.concatenate( (add_vmax,MAEs[ist][imp][ida]['vmax'] ),axis=0)
            # load the dictionary
            mae_add[imp][ida]['track'] = add_track
            mae_add[imp][ida]['mslp'] = add_mslp
            mae_add[imp][ida]['vmax'] = add_vmax
    
    return mae_add



def stat_test( mae_sts ):

    var = ['track','mslp']
    pairs = [['CONV','IR'],['IR','MW'],['CONV','MW']]
    for imp in MP:
        print('MPS: '+imp)
        for iv in var:
            print('')
            print('MPS-var: '+imp+'-'+iv)
            for ip in pairs:
                print(ip)
                if if_anderson:
                    print('Anderson test...')
                    for ipp in ip:
                        result = anderson(mae_sts[imp][ipp][iv])
                        print("Anderson-Darling Test Statistic:", result.statistic)
                        print("Critical Values:", result.critical_values)
                        print("Significance Levels:", result.significance_level)

                if independent_t:
                    print('T test...')
                    Result = ttest_ind(mae_sts[imp][ip[0]][iv], mae_sts[imp][ip[1]][iv], equal_var=True)
                    print( Result )
                    if Result.pvalue < 0.05:
                        print("The errors are significantly different!")

                if pair_t:
                    print('Paired-T test...')
                    Result = ttest_rel(mae_sts[MP[0]][ip[0]][iv], mae_sts[MP[1]][ip[1]][iv] )
                    print( Result )
                    if Result.pvalue < 0.01:
                    print("The errors are significantly different!")

                if if_wilcoxon:
                    print('Wilcoxon test...')
                    Result = wilcoxon(mae_sts[imp][ip[0]][iv], mae_sts[imp][ip[1]][iv])
                    print( Result )
                    if Result.pvalue < 0.01:
                        print("The errors are significantly different!")




if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/Clean_results/'

    #--------Configuration------------
    Storms =  ['HARVEY','IRMA','JOSE','MARIA']
    DA = ['CONV','IR','MW']
    MP = ['WSM6','THO'] #

    # if operate over the same number of samples for all forecasts
    sameNum_sample = False
    if sameNum_sample:
        fc_run_hrs = 60

    MAE_only = False
    MAE_wrt_time = True
    if_scale = True
    if if_scale:
        min_max_scale = True
        z_score = False
    # Significant test
    independent_t = False
    pair_t = True
    if_wilcoxon = False
    if_anderson = False # check if the data is normal distribution
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

    # Calculate absolute error at each synoptic time for a forecast
    Exper_error = calculate_abs_error()

    # Plot
    if MAE_only:
        test_plot()

    if MAE_wrt_time:

        # Assemble the MAE array
        aerr = Assemble_err_wrt_time( Exper_error )

        if not if_scale:
            # Calculate MAE with same samples
            fc_srt_len, MAEs = calculate_MAEs_wrt_time( aerr )
            mae_sts = MAE_overStorms()

            # Perform test
            mae_add = MAE_addStorms( MAEs )
            stat_test( mae_sts )

            #plot_MSLP_MAE()
            #plot_2by2_MAEs()
            #plot_2by2_MAEs_each()
        else:
            aerr = Assemble_err_wrt_time( Exper_error )
            fc_srt_len, MAEs = calculate_scaled_MAEs_wrt_time( aerr )
            mae_sts = MAE_overStorms()

            # Perform test
            mae_add = MAE_addStorms( MAEs )
            stat_test( mae_sts )

            #plot_2by2_scaled_MAEs()


    #plot_4by2():unable to plot it with 6.5 inch wide and 8.5 inch tall
    # layout:
    # WSM6_normMAE, WSM6_FSP
    # THO_normMAE, THO_FSP
    # WSM6_MAE_track, WSM6_MAE_MSLP
    # THO_MAE_track, THO_MAE_MSLP











import os
import glob
import numba
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib import ticker
import math
from datetime import datetime, timedelta
import time
from fast_histogram import histogram2d as hist2d
import time
import matplotlib.colors as mcolors
from scipy import interpolate
import matplotlib.path as mpath 
import matplotlib.patches as mpatches 

import Util_data as UD
import Util_Vis

# Generate time series
def generate_times( Storms, start_time_str, end_time_str, interval ):

    dict_times = {}
    for ist in Storms:
        time_diff = datetime.strptime(end_time_str[ist],"%Y%m%d%H%M") - datetime.strptime(start_time_str[ist],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[ist],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+interval, interval))]
        dict_times[ist] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times

# Read variables at obs resolution/location for one experiment 
def read_Tbs_obsRes_oneExper(big_dir,istorm,imp,ida,exp_name,d_times,sensor):

    Hx_dir = big_dir+istorm+'/'+exp_name+'/Obs_Hx/IR/'
    dict_allTb = {}

    Yo_all = []
    meanYb_all = []
    meanYa_all = []

    for DAtime in d_times[istorm]:

        Tb_file = Hx_dir+DAtime+'/mean_obs_res_d03_' + DAtime + '_' +  sensor + '.txt'
        Yo_obs = []
        meanYb_obs = []
        meanYa_obs = []

        # Read records
        print('Reading ', Tb_file)
        with open(Tb_file) as f:
            next(f)
            all_lines = f.readlines()

        for line in all_lines:
            split_line = line.split()
            Yo_obs.append( float(split_line[3]) )
            meanYb_obs.append( float(split_line[4]) )
            meanYa_obs.append( float(split_line[5]) )

        # add up records
        Yo_all.extend( Yo_obs )
        meanYb_all.extend( meanYb_obs )
        meanYa_all.extend( meanYa_obs )

    dict_allTb['All_times'] = {'Yo_obs':Yo_all, 'meanYb_obs':meanYb_all, 'meanYa_obs':meanYa_all}

    return dict_allTb

def combine_storms_allTimes( Storms, Exper_Tb ):

    # Initialize
    allStorm_tb = {}
    for imp in MP:
        allStorm_tb[imp] = {}

    # Gather data
    for imp in MP:
        all_Yo = []
        all_meanYb = []
        all_meanYa = []

        for istorm in Storms:
            tb_case = Exper_Tb[istorm][imp]['All_times']
            all_Yo.extend( tb_case['Yo_obs'] )
            all_meanYb.extend( tb_case['meanYb_obs'] )
            all_meanYa.extend( tb_case['meanYa_obs'] )

        allStorm_tb[imp]['Yo_obs'] = np.array( all_Yo )
        allStorm_tb[imp]['meanYb_obs'] = np.array( all_meanYb )
        allStorm_tb[imp]['meanYa_obs'] = np.array( all_meanYa )

    return allStorm_tb



# Calculate the conditional mean of y for each bin of x
def cal_conditional_mean( x,y,xedges ): 

    conditional_means = []

    for i in range(len(xedges) - 1):
        # Find all OMB values within this bin of predictor
        in_bin = (x >= xedges[i]) & (x < xedges[i+1])
        if np.any(in_bin):
            mean_omb = np.nanmean(y[in_bin])
        else:
            mean_omb = np.nan  # If no points in this bin
        conditional_means.append(mean_omb)

    conditional_means = np.array(conditional_means)
    return conditional_means 

def OMB_predictor_2D_bins( Storms_Tb,imp,ipd,omb_val ):

    omb_pdf = {}
    omb_pdf['global_mean'] = np.nanmean( omb_val )

    # Count
    if ipd == 'obs':
        x = Storms_Tb[imp]['Yo_obs']
    elif ipd == 'hxb':
        x = Storms_Tb[imp]['meanYb_obs']
    elif ipd == 'obs+hxb':
        x = 0.5*(Storms_Tb[imp]['meanYb_obs']+Storms_Tb[imp]['Yo_obs'])

    tmp = hist2d(x,omb_val,range=[[min_Tb_rg,max_Tb_rg],[min_Tbdiff_rg,max_Tbdiff_rg]],bins=number_bins)
    #meanYa_Yo = hist2d(Storms_Tb['Yo_obs'],Storms_Tb['meanYa_obs'],range=[[min_Tb_rg,max_Tb_rg],[min_Tb_rg,max_Tb_rg]],bins=number_bins)
    # Compute log10 of non-zero elements, return NAN where elements are ZERO
    with np.errstate(divide='ignore', invalid='ignore'):
        omb_pdf['2d'] = np.where(tmp != 0, np.log10(tmp), np.nan)
        #dcount['meanYa_Yo'] = np.where(meanYa_Yo != 0, np.log10(meanYa_Yo), np.nan)
    # Calculate the conditional mean
    xedges = np.linspace(min_Tb_rg, max_Tb_rg, number_bins + 1)
    xcenter = 0.5 * (xedges[:-1] + xedges[1:])
    cm = cal_conditional_mean( x,omb_val,xedges )
    omb_pdf['x_center'] = xcenter
    omb_pdf['con_mean'] = cm


    return omb_pdf

# Calculate A and c
def Prepare_for_bc( Storms_Tb,imp,ipd,n_order ):

    # omb
    omb = Storms_Tb[imp]['Yo_obs']-Storms_Tb[imp]['meanYb_obs']

    # Idetify samples of the predictor
    if ipd == 'obs':
        p = Storms_Tb[imp]['Yo_obs']
    elif ipd == 'hxb':
        p = Storms_Tb[imp]['meanYb_obs']
    elif ipd == 'obs+hxb':
        p = 0.5*(Storms_Tb[imp]['meanYb_obs']+Storms_Tb[imp]['Yo_obs'])

    # Find the c (define it to be the mean of the predictor value)
    c = np.nanmean( p )

    # Fill in the matrix A
    if n_order == 0:
        A = np.ones( (len(p),1) )        
    else:
        A = np.zeros( (len(p),n_order+1) )
        A[:] = np.nan
        for i in range( len(p) ):
            for j in range( n_order+1):
                A[i,j] = (p[i] - c)**j

    d_forBC = {'omb':omb,'p_vals':p,'c':c,'A':A}
    return d_forBC

# Find the BC coefficients for the nonlinear relationship
# Eqn 11 in Otkin et al., 2018
def find_BCbeta( d_forBC,imp,ipd,N_order ):

    A = d_forBC['A']
    omb = d_forBC['omb']
    # Check NaNs or Infs
    #print('Check nan value for A and omb:',np.any(np.isnan(A)), np.any(np.isnan(omb)))
    #print('Check infinite for A and omb:',np.any(np.isinf(A)), np.any(np.isinf(omb)))
    #np.linalg.cond(A)  # Large condition number (e.g. > 1e10) indicates instability

    # Keep non-nan samples
    print('Cleaning data... Kick out nan samples...')
    ## Step 1: Remove rows from A and b where b has NaN
    mask_omb = ~np.isnan(omb)   # b is (m,1); flatten or index column
    A = A[mask_omb]
    omb = omb[mask_omb]
    ## Remove rows from A and b where any element in a row of A is NaN
    mask_A = ~np.isnan(A).any(axis=1)  # True for rows with no NaNs
    A_clean = A[mask_A]
    omb_clean = omb[mask_A]

    # Calculate the beta
    alpha = 1e-9
    I = np.eye(N_order+1)
    beta_opt = np.linalg.inv(alpha*I + A_clean.T @ A_clean) @ A_clean.T @ omb_clean
    #print( beta_opt )

    return beta_opt

# Bias-correct OMB
def bc_omb( imp,ipd,d_forBC,bc_opt,n_order ):

    # omb
    omb = d_forBC['omb']
    #p = d_forBC['p_vals'] # sample values for one predictor
    A = d_forBC['A']

    omb_BCed = omb - A @ bc_opt[imp][ipd][:n_order+1]
    return omb_BCed

# ------------------------------------------------------------------------------------------------------
#            Operation: Plot
# ------------------------------------------------------------------------------------------------------
def set_figure_background_2by3():

    # Set up figure
    fig = plt.figure( figsize=(6.5,5),dpi=200) # standard: 6.5,8.5
    grids = fig.add_gridspec(ncols=3,nrows=2,top=0.95,left=0.1,hspace=0.03,bottom=0.2,)

    ax = {}
    # gridspec inside gridspec
    for imp in MP:
        ax[imp]= {}
        for ipd in preds:
            ir = preds.index( ipd )
            ax[imp][ipd] = fig.add_subplot( grids[MP.index(imp),ir] )

    return fig,ax

def set_figure_text_2by3( fig,ax,show ):

    # Add color bar below the plot
    cbar_ax = fig.add_axes([0.10, 0.1, 0.78, 0.015])
    color_bar = fig.colorbar(show,cax=cbar_ax,orientation='horizontal')#ticks=bounds)
    color_bar.ax.xaxis.set_major_locator(ticker.FixedLocator([0,1,2,3,4]))
    color_bar.ax.xaxis.set_major_formatter(ticker.FixedFormatter(['$10^0$', '$10^1$', '$10^2$', '$10^3$','$10^4$']))
    color_bar.ax.tick_params(labelsize=10)
    color_bar.ax.set_xlabel('Number of Samples (#)',fontsize=8)

    # ticks and labels
    xtick = list(range(min_Tb_rg,max_Tb_rg+1,10))
    ytick = list(range(min_Tbdiff_rg,max_Tbdiff_rg+1,10))
    for imp in MP:
        for ipd in preds:
            # ticks
            ax[imp][ipd].set_xlim(xmin=min_Tb_rg,xmax=max_Tb_rg)
            ax[imp][ipd].set_ylim(ymin=min_Tbdiff_rg,ymax=max_Tbdiff_rg)
            ax[imp][ipd].set_xticks( xtick )
            ax[imp][ipd].set_yticks( ytick )
            ax[imp][ipd].set_xticklabels( [str(it) for it in xtick],fontsize=8,rotation=30, )
            ax[imp][ipd].set_yticklabels( [str(it) for it in ytick],fontsize=9 )
            # lables
            ax[imp][ipd].set_xlabel('Tbs (K)',fontsize='9')

    # Add sub title
    for ipd in preds:
        if preds.index(ipd) == 0 and ipd == 'obs':
            fig.text(0.21,0.943,'obs', fontsize=12, ha='center', va='center')
        elif preds.index(ipd) == 1 and ipd == 'hxb':
            fig.text(0.50,0.943,r'$\overline{H(X_b)}$', fontsize=12, ha='center', va='center')
        elif preds.index(ipd) == 2 and ipd == 'obs+hxb':
            fig.text(0.78,0.943,r'$\frac{1}{2}\left(obs+\overline{H(X_b)}\right)$', fontsize=12, ha='center', va='center')

    # Add y label
    fig.text(0.02,0.76,'WSM6', fontsize=10, ha='center', va='center' ,rotation='vertical')
    fig.text(0.02,0.38,'THO', fontsize=10, ha='center', va='center',rotation='vertical')


# Plot 2D histogram of OMB dependent on predictors and conditional bias
# row-microphysics: WSM6, THO
# column-predictors: observation, hxb, symmetric
def plot_2by3_OMB( pred_omb ):

    fig,ax = set_figure_background_2by3()

    # Customize the colormap
    color_intervals = [0,0.5,1,1.5,2,2.5,3,3.5,4.0]
    exist_cmap = plt.cm.jet.reversed()
    colors = exist_cmap(np.linspace(0,1,len(color_intervals)-1))
    new_map = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals)-1)

    # Plot 2d historgram
    xedges = np.linspace(min_Tb_rg, max_Tb_rg, number_bins + 1)
    for imp in MP:
        print('Global mean of OMB for',imp,pred_omb[imp]['global_mean'])
        for ipd in preds:
            show = ax[imp][ipd].imshow(np.transpose(pred_omb[imp][ipd]['2d']),cmap=new_map,
                        origin='lower',vmin=0,vmax=4,extent=[min_Tb_rg,max_Tb_rg,min_Tbdiff_rg,max_Tbdiff_rg],aspect='equal')
            for i in range(len(xedges) - 1):
                ax[imp][ipd].hlines(y=pred_omb[imp][ipd]['con_mean'][i], xmin=xedges[i], xmax=xedges[i+1], colors='black', linewidth=1)
            # plot line where y=0 is
            ax[imp][ipd].axhline( y = 0, color='grey', linestyle='-', linewidth=0.8)
            # Plot the global mean of bias
            ax[imp][ipd].axhline( y=pred_omb[imp]['global_mean'], color='magenta', linestyle='--', linewidth=0.8)

    # Add details
    set_figure_text_2by3( fig,ax,show )

    # Add sub y label
    fig.text(0.05,0.76,'OMB (K)', fontsize=8, ha='center', va='center' ,rotation='vertical')
    fig.text(0.05,0.38,'OMB (K)', fontsize=8, ha='center', va='center' ,rotation='vertical')

    # Save
    des_name = 'IR_2dPDF_OMB_'+str(len(dict_times[Storms[0]]))+'cycles_'+DA+'.png'
    plt.savefig( des_name )
    print( 'Saving the figure: ', des_name )
    plt.close()

# Plot 2D histogram of bias-corrected OMB dependent on predictors and conditional bias
# row-microphysics: WSM6, THO
# column-predictors: observation, hxb, symmetric
def plot_2by3_OMB_BCed( pred_omb_bced,n_order ):

    fig,ax = set_figure_background_2by3()

    # Customize the colormap
    color_intervals = [0,0.5,1,1.5,2,2.5,3,3.5,4.0]
    exist_cmap = plt.cm.jet.reversed()
    colors = exist_cmap(np.linspace(0,1,len(color_intervals)-1))
    new_map = mcolors.LinearSegmentedColormap.from_list('custom_colormap',colors,N=len(color_intervals)-1)

    # Plot 2d historgram
    xedges = np.linspace(min_Tb_rg, max_Tb_rg, number_bins + 1)
    for imp in MP:
        #print('Global mean of OMB for',imp,pred_omb_bced[imp]['global_mean'])
        for ipd in preds:
            show = ax[imp][ipd].imshow(np.transpose(pred_omb_bced[imp][ipd]['2d']),cmap=new_map,
                        origin='lower',vmin=0,vmax=4,extent=[min_Tb_rg,max_Tb_rg,min_Tbdiff_rg,max_Tbdiff_rg],aspect='equal')
            for i in range(len(xedges) - 1):
                ax[imp][ipd].hlines(y=pred_omb_bced[imp][ipd]['con_mean'][i], xmin=xedges[i], xmax=xedges[i+1], colors='black', linewidth=1)
            # Plot the global mean of bias
            ax[imp][ipd].axhline( y=pred_omb_bced[imp][ipd]['global_mean'], color='magenta', linestyle='--', linewidth=0.8)

    # Add details
    set_figure_text_2by3( fig,ax,show )

    # Add sub y label
    fig.text(0.05,0.76,'OMB_BCed (K)', fontsize=8, ha='center', va='center' ,rotation='vertical')
    fig.text(0.05,0.38,'OMB_BCed (K)', fontsize=8, ha='center', va='center' ,rotation='vertical')

    # Save
    if n_order == 0:
        des_name = 'IR_2dPDF_OMB_BCconst_'+str(len(dict_times[Storms[0]]))+'cycles_'+DA+'.png'
    elif n_order == 1:
        des_name = 'IR_2dPDF_OMB_BC1st_'+str(len(dict_times[Storms[0]]))+'cycles_'+DA+'.png'
    elif n_order == 2:
        des_name = 'IR_2dPDF_OMB_BC2nd_'+str(len(dict_times[Storms[0]]))+'cycles_'+DA+'.png'
    elif n_order == 3:
        des_name = 'IR_2dPDF_OMB_BC3rd_'+str(len(dict_times[Storms[0]]))+'cycles_'+DA+'.png'
    elif n_order == 4:
        des_name = 'IR_2dPDF_OMB_BC4th_'+str(len(dict_times[Storms[0]]))+'cycles_'+DA+'.png'


    plt.savefig( des_name )
    print( 'Saving the figure: ', des_name )
    plt.close()

# ------------------------------------------------------------------------------------------------------
#            Operation: Workflow to correct bias 
# ------------------------------------------------------------------------------------------------------

def correct_bias( n_order,Storms_Tb,bc_opt ):
    # Make bins
    pred_omb_bced = {}
    for imp in MP:
        pred_omb_bced[imp] = {}
        print('OMB minus constant bias for ',imp,' experiment...')
        for ipd in preds:
            # Calculate polynomials and other parameters
            d_forBC = Prepare_for_bc( Storms_Tb,imp,ipd,n_order )
            # Make correction on OMB
            omb_val = bc_omb( imp,ipd,d_forBC,bc_opt,n_order )
            # Make bins
            pred_omb_bced[imp][ipd] = OMB_predictor_2D_bins( Storms_Tb,imp,ipd,omb_val )
    # Plot
    plot_2by3_OMB_BCed( pred_omb_bced, n_order)


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Clean_Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/Clean_results/'

    #--------Configuration------------
    Storms = ['IRMA','JOSE','MARIA','HARVEY']#['HARVEY','IRMA','JOSE','MARIA']
    DA = 'CONV'
    MP = ['WSM6','THO']
    sensor = 'abi_gr'

    start_time_str = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}
    end_time_str = {'HARVEY':'201708231200','IRMA':'201709040000','JOSE':'201709060000','MARIA':'201709170000'}
    Consecutive_times = True

    # Predictor for OMB 
    preds = ['obs','hxb','obs+hxb']

    # Order of the nonlinear 
    N_order = 4

    # Range
    min_Tb_rg = 180
    max_Tb_rg = 260
    min_Tbdiff_rg = -40
    max_Tbdiff_rg = 40

    # number of bins
    scale_factor = 1 # 1: 1k per bin; 2: 0.5k per bin
    number_bins = (max_Tb_rg-min_Tb_rg)*scale_factor

    # Actions
    minimize_sbeta = True # least square fitting to find the optimal value for BC coefficient   
    omb_original = False
    
    omb_BCconst = True
    omb_BC1st = True
    omb_BC2nd = True
    omb_BC3rd = True
    omb_BC4th = True

    #------------------------------------

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm]= {}
        for imp in MP:
            Exper_names[istorm][imp] = UD.generate_one_name( istorm,DA,imp )

    # Identify DA times in the period of interest
    dict_times = generate_times( Storms, start_time_str, end_time_str, 1 )

    # Read obs, Hxb, and Hxa of all files
    Exper_Tb = {}
    for istorm in Storms:
        Exper_Tb[istorm] = {}
        for imp in MP:
            Exper_Tb[istorm][imp] = read_Tbs_obsRes_oneExper(big_dir,istorm,imp,DA,Exper_names[istorm][imp],dict_times,sensor)

    # Combine data for all storms
    Storms_Tb = combine_storms_allTimes( Storms, Exper_Tb )

    # PDF for OMB
    if omb_original:
        # Make bins
        pred_omb = {}
        for imp in MP:
            pred_omb[imp] = {}
            omb_val = Storms_Tb[imp]['Yo_obs']-Storms_Tb[imp]['meanYb_obs']
            print('Predicting OMB based on predictors for ',imp,' experiment...')
            for ipd in preds:
                pred_omb[imp][ipd] = OMB_predictor_2D_bins( Storms_Tb,imp,ipd,omb_val )
        # Plot
        plot_2by3_OMB( pred_omb ) #OMB

    # Find the BC coefficients
    if minimize_sbeta:
        bc_opt = {}
        for imp in MP:
            bc_opt[imp] = {}
            for ipd in preds:
                d_forBC = Prepare_for_bc( Storms_Tb,imp,ipd,N_order )
                bc_opt[imp][ipd] = find_BCbeta( d_forBC,imp,ipd,N_order ) 

    # Remove constant value
    if omb_BCconst:
        n_order = 0
        correct_bias( n_order,Storms_Tb,bc_opt ) 

    # Remove first order of bias
    if omb_BC1st:
        n_order = 1
        correct_bias( n_order,Storms_Tb,bc_opt ) 

    # Remove 2nd order of bias
    if omb_BC2nd:
        n_order = 2
        correct_bias( n_order,Storms_Tb,bc_opt )

    # Remove 3rd order of bias
    if omb_BC2nd:
        n_order = 3
        correct_bias( n_order,Storms_Tb,bc_opt )

    # Remove 4th order of bias
    if omb_BC2nd:
        n_order = 4
        correct_bias( n_order,Storms_Tb,bc_opt )

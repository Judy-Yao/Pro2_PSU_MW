
import glob
import numpy as np
import re
from matplotlib import pyplot as plt
from datetime import datetime

import Util_data as UD
import Obspace_compare_IR_txt_bin as IR_obs


def RMSE(simu, obs):
    return np.sqrt( ((simu - obs) ** 2).mean() )

def Bias(simu, obs):
    return  np.sum((simu - obs),0)/np.size(obs,0)

def mean_Yo_Hx(simu, obs):
    return  np.sum((obs - simu),0)/np.size(obs,0)


# Read variables at obs resolution/location
def read_EnsTb( Tb_files ):

    dict_IR = {}
    dict_IR['Ya'] = {}
    dict_IR['Yb'] = {}

    print('Reading the ensemble...')
    for ifile in Tb_files:

        Yb_obs = []
        Ya_obs = []

        with open( ifile ) as f:
            next(f)
            all_lines = f.readlines()
        for line in all_lines:
            split_line = line.split()
            Yb_obs.append( float(split_line[4]) )
            Ya_obs.append( float(split_line[5]) )
        
        # Use regex to extract the member id
        match = re.search(r"mem(\d{3})", ifile)
        id_ = match.group(1) if match else None
        
        dict_IR['Yb'][id_] =  np.array( Yb_obs )
        dict_IR['Ya'][id_] =  np.array( Ya_obs )

    return dict_IR


# calculate the bias for the whole ensemble
def Ens_bias( simus ,obs ):

    bias = []
    for key in simus.keys():
        bias.append( Bias(simus[key],obs) )
    return bias

# calculate the rmse for the whole ensemble
def Ens_rmse( simus ,obs ):

    rmse = []
    for key in simus.keys():
        rmse.append( RMSE(simus[key],obs) )
    return rmse

def IR_metric_one_cycle( ist,imp ):

    bias = {}
    bias['mean'] = {}
    bias['mem'] = {}
    rmse = {}
    rmse['mean'] = {}
    rmse['mem'] = {}

    # Get metrics for the CONV experiment and IR obs
    conv_name = UD.generate_one_name( ist,'CONV',imp )
    # conv mean
    conv_mean_file = big_dir+ist+'/'+conv_name+'/Obs_Hx/IR/'+DAtime[ist]+"/mean_obs_res_d03_"+DAtime[ist]+'_'+ sensor+'.txt'
    d_conv = IR_obs.read_Tb_obsRes(conv_mean_file, sensor)

    bias['mean']['xb'] = Bias(d_conv['meanYb_obs'], d_conv['Yo_obs'] )
    bias['mean']['conv'] = Bias(d_conv['meanYa_obs'], d_conv['Yo_obs'] )
    rmse['mean']['xb'] = RMSE(d_conv['meanYb_obs'], d_conv['Yo_obs'] )
    rmse['mean']['conv'] = RMSE(d_conv['meanYa_obs'], d_conv['Yo_obs'] )
    # conv for members
    conv_mem_files = sorted( glob.glob( big_dir+ist+'/'+conv_name+'/Obs_Hx/IR/'+DAtime[ist]+'/Interp_Tb_mem0*.txt') )
    d_conv_mem = read_EnsTb( conv_mem_files )

    bias['mem']['xb'] = Ens_bias( d_conv_mem['Yb'], d_conv['Yo_obs'] )
    bias['mem']['conv'] = Ens_bias( d_conv_mem['Ya'], d_conv['Yo_obs'] )
    rmse['mem']['xb'] = Ens_rmse( d_conv_mem['Yb'], d_conv['Yo_obs'] )
    rmse['mem']['conv'] = Ens_rmse( d_conv_mem['Ya'], d_conv['Yo_obs'] )

    # Get metrics for the IR experiment 
    ir_name = UD.generate_one_name( ist,'IR',imp )
    ir_file = big_dir+ist+'/'+ir_name+'/Obs_Hx/IR/'+DAtime[ist]+"/mean_obs_res_d03_"+DAtime[ist]+'_'+ sensor+'.txt'
    d_ir = IR_obs.read_Tb_obsRes(ir_file, sensor )

    bias['mean']['ir'] = Bias(d_ir['meanYa_obs'], d_ir['Yo_obs'] )
    rmse['mean']['ir'] = RMSE(d_ir['meanYa_obs'], d_ir['Yo_obs'] )
    # ir for members
    ir_mem_files = sorted( glob.glob( big_dir+ist+'/'+ir_name+'/Obs_Hx/IR/'+DAtime[ist]+'/Interp_Tb_mem0*.txt') )
    d_ir_mem = read_EnsTb( ir_mem_files )
    bias['mem']['ir'] = Ens_bias( d_ir_mem['Ya'], d_conv['Yo_obs'] )
    rmse['mem']['ir'] = Ens_rmse( d_ir_mem['Ya'], d_conv['Yo_obs'] )


    return bias,rmse

# Plot figures
def Plot_metric():

    fig,ax = plt.subplots(1, 2, dpi=300, figsize=(10,7)) #8) ) #figsize=(8,8)
    colors = {'WSM6':'red','THO':'blue'}
    lines = {'HARVEY':'solid','IRMA':'dotted','JOSE':'dashed','MARIA':'dashdot'}
    #x_label = ['Xb:mean','Xb:HARVEY','Xb:IRMA','Xb:JOSE','Xb:MARIA',
    #     'Xb-CONV:mean','Xa-CONV:HARVEY','Xa-CONV:IRMA','Xa-CONV:JOSE','Xa-CONV:MARIA',
    #     'Xb-IR:mean','Xa-IR:HARVEY','Xa-IR:IRMA','Xa-IR:JOSE','Xa-IR:MARIA',        
    #     ]
    x_label = ['Xb:HARVEY','Xb:IRMA','Xb:mean','Xb:JOSE','Xb:MARIA',
               'Xa-CONV:HARVEY','Xa-CONV:IRMA','Xa-CONV:mean','Xa-CONV:JOSE','Xa-CONV:MARIA',
               'Xa-IR:HARVEY','Xa-IR:IRMA','Xa-IR:mean','Xa-IR:JOSE','Xa-IR:MARIA', ]


    for ist in Storms:
        
        for imp in MP:
            x = [2,7,12]  # Keys will be used for x-axis
            y_bias = d_bias[ist][imp]['mean']
            y = [y_bias['xb'],y_bias['conv'],y_bias['ir']]
            ax[0].plot(x, y, linestyle=lines[ist], color=colors[imp], marker='s', linewidth=3)
            y_rmse = d_rmse[ist][imp]['mean']
            y = [y_rmse['xb'],y_rmse['conv'],y_rmse['ir']]
            ax[1].plot( x,y,linestyle=lines[ist],color=colors[imp],marker='s',linewidth=3 )

    # Plot ensembles for Xb
    for ist in Storms:
        for imp in MP:
            x_idx = x_label.index( 'Xb:'+ist )
            y = d_bias[ist][imp]['mem']['xb']
            ax[0].scatter( [x_idx]*len(y),y,s=2,color=colors[imp],alpha=0.5)
            y = d_rmse[ist][imp]['mem']['xb']
            ax[1].scatter( [x_idx]*len(y),y,s=2,color=colors[imp],alpha=0.5)

    # Plot ensembles for Xa-CONV
    for ist in Storms:
        for imp in MP:
            x_idx = x_label.index( 'Xa-CONV:'+ist )
            y = d_bias[ist][imp]['mem']['conv']
            ax[0].scatter( [x_idx]*len(y),y,s=2,color=colors[imp],alpha=0.8)
            y = d_rmse[ist][imp]['mem']['conv']
            ax[1].scatter( [x_idx]*len(y),y,s=2,color=colors[imp],alpha=0.8)

    # Plot ensembles for Xa-CONV
    for ist in Storms:
        for imp in MP:
            x_idx = x_label.index( 'Xa-IR:'+ist )
            y = d_bias[ist][imp]['mem']['ir']
            ax[0].scatter( [x_idx]*len(y),y,s=2,color=colors[imp],alpha=0.8)
            y = d_rmse[ist][imp]['mem']['ir']
            ax[1].scatter( [x_idx]*len(y),y,s=2,color=colors[imp],alpha=0.8)

    # lines to separate sections
    for i in range(2):
        ax[i].axvline( x = 4.5, color='gray', linestyle='-',linewidth=1,alpha=0.5)
        ax[i].axvline( x = 9.5, color='gray', linestyle='-',linewidth=1,alpha=0.5)
    # y = 0 line
    ax[0].axhline(y=0, color='gray', linestyle='-',linewidth=1,alpha=0.5)

    # Add storm label
    H_coords = [0,5,10]
    y0_coords = [-9,-9,-9]
    y1_coords = [4.5,4.5,4.5]
    
    for x,y in zip(H_coords,y0_coords):
        ax[0].text( x,y,'H',fontsize='10',ha='center', fontweight='bold')
    for x,y in zip(H_coords,y1_coords):
        ax[1].text( x,y,'H',fontsize='10',ha='center', fontweight='bold') 

    I_coords = [1,6,11]
    for x,y in zip(I_coords,y0_coords):
        ax[0].text( x,y,'I',fontsize='10',ha='center', fontweight='bold') 
    for x,y in zip(I_coords,y1_coords):
        ax[1].text( x,y,'I',fontsize='10',ha='center', fontweight='bold') 

    J_coords = [3,8,13]
    for x,y in zip(J_coords,y0_coords):
        ax[0].text( x,y,'J',fontsize='10',ha='center', fontweight='bold') 
    for x,y in zip(J_coords,y1_coords):
        ax[1].text( x,y,'J',fontsize='10',ha='center', fontweight='bold') 

    M_coords = [4,9,14]
    for x,y in zip(M_coords,y0_coords):
        ax[0].text( x,y,'M',fontsize='10',ha='center', fontweight='bold') 
    for x,y in zip(M_coords,y1_coords):
        ax[1].text( x,y,'M',fontsize='10',ha='center', fontweight='bold') 

    # labels
    ax[0].set_ylim([-10,10])
    ax[1].set_ylim([4,18])
    for i in range(2):
        ax[i].set_xticks([2,7,12])
        ax[i].tick_params(axis='both', which='major', labelsize=13) 
        ax[i].set_xticklabels(['Xb','Xa:CONV','Xa:CONV+IR'])
    ax[0].set_ylabel('Bias (K)',fontsize=15) 
    ax[1].set_ylabel('RMSE (K)',fontsize=15)

    # legend: bias
    proxy_artists = [plt.Line2D([0], [0], color='red', linestyle=lw) for lw in list(lines.values()) ]
    ax[0].legend(proxy_artists,list(lines.keys()),fontsize='12',loc='upper center',ncol=2)
    # legend: rmse
    proxy_artists = [plt.Line2D([0], [0], color=color) for color in list(colors.values()) ]
    ax[1].legend(proxy_artists,list(colors.keys()),fontsize='12',loc='upper center',ncol=2)

    # title
    ax[0].set_title( 'IR Bias: mean{'+r'$\mathbf{\overline{H(X)}}$'+' - Obs}',fontsize=15 )
    ax[1].set_title( 'IR RMSE: mean{('+r'$\mathbf{\overline{H(X)}}$'+' - Obs)**2}',fontsize=15 )
    fig.suptitle('1st WRF-EnKF cycle',fontsize=18)

    des_name = small_dir+'/SYSTEMS/Vis_analyze/Tb/IR_bias_rmse_1stCycle.png'
    plt.savefig( des_name )
    print( 'Saving the figure to '+des_name )


if __name__ == '__main__':


    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir =  '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ---------- Configuration -------------------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    MP = ['THO','WSM6']

    sensor = 'abi_gr'
    Num_ens = 60

    DAtime = {'HARVEY':'201708221200','IRMA':'201709030000','JOSE':'201709050000','MARIA':'201709160000'}

    # limitations
    limit = False
    # ------------------------------------------------------

    # Calculate metrics
    d_bias = {}
    d_rmse = {}
    for ist in Storms:
        d_bias[ist] = {}
        d_rmse[ist] = {}
        for imp in MP:
            bias,rmse = IR_metric_one_cycle( ist,imp )
            d_bias[ist][imp] = bias
            d_rmse[ist][imp] = rmse

    # Plot IR metric
    Plot_metric()




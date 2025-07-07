
import os,fnmatch # functions for interacting with the operating system
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
#from matplotlib. import gridspec
from copy import deepcopy

def plot_truth_obs(  ):

    # Generate x and y values
    x_vals = np.linspace(0, 60, 300)
    ir_vals = perfect_CloudtoIR(x_vals)
    w_vals = perfect_CloudtoWSD(x_vals)

    # Plot
    fig = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    gs = GridSpec(1, 1, figure=fig)
    ax = plt.subplot(gs[0,0])

    ax.axvline(x=xt,color='red',linestyle='--',label='Observation')
    #---------------------------------
    # Plot on left y-axis (IR)
    ax.plot(x_vals, ir_vals, color='grey', linewidth=2, zorder=0)

    # truth
    ax.scatter(xt,yt_ir-0.5,s=70,marker='s',color='#b4285f') ##e58a8a 

    # observation
    ax.scatter(0.1,yo_ir,s=100,marker='*',color='#b4285f')
    ax.text(0.8, yo_ir-0.5, 'obs', fontsize=15,color='#b4285f')

    ax.text(12, 210, 'Perfect Model: Cloud and IR',
        fontsize=12, rotation=-10,
        bbox=dict(facecolor='lightgray', edgecolor='black', boxstyle='round',alpha=0.3))

    ax.set_xlim(0,30)
    ax.set_ylim(190,250)
    #ax.set_yticks(range(190, 250, 5))
    ax.set_xlabel('Cloud (gram m$^{-2}$)',fontsize=12)
    ax.set_ylabel('IR Tbs (K)',fontsize=12)

    #---------------------------------
    # Plot on right y-axis (Wind speed)
    ax1 = ax.twinx()
    ax1.plot(x_vals, w_vals, color='grey', linewidth=2, zorder=1)

    # truth
    ax1.scatter(xt,yt_wsd,s=70,marker='s',color='#1c1b78') #8b8ae5

    # observation
    ax1.scatter(29.9,yo_wsd,s=100,marker='*',color='#1c1b78')
    ax1.text(27, yo_wsd-0.5, 'obs', fontsize=15,color='#1c1b78')

    ax1.text(12, 30, 'Perfect Model: Cloud and Wind',
        fontsize=12,rotation=30,
        bbox=dict(facecolor='lightgray', edgecolor='black', boxstyle='round',alpha=0.3))

    ax1.set_ylim(0,55)
    ax1.set_yticks(range(0, 56, 5))
    ax1.set_ylabel('10-meter Wind (m/s)',fontsize=12)
    #plt.tight_layout()
 
    ax.set_title('Truth and Observations: yt_ir = '+str(yt_ir)+' K',fontsize=15)
    plt.savefig("truth_obs.png", dpi=300, )
    print('truth_obs.png')

def plot_background():

    # Calculate means
    xb_mean = np.mean( xb_ens )
    hxb_ir_mean = np.mean( hxb_ens_ir )
    hxb_wsd_mean = np.mean( hxb_ens_wsd )

    # Generate x and y values
    x_vals = np.linspace(0, 60, 300)
    ir_vals = perfect_CloudtoIR(x_vals)
    w_vals = perfect_CloudtoWSD(x_vals)

    # Plot
    fig = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    gs = GridSpec(1, 1, figure=fig)
    ax = plt.subplot(gs[0,0])

    ax.axvline(x=xt,color='red',linestyle='--',label='Observation')
    #---------------------------------
    # Plot on left y-axis (IR)
    ax.plot(x_vals, ir_vals, color='grey', linewidth=2, zorder=0)

    # truth
    ax.scatter(xt,yt_ir-0.5,s=70,marker='s',color='#b4285f') ##e58a8a 

    # observation
    ax.scatter(0.1,yo_ir,s=100,marker='*',color='#b4285f')
    ax.text(0.8, yo_ir-0.5, 'obs', fontsize=15,color='#b4285f')

    # background ensemble in model space
    ax.scatter(xb_ens,[190.5 for i in xb_ens],s=12,marker='o',facecolors='none',edgecolors='black',)
   
    # background ensemble
    ax.scatter(xb_ens,hxb_ens_ir,s=15,marker='o',color='#b4285f',)
    # background ensemble in IR space
    ax.scatter([0.2 for i in xb_ens],hxb_ens_ir,s=12,marker='o',facecolors='none',edgecolors='#b4285f',) 

    # simulated ensemble mean
    ax.scatter(xb_mean,hxb_ir_mean,s=50,marker='o',color='#b4285f')
    ax.text(4, 245, 'H(X) to IR Tbs', fontsize=15,color='#b4285f')
    #ax.text(20, 225, 'Perfect Model: Cloud and Wind',
    #    fontsize=12,rotation=30,
    #    bbox=dict(facecolor='lightgray', edgecolor='black', boxstyle='round',alpha=0.3))

    ax.set_xlim(0,30)
    ax.set_ylim(190,250)
    #ax.set_yticks(range(190, 250, 5))
    ax.set_xlabel('Cloud (gram m$^{-2}$)',fontsize=12)
    ax.set_ylabel('IR Tbs (K)',fontsize=12)

    #---------------------------------
    # Plot on right y-axis (Wind speed)
    ax1 = ax.twinx()
    ax1.plot(x_vals, w_vals, color='grey', linewidth=2, zorder=1)

    # truth
    ax1.scatter(xt,yt_wsd,s=70,marker='s',color='#1c1b78') #8b8ae5

    # observation
    ax1.scatter(29.9,yo_wsd,s=100,marker='*',color='#1c1b78')
    ax1.text(27, yo_wsd-0.5, 'obs', fontsize=15,color='#1c1b78')

    # background ensemble
    ax1.scatter(xb_ens,hxb_ens_wsd,s=15,marker='o',color='#1c1b78',)
    # background ensemble in IR space
    ax1.scatter([29.8 for i in xb_ens],hxb_ens_wsd,s=12,marker='o',facecolors='none',edgecolors='#1c1b78',)

    # simulated ensemble mean
    ax1.scatter(xb_mean,hxb_wsd_mean,s=50,marker='o',color='#1c1b78')
    ax.text(4, 242, 'H(X) to Wind', fontsize=15,color='#1c1b78')

    #ax1.text(20, 12, 'Perfect Model: Cloud and IR',
    #    fontsize=12, rotation=-5,
    #    bbox=dict(facecolor='lightgray', edgecolor='black', boxstyle='round',alpha=0.3))

    ax1.set_ylim(0,55)
    ax1.set_yticks(range(0, 56, 5))
    ax1.set_ylabel('10-meter Wind (m/s)',fontsize=12)
    #plt.tight_layout()

    ax.set_title('Ensemble background',fontsize=15)
    plt.savefig("background.png", dpi=300, )
    print('background.png')

# show prior and posterior during EnSRF
def plot_only_assimilate_yo_wsd( d_da ):
   
    # Calculate means
    xb_mean = np.mean( xb_ens )
    hxb_wsd_mean = np.mean( hxb_ens_wsd )

    # Generate x and y values
    x_vals = np.linspace(0, 60, 300)
    ir_vals = perfect_CloudtoIR(x_vals)
    w_vals = perfect_CloudtoWSD(x_vals)

    # Plot
    fig1 = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    gs = GridSpec(1, 1, figure=fig1)
    ax = plt.subplot(gs[0,0])

    #---------------------------------
    # Plot on right y-axis (Wind speed)
    ax1 = ax.twinx()
    ax1.plot(x_vals, w_vals, color='grey', linewidth=2, zorder=1)

    # truth
    #ax1.scatter(xt,yt_wsd,s=70,marker='s',color='#1c1b78')

    # observation
    ax1.scatter(29.9,yo_wsd,s=100,marker='*',color='#1c1b78')
    ax1.text(27, yo_wsd-0.5, 'obs', fontsize=15,color='#1c1b78')
    ax1.axhline(y=yo_wsd,color='#1c1b78',linestyle='--',)

    # background ensemble
    ax1.scatter(xb_ens,hxb_ens_wsd,s=15,marker='o',color='#1c1b78',)
    # background ensemble mean
    ax1.scatter(xb_mean,hxb_wsd_mean,s=50,marker='o',facecolors='#1c1b78',edgecolors='black',)
    # background ensemble in model space
    ax1.scatter(xb_ens,[0.5 for i in xb_ens],s=12,marker='o',facecolors='none',edgecolors='#1c1b78',)

    ax1.set_xlim(0,30)
    ax1.set_ylim(0,55)
    ax1.set_yticks(range(0, 56, 5))
    ax1.set_ylabel('10-meter Wind (m/s)',fontsize=12)
    ax.set( yticks=[] )
    ax.set_xlabel('Cloud (gram m$^{-2}$)',fontsize=12)

    ax.set_title('Before EnSRF Assimilates ONE Observation',fontsize=15)
    plt.savefig("before_assmi_only_oneobs.png", dpi=300, )
    print('before_assmi_only_oneobs.png')

    ####################################################################
    # Second plot
    fig2 = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    gs = GridSpec(1, 1, figure=fig2)
    ax = plt.subplot(gs[0,0])

    #---------------------------------
    # Plot on right y-axis (Wind speed)
    ax1 = ax.twinx()
    ax1.plot(x_vals, w_vals, color='grey', linewidth=2, zorder=1)

    # truth
    #ax1.scatter(xt,yt_wsd,s=70,marker='s',color='#1c1b78')

    # observation
    ax1.scatter(29.9,yo_wsd,s=100,marker='*',color='#1c1b78')
    ax1.text(27, yo_wsd-0.5, 'obs', fontsize=15,color='#1c1b78')
    ax1.axhline(y=yo_wsd,color='#1c1b78',linestyle='--',)

    # background ensemble
    ax1.scatter(xb_ens,hxb_ens_wsd,s=15,marker='o',color='#1c1b78',)
    # background ensemble mean
    ax1.scatter(xb_mean,hxb_wsd_mean,s=50,marker='o',facecolors='#1c1b78',edgecolors='black')
    # background ensemble in model space
    ax1.scatter(xb_ens,[0.5 for i in xb_ens],s=12,marker='o',facecolors='none',edgecolors='#1c1b78',)

    # analysis ensemble in model space
    ax1.scatter(d_da['1st_xa_ens'],[0.5 for i in d_da['1st_xa_ens']],s=12,marker='o',facecolors='none',edgecolors='#8b8ae5',)

    ax1.set_xlim(0,30)
    ax1.set_ylim(0,55)
    ax1.set_yticks(range(0, 56, 5))
    ax1.set_ylabel('10-meter Wind (m/s)',fontsize=12)
    ax.set( yticks=[] )
    ax.set_xlabel('Cloud (gram m$^{-2}$)',fontsize=12)

    ax.set_title('After EnSRF Assimilates ONE Observation',fontsize=15)
    plt.savefig("after_assmi_only_oneobs.png", dpi=300, )
    print('after_assmi_only_oneobs.png')

    ####################################################################
    # Third plot
    fig3 = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    gs = GridSpec(1, 1, figure=fig3)
    ax = plt.subplot(gs[0,0])

    #---------------------------------
    # Plot on right y-axis (Wind speed)
    ax1 = ax.twinx()
    ax1.plot(x_vals, w_vals, color='grey', linewidth=2, zorder=1)

    # truth
    #ax1.scatter(xt,yt_wsd,s=70,marker='s',color='#1c1b78')

    # observation
    ax1.scatter(29.9,yo_wsd,s=100,marker='*',color='#1c1b78')
    ax1.text(27, yo_wsd-0.5, 'obs', fontsize=15,color='#1c1b78')
    ax1.axhline(y=yo_wsd,color='#1c1b78',linestyle='--',)

    # background ensemble
    ax1.scatter(xb_ens,hxb_ens_wsd,s=15,marker='o',color='#1c1b78',)
    # background ensemble mean
    ax1.scatter(xb_mean,hxb_wsd_mean,s=50,marker='o',facecolors='#1c1b78',edgecolors='black')
    # background ensemble in model space
    ax1.scatter(xb_ens,[0.5 for i in xb_ens],s=12,marker='o',facecolors='none',edgecolors='#1c1b78',)

    # analysis ensemble in model space
    ax1.scatter(d_da['1st_xa_ens'],[0.5 for i in d_da['1st_xa_ens']],s=12,marker='o',color='#8b8ae5',)
    # analysis ensemble
    hxa_ens_wsd  = imperfect_CloudtoWSD( d_da['1st_xa_ens'] )
    hxa_wsd_mean = np.mean( hxa_ens_wsd )
    ax1.scatter(d_da['1st_xa_ens'],hxa_ens_wsd,s=15,marker='o',color='#8b8ae5',)
    # analysis ensemble mean
    ax1.scatter(d_da['1st_xa_mean'],hxa_wsd_mean,s=50,marker='o',facecolors='#8b8ae5',edgecolors='black')

    ax1.set_xlim(0,30)
    ax1.set_ylim(0,55)
    ax1.set_yticks(range(0, 56, 5))
    ax1.set_ylabel('10-meter Wind (m/s)',fontsize=12)
    ax.set( yticks=[] )
    ax.set_xlabel('Cloud (gram m$^{-2}$)',fontsize=12)

    # Add a wider upward arrow
    #ax1.arrow(x=12, y=2, dx=0, dy=20, width=1, head_length=1, fc='k', ec='k')
    ax1.annotate(
        '', xy=(12, 28), xytext=(12, 2),
        arrowprops=dict(facecolor='black', width=6, headwidth=12,alpha=0.8)
    )
    ax1.text(14, 20, 'H(X) to Wind', fontsize=12,color='black')

    ax.set_title('H(X) on Updated Model-Space Ensemble',fontsize=15)
    plt.savefig("after_assmi_only_oneobs_perform_Hx.png", dpi=300, )
    print('after_assmi_only_oneobs_perform_Hx.png')

# 
def plot_1st_assimilate_yo_wsd(  d_da ):

    # Calculate means
    xb_mean = np.mean( xb_ens )
    hxb_wsd_mean = np.mean( hxb_ens_wsd )

    # Generate x and y values
    x_vals = np.linspace(0, 60, 300)
    ir_vals = perfect_CloudtoIR(x_vals)
    w_vals = perfect_CloudtoWSD(x_vals)

    # Plot
    fig1 = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    gs = GridSpec(1, 1, figure=fig1)
    ax = plt.subplot(gs[0,0])

    # observation
    ax.scatter(0.1,yo_ir,s=100,marker='*',color='#b4285f')
    ax.text(0.8, yo_ir-0.5, 'obs', fontsize=15,color='#b4285f')

    ax.set_ylim(190,250)
    ax.set_xlabel('Cloud (gram m$^{-2}$)',fontsize=12)
    ax.set_ylabel('IR Tbs (K)',fontsize=12)
    #---------------------------------
    # Plot on right y-axis (Wind speed)
    ax1 = ax.twinx()
    ax1.plot(x_vals, w_vals, color='grey', linewidth=2, zorder=1)

    # truth
    #ax1.scatter(xt,yt_wsd,s=70,marker='s',color='#1c1b78')

    # observation
    ax1.scatter(29.9,yo_wsd,s=100,marker='*',color='#1c1b78')
    ax1.text(27, yo_wsd-0.5, 'obs', fontsize=15,color='#1c1b78')
    ax1.axhline(y=yo_wsd,color='#1c1b78',linestyle='--',)

    # background ensemble
    ax1.scatter(xb_ens,hxb_ens_wsd,s=15,marker='o',color='#1c1b78',)
    # background ensemble mean
    ax1.scatter(xb_mean,hxb_wsd_mean,s=50,marker='o',color='#1c1b78',)
    ax.text(4, 242, 'H(X) to Wind', fontsize=15,color='#1c1b78')
    # background ensemble in model space
    ax1.scatter(xb_ens,[0.5 for i in xb_ens],s=12,marker='o',facecolors='none',edgecolors='#1c1b78',)

    ax1.set_xlim(0,30)
    ax1.set_ylim(0,55)
    ax1.set_yticks(range(0, 56, 5))
    ax1.set_ylabel('10-meter Wind (m/s)',fontsize=12)

    ax.set_title('Before EnSRF Assimilates 1st Observation',fontsize=15)
    plt.savefig("step1.png", dpi=300, )
    print('step1.png')

    # observation
    ax.scatter(0.1,yo_ir,s=100,marker='*',color='#b4285f')
    ax.text(0.8, yo_ir-0.5, 'obs', fontsize=15,color='#b4285f')
    ####################################################################
    # Second plot
    fig2 = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    gs = GridSpec(1, 1, figure=fig2)
    ax = plt.subplot(gs[0,0])

    #---------------------------------
    # Plot on right y-axis (Wind speed)
    ax1 = ax.twinx()
    #ax1.plot(x_vals, w_vals, color='grey', linewidth=2, zorder=1)

    # truth
    #ax1.scatter(xt,yt_wsd,s=70,marker='s',color='#1c1b78')

    # observation
    ax1.scatter(29.9,yo_wsd,s=100,marker='*',color='#1c1b78')
    ax1.text(27, yo_wsd-0.5, 'obs', fontsize=15,color='#1c1b78')

    # background ensemble
    #ax1.scatter(xb_ens,hxb_ens_wsd,s=15,marker='o',color='#1c1b78',)
    # background ensemble mean
    #ax1.scatter(xb_mean,hxb_wsd_mean,s=50,marker='o',facecolors='#1c1b78',edgecolors='black')
    # background ensemble in model space
    ax1.scatter(xb_ens,[0.5 for i in xb_ens],s=12,marker='o',facecolors='none',edgecolors='#1c1b78',)

    # analysis ensemble in model space
    ax1.scatter(d_da['1st_xa_ens'],[1.5 for i in d_da['1st_xa_ens']],s=12,marker='o',facecolors='none',edgecolors='#8b8ae5',)

    ax1.set_ylim(0,55)
    ax1.set_yticks(range(0, 56, 5))
    ax1.set_ylabel('10-meter Wind (m/s)',fontsize=12)

    #---------------------------------
    # Plot on left y-axis (IR Tbs)
    ax.plot(x_vals, ir_vals, color='grey', linewidth=2, zorder=0)
    ax.scatter(0.1,yo_ir,s=100,marker='*',color='#b4285f')
    ax.text(0.8, yo_ir+1, 'obs', fontsize=15,color='#b4285f')
    ax.axhline(y=yo_ir,color='#b4285f',linestyle='--',)
    #print(d_da)
    # EnSRF-updated Hx in IR
    ax.scatter(d_da['1st_xa_ens'],d_da['1st_hxa_ens'][1,:],s=15,marker='o',facecolors='#b4285f',edgecolors='black')
    # analysis ensemble mean
    ax.scatter(d_da['1st_xa_mean'],d_da['1st_hxa_mean'][1],s=50,marker='o',facecolors='#b4285f',edgecolors='black',)
    ax.text(4, 245, 'EnSRF-updated simulated obs', fontsize=15,color='#b4285f')
    # analysis ensemble in IR space
    ax.scatter([0.2 for i in xb_ens],d_da['1st_hxa_ens'][1,:],s=12,marker='o',facecolors='none',edgecolors='#b4285f',)

    ax.set_xlim(0,30)
    ax.set_ylim(190,250)
    ax.set_xlabel('Cloud (gram m$^{-2}$)',fontsize=12)
    ax.set_ylabel('IR Tbs (K)',fontsize=12)

    ax.set_title('After EnSRF Assimilates 1st Observation',fontsize=15)
    plt.savefig("step2.png", dpi=300, )
    print('step2.png')

    ####################################################################
    # Third plot
    # Plot
    fig3 = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    gs = GridSpec(1, 1, figure=fig3)
    ax = plt.subplot(gs[0,0])

   #---------------------------------
    # Plot on right y-axis (Wind speed)
    ax1 = ax.twinx()
    #ax1.plot(x_vals, w_vals, color='grey', linewidth=2, zorder=1)

    # truth
    #ax1.scatter(xt,yt_wsd,s=70,marker='s',color='#1c1b78')

    # observation
    ax1.scatter(29.9,yo_wsd,s=100,marker='*',color='#1c1b78')
    ax1.text(27, yo_wsd-0.5, 'obs', fontsize=15,color='#1c1b78')

    # analysis ensemble in model space
    ax1.scatter(d_da['1st_xa_ens'],[0.5 for i in d_da['1st_xa_ens']],s=12,marker='o',facecolors='none',edgecolors='#b4285f',)

    ax1.set_ylim(0,55)
    ax1.set_yticks(range(0, 56, 5))
    ax1.set_ylabel('10-meter Wind (m/s)',fontsize=12)
    #---------------------------------
    # Plot on left y-axis (IR Tbs)
    ax.plot(x_vals, ir_vals, color='grey', linewidth=2, zorder=0)
    ax.scatter(0.1,yo_ir,s=100,marker='*',color='#b4285f')
    ax.text(0.8, yo_ir+1, 'obs', fontsize=15,color='#b4285f')
    ax.axhline(y=yo_ir,color='#b4285f',linestyle='--',)

    # EnSRF-updated Hx in IR
    ax.scatter(d_da['1st_xa_ens'],d_da['1st_hxa_ens'][1,:],s=15,marker='o',facecolors='#b4285f',edgecolors='black')
    # background (analysis after 1st update) ensemble mean
    ax.scatter(d_da['1st_xa_mean'],d_da['1st_hxa_mean'][1],s=50,marker='o',facecolors='#b4285f',edgecolors='black')
    ax.text(4, 245, 'EnSRF-updated simulated obs', fontsize=15,color='#b4285f')
    # background (analysis after 1st update) ensemble in IR space
    ax.scatter([0.2 for i in xb_ens],d_da['1st_hxa_ens'][1,:],s=12,marker='o',facecolors='none',edgecolors='#b4285f',)

    ax.set_xlim(0,30)
    ax.set_ylim(190,250)
    ax.set_xlabel('Cloud (gram m$^{-2}$)',fontsize=12)
    ax.set_ylabel('IR Tbs (K)',fontsize=12)

    ax.set_title('Before EnSRF Assimilates 2nd Observation',fontsize=15)
    plt.savefig("step3.png", dpi=300, )
    print('step3.png')


def plot_2nd_assimilate_yo_ir( d_da ):

   # Generate x and y values
    x_vals = np.linspace(0, 60, 300)
    ir_vals = perfect_CloudtoIR(x_vals)
    #w_vals = perfect_CloudtoWSD(x_vals)

    # Plot
    fig1 = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    gs = GridSpec(1, 1, figure=fig1)
    ax = plt.subplot(gs[0,0])

   #---------------------------------
    # Plot on right y-axis (Wind speed)
    ax1 = ax.twinx()
    #ax1.plot(x_vals, w_vals, color='grey', linewidth=2, zorder=1)

    # truth
    #ax1.scatter(xt,yt_wsd,s=70,marker='s',color='#1c1b78')

    # observation
    ax1.scatter(29.9,yo_wsd,s=100,marker='*',color='#1c1b78')
    ax1.text(27, yo_wsd-0.5, 'obs', fontsize=15,color='#1c1b78')

    # analysis ensemble in model space
    ax1.scatter(d_da['1st_xa_ens'],[0.5 for i in d_da['1st_xa_ens']],s=12,marker='o',facecolors='none',edgecolors='#b4285f',)

    ax1.set_ylim(0,55)
    ax1.set_yticks(range(0, 56, 5))
    ax1.set_ylabel('10-meter Wind (m/s)',fontsize=12)
    #---------------------------------
    # Plot on left y-axis (IR Tbs)
    # truth
    ax.scatter(xt,yt_ir-0.5,s=70,marker='s',color='#b4285f')

    ax.plot(x_vals, ir_vals, color='grey', linewidth=2, zorder=0)
    ax.scatter(0.1,yo_ir,s=100,marker='*',color='#b4285f')
    ax.text(0.8, yo_ir+1, 'obs', fontsize=15,color='#b4285f')
    ax.axhline(y=yo_ir,color='#b4285f',linestyle='--',)
   
    # analysis ensemble in model space
    hxa_ens_ir  = imperfect_CloudtoIR( d_da['2nd_xa_ens'] )
    hxa_ir_mean = np.mean( hxa_ens_ir )
    ax1.scatter(d_da['2nd_xa_ens'],[1.5 for i in d_da['2nd_xa_ens']],s=12,marker='o',facecolors='none',edgecolors='#e58a8a',)
    ax.scatter(d_da['2nd_xa_mean'],hxa_ir_mean,s=50,marker='o',facecolors='#e58a8a',edgecolors='#e58a8a',)
    ax.text(4, 245, 'H(X) to IR Tbs', fontsize=15,color='#e58a8a')
    # background (analysis after 1st update) ensemble mean
    ax.scatter(d_da['1st_xa_mean'],d_da['1st_hxa_mean'][1],s=50,marker='o',facecolors='#b4285f',edgecolors='black')
    ax.text(4, 242, 'EnSRF-updated simulated obs', fontsize=15,color='#b4285f')

    # EnSRF-updated Hx in IR
    #ax.scatter(d_da['1st_xa_ens'],d_da['1st_hxa_ens'][1,:],s=15,marker='o',color='#b4285f',)
    # background (analysis after 1st update) ensemble mean
    #ax.scatter(d_da['1st_xa_mean'],d_da['1st_hxa_mean'][1],s=50,marker='o',facecolors='#b4285f',edgecolors='black')
    # background (analysis after 1st update) ensemble in IR space
    #ax.scatter([0.2 for i in xb_ens],d_da['1st_hxa_ens'][1,:],s=12,marker='o',facecolors='none',edgecolors='#b4285f',)

    ax.set_xlim(0,30)
    ax.set_ylim(190,250)
    ax.set_xlabel('Cloud (gram m$^{-2}$)',fontsize=12)
    ax.set_ylabel('IR Tbs (K)',fontsize=12)

    ax.set_title('After EnSRF Assimilates 2nd Observation',fontsize=15)
    plt.savefig("step4.png", dpi=300, )
    print('step4.png')


# Calculate the dd^
def calculate_var_inno_single( y, mean_xb ):

    return ( y - mean_xb )**2

# Generate N samples
def generate_ensemble():

    # File path for storing the sequence
    file_path = "random_X1Y2/generate_ensemble.npy"

    # Check if the file already exists
    if os.path.exists(file_path):
        # Load the existing sequence
        random_seq = np.load(file_path)
        print("Loaded existing sequence: "+ file_path)
    else:
        # Set a fixed seed
        np.random.seed(5)
        random_seq = np.random.normal(loc=0, scale=sigma_ens_xb, size=Nens)
        np.save(file_path, random_seq)
        print("Generated and saved new sequence:"+ file_path)

    xb_ens = ens_mean_xb + random_seq

    return xb_ens

# Calculate HBH^ 
def calculate_HBH( hx_ens, hx_mean ):

    sum_var = np.sum(( hx_ens - hx_mean ) ** 2)
    return sum_var/(Nens-1)



# Perfect model mapping between cloud amount and IR Tbs
def perfect_CloudtoIR( x ):

    hx = [ 50*np.exp(-0.1304*ix)+195 for ix in x]
    return hx

def perfect_IRtoCloud( hx ):

    if isinstance( hx,(int,float) ):
        x = -1 / 0.1304 * np.log((hx - 195) / 50)
    else:
        x = [-1 / 0.1304 * np.log((iy - 195) / 50) for iy in hx]
    return x

# Perfect model mapping between cloud amount and vertical velocity
def perfect_CloudtoWSD( x ):

    if isinstance( x,(int,float) ):
        hx = 75*( 1 - np.exp( -0.05416*x )) 
    else:
        hx = [ 75*( 1 - np.exp( -0.05416*ix )) for ix in x]
    return hx

def perfect_WSDtoCloud( hx ):

    if isinstance( hx,(int,float) ):
        x = -18.460*np.log( 1- hx/75 )
    else:
        x = [ -18.460*np.log( 1- iy/75 ) for iy in hx]
    return x



# Imperfect model mapping from cloud amount to vertical velocity
def imperfect_CloudtoWSD( x ):

    # File path for storing the sequence
    file_path = "random_X1Y2/imperfect_CloudtoWSD.npy"

    # Check if the file already exists
    if os.path.exists(file_path):
        # Load the existing sequence
        random_seq = np.load(file_path)
        print("Loaded existing sequence:"+ file_path)
    else:
        # Set a fixed seed
        np.random.seed(44)
        random_seq = [np.random.normal(loc=0, scale=sigma_x_inobs) for ix in x]
        np.save(file_path, random_seq)
        print("Generated and saved new sequence:"+ file_path)

    # Map from x to Hx
    hx = []
    for i in range(len(x)):
        hx.append(  75*( 1 - np.exp( -0.05416*x[i] )) + random_seq[i] )

    return hx


# Imperfect model mapping from cloud amount to IR Tbs
def imperfect_CloudtoIR( x ):

    # File path for storing the sequence
    file_path = "random_X1Y2/imperfect_CloudtoIR.npy"

    # Check if the file already exists
    if os.path.exists(file_path):
        # Load the existing sequence
        random_seq = np.load(file_path)
        print("Loaded existing sequence:"+ file_path)
    else:
        # Set a fixed seed
        np.random.seed(42)
        random_seq = [np.random.normal(loc=0, scale=sigma_x_inobs) for ix in x]
        np.save(file_path, random_seq)
        print("Generated and saved new sequence:"+ file_path)

    # Map from x to Hx
    hx = []
    for i in range(len(x)):
        hx.append( 50*np.exp(-0.1304*x[i] )+195 + random_seq[i] )

    return hx



# IR: truth to observation
def IR_yt_to_yo( yt ):

    # File path for storing the sequence
    file_path = "random_X1Y2/IR_yt_to_yo.npy"

    # Check if the file already exists
    if os.path.exists(file_path):
        # Load the existing sequence
        random_seq = np.load(file_path)
        print("Loaded existing sequence:"+ file_path)
    else:
        # Set a fixed seed
        np.random.seed(5)
        random_seq = np.random.normal(loc=0, scale=sigma_ot_ir)
        np.save(file_path, random_seq)
        print("Generated and saved new sequence:"+ file_path)
    
    # truth to observed in IR space
    yo = yt + random_seq

    return yo

# W: truth to observation
def WSD_yt_to_yo( yt ):

    # File path for storing the sequence
    file_path = "random_X1Y2/WSD_yt_to_yo.npy"

    # Check if the file already exists
    if os.path.exists(file_path):
        # Load the existing sequence
        random_seq = np.load(file_path)
        print("Loaded existing sequence:"+ file_path)
    else:
        # Set a fixed seed
        np.random.seed(5)
        random_seq = np.random.normal(loc=0, scale=sigma_ot_wsd)
        np.save(file_path, random_seq)
        print("Generated and saved new sequence:"+ file_path)

    # truth to observed in IR space
    yo = yt + random_seq

    return yo


# Covariance between cloud amount and IR Tbs
def covxy( x,y ):
    
    mean_x = sum(x) / Nens
    mean_y = sum(y) / Nens

    cov_xy = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(Nens)) / (Nens - 1)
    return cov_xy


def EnSRF( ):

    # Create independent copies
    xm = np.mean( xb_ens )
    x = xb_ens - xm # perturbation
    hx = deepcopy( hxb_ens  )
    hx_m = np.mean( hx, axis=1 )

    # Loop thru obs
    for iob in range(len( yo )):
        
        print(' ')
        print('Assimilating obs: '+str(iob+1))
        print('yo: ', yo)
        print('hx_m: ',hx_m)
        print('hx[1]: ',hx[1,:])

        # Innovation
        y_hxm = yo[iob] - hx_m[iob]

        # Variance in modeled observation space
        var_hx = 0
        hx_pert = hx[iob,:] - hx_m[iob] # perturbation
        var_hx = sum( hx_pert ** 2 )

        # Total observation variance
        d = var_hx/(Nens-1) + sigma_ot[iob]**2

        # Alpha used in EnSRF
        alpha = 1/(1+np.sqrt(sigma_ot[iob]**2/d))

        # Compute covariance
        km = 0
        for im in range(Nens):
            km = km + x[im]*hx_pert[im] # covariance !!!!
        km = ( km/(Nens-1) )/d

        #---------------------------------
        # Update model
        #---------------------------------
        # update model member (perturbation!)
        for im in range(Nens):
            x[im] = x[im] - km * alpha * hx_pert[im]

        #print('After: x:', x) 
        # update model mean
        xm = xm + km * y_hxm
        #print('After: xm:', xm) 

        #---------------------------------
        # Deallocate Kalman gain
        #---------------------------------
        km = 0
        #print('Before: hx', hx)
        #print('Before: hx_m', hx_m)

        #---------------------------------
        # Update simulated observations
        #---------------------------------
        if iob == len( yo ) - 1:
            print('No need to update simulated observation ensemble for the next obs!')
        else:

            ya_mean = deepcopy( hx_m[iob] )

            # loop thru next obs
            for inob in range(len(yo)):   # inob: index of next obs
            
                pert_inob = np.zeros( (Nens) )
                cov_obs = 0
                for im in range(Nens):
                    pert_inob[im] = hx[inob,im] - hx_m[inob]
                    cov_obs = cov_obs + hx_pert[im]*( hx[inob,im] - hx_m[inob] )

                ######## current EnSRF: update ensemble member directly
                ## update ensemble member
                #for im in range(Nens):
                #    hx[inob,im] = hx[inob,im] - (cov_obs / (Nens-1) / d)* alpha * hx_pert[im]
                #    hx[inob,im] = hx[inob,im] + cov_obs / (Nens-1) / d * ( yo[iob] - ya_mean )

                ## update mean
                #hx_m[inob] = hx_m[inob] + cov_obs / (Nens-1) / d * ( yo[iob] - ya_mean )

                ####### If update perturbation instead
                # update mean
                hx_m[inob] = hx_m[inob] + cov_obs / (Nens-1) / d * ( yo[iob] - ya_mean )

                # update ensemble member
                for im in range(Nens):
                    hx[inob,im] = (pert_inob[im] - (cov_obs / (Nens-1) / d)* alpha * hx_pert[im]) + hx_m[inob]
            print('After EnSRF update: hx', hx)
            print('After EnSRF update: hx_m', hx_m)

        
        if iob == 0:
            d_da = {'1st_xa_ens':x+xm,'1st_xa_mean':xm,'1st_hxa_ens':hx,'1st_hxa_mean':hx_m}
        elif iob == 1:
            d_da['2nd_xa_ens'] = x+xm
            d_da['2nd_xa_mean'] = xm
            d_da['2nd_hxa_ens'] = hx
            d_da['2nd_hxa_mean'] = hx_m

        if len(yo) == 1:
            plot_only_assimilate_yo_wsd( d_da )
        elif len(yo) == 2 and iob == 0:
            plot_1st_assimilate_yo_wsd( d_da )        
        elif  len(yo) == 2 and iob == 1:
            plot_2nd_assimilate_yo_ir( d_da )


if __name__ == '__main__':

    #------ Default Setup -------
    # number of ensemble members
    Nens = 30
    # truth
    yt_ir = 230  #200                     # truth: IR
    xt = perfect_IRtoCloud( yt_ir ) # truth: cloud
    yt_wsd =  perfect_CloudtoWSD( xt )  # truth: 10-meter wind
    # instrument observation error
    sigma_ot_ir = 3
    sigma_ot_wsd = 2
    sigma_ot = np.array( [sigma_ot_wsd,sigma_ot_ir] )

    # model error (also assume Gaussian)
    sigma_x_inobs = 5 # model error in observation space

    # prior ensemble mean in cloud amount
    ens_mean_xb = 15 #5
    # standard deviation of prior ensemble in cloud amount
    sigma_ens_xb = 3

    # only assimilates one observation or two
    only_oneobs = False

    #------------------------------------

    # Step 1: Make observations from the truth
    yo_ir =  IR_yt_to_yo( yt_ir )
    yo_wsd = WSD_yt_to_yo( yt_wsd )
    if only_oneobs:
        yo = np.array( [yo_wsd, ] )
    else:
        yo = np.array( [yo_wsd, yo_ir] )

    # Step 2: Generate model-state ensemble
    xb_ens = generate_ensemble()
    xb_ens = np.array( xb_ens )

    # Step 3: Simulate observations
    hxb_ens_wsd  = imperfect_CloudtoWSD( xb_ens )
    hxb_ens_ir =  imperfect_CloudtoIR( xb_ens )
    if only_oneobs:
        hxb_ens = np.array( [hxb_ens_wsd,] )
    else:
        hxb_ens = np.array( [hxb_ens_wsd,hxb_ens_ir] )

    # Step 4: Perform the EnSRF
    EnSRF( )

    plot_truth_obs()
    plot_background()
    





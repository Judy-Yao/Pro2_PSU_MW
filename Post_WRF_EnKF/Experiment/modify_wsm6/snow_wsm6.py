#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Equation 6 in HDC
def HDC6():
    
    # Constants
    n0smax = 1e11
    n0s = 2e6 # m-3 m-1
    T0 = 273.15 #K

    Tcold = 223.15
    Tx = np.arange(0, T0-Tcold+1, 1)

    # Factor to N0s (temperature dependent)
    alpha = np.arange(0.05,0.55,0.05)
    alpha = np.insert(alpha,2,0.12)
    itc_fac = np.zeros((len(alpha),len(Tx)))
    for ia in range(len(alpha)):
        itc_fac[ia,:] = np.exp(alpha[ia]*Tx)
    # constraint: higher end
    th_high = n0smax/n0s
    mask = itc_fac >= th_high
    itc_fac[mask] = th_high
    # constraint: lower end
    th_low = 1
    mask = itc_fac <= th_low
    itc_fac[mask] = th_low

    # Intercept parameter of snow
    itc_snow = n0s*itc_fac

    # ------------Plot-------------------------
    fig = plt.figure( figsize=(12,6), dpi=300 )
    ax = plt.subplot(1,1,1)
    
    # Set color for lines
    # define the colormap you want to sample from
    colormap = plt.get_cmap('plasma_r')  # Replace 'viridis' with the desired colormap
    # number of evenly spaced samples
    num_samples = len(alpha) 
    # create an array of evenly spaced values between 0 and 1
    values = np.linspace(0, 1, num_samples)
    # map the values to colors from the colormap
    colors = [colormap(value) for value in values]
    # black for original value
    colors[2] = (0,0,0,1)

    # Set the y-axis to use a logarithmic scale
    ax.set_yscale('log')

    for ia in range(len(alpha)):
        ax.plot(Tx,itc_snow[ia,:],label="{0:.2f}".format(alpha[ia]),color=colors[ia],linewidth=3)

    ax.grid(True,linestyle='--',alpha=0.5)

    # Set lables
    ax.legend(loc='upper left',fontsize='15')
    ax.set_xlim(xmin=0,xmax=T0-Tcold)
    ax.set_ylim(ymax=1e12)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_xlabel('Supercool: 273.15-T (k)',fontsize=15)
    ax.set_ylabel('N0s (# m-4)',fontsize=15)
    
    ax.set_title( 'Intercept parameter for snow (temperature dependent)',fontweight="bold",fontsize=15 )
    
    plt.savefig( 'N0s.png')
    #plt.grid(True)
    #plt.legend()
    plt.close()


if __name__ == '__main__':

   HDC6() 

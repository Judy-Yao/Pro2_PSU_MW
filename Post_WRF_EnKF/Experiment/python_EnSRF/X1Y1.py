
import os,fnmatch # functions for interacting with the operating system
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
#from matplotlib. import gridspec
from copy import deepcopy

def plot_test():

    # Generate x and y values
    x_vals = np.linspace(0, 60, 300)
    y_vals = perfect_XtoHx(x_vals)

    # Plot
    fig = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    gs = GridSpec(1, 1, figure=fig)

    # Main plots
    ax = plt.subplot(gs[0,0])

    # Truth
    ax.scatter(xt, yt, s=30, color='green', alpha=0.7, edgecolors='green')

    # Observation
    ax.scatter(0,yo,s=30,marker='*',color='red')
    ax.axhline(y=yo,color='red',linestyle='--',label='Observation')

    # Perfect model 
    ax.plot(x_vals, y_vals, label='Perfect Model: Cloud and IR', linewidth=2, zorder=0)

    # Imperfect background ensemble
    ax.scatter(xb_ens, hxb_ens, s=5, color='black', alpha=0.7, edgecolors='black')
    ax.scatter(xb_mean, hxb_mean, s=30, color='black', alpha=0.7, edgecolors='black')

    # Imperfect analysed ensemble
    ax.scatter(d_da['xa_ens'],d_da['hxa_ens'], s=5, color='purple', alpha=0.7, edgecolors='purple')
    ax.scatter(d_da['xa_mean'],d_da['hxa_mean'], s=30, color='purple', alpha=0.7, edgecolors='purple')

    # Create histogram on top
    #ax_bottom = fig.add_axes([0.13, 0.08, 0.62, 0.1])
    #ax_bottom.hist(xb_ens, bins=5, color='gray', alpha=0.7)
    #ax_bottom.axis('off')  # Hide axis for a cleaner look

    # Create left marginal histogram (flattened)
    #ax_left = fig.add_axes([0.05, 0.2, 0.1, 0.62])
    #ax_left.hist(hxb_ens, bins=5, orientation='horizontal', color='gray', alpha=0.7)
    #ax_left.axis('off')  # Hide axis for a cleaner look

    #ax.axhline(y=195, color='gray', linestyle='--', label='Asymptote: y = 195')

    # background
    #ax.scatter(xb_mean, hxb_mean, s=5, color='black', alpha=0.7, edgecolors='black')


    # analysis
    #ax.scatter(xa  )

    ax.set_xlim(0,50)
    ax.set_ylim(190,250)
    ax.set_yticks(range(190, 250, 5))
    #ax.title('Exponential Decay Curve: y = 50e^{-0.1304x} + 195')
    ax.set_xlabel('Cloud ')
    ax.set_ylabel('IR Tbs (K)')
    ax.grid(True)
    ax.legend()
    #plt.tight_layout()


    plt.savefig("test.png", dpi=300, )





# Calculate the dd^
def calculate_var_inno_single( y, mean_xb ):

    return ( y - mean_xb )**2

# Generate N samples
def generate_ensemble():

    # File path for storing the sequence
    file_path = "random_X1Y1/generate_ensemble.npy"

    # Check if the file already exists
    if os.path.exists(file_path):
        # Load the existing sequence
        random_seq = np.load(file_path)
        print("Loaded existing sequence: "+ file_path)
    else:
        # Set a fixed seed
        np.random.seed(5)
        random_seq = np.random.normal(loc=0, scale=stddev_ens_xb, size=Nens)
        np.save(file_path, random_seq)
        print("Generated and saved new sequence:"+ file_path)

    xb_ens = ens_mean_xb + random_seq

    return xb_ens

# Calculate HBH^ 
def calculate_HBH( hx_ens, hx_mean ):

    sum_var = np.sum(( hx_ens - hx_mean ) ** 2)
    return sum_var/(Nens-1)



# Perfect model mapping from cloud amount to IR Tbs
def perfect_XtoHx( x ):

    hx = [ 50*np.exp(-0.1304*ix)+195 for ix in x]
    return hx

def perfect_HxtoX( hx ):

    if isinstance( hx,int ):
        x = -1 / 0.1304 * np.log((hx - 195) / 50)
    else:
        x = [-1 / 0.1304 * np.log((iy - 195) / 50) for iy in hx]
    return x


# Imperfect model mapping from cloud amount to IR Tbs
def imperfect_XtoHx( x ):

    # File path for storing the sequence
    file_path = "random_X1Y1/imperfect_XtoHx.npy"

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
def yt_to_yo( yt ):

    # File path for storing the sequence
    file_path = "random_X1Y1/yt_to_yo.npy"

    # Check if the file already exists
    if os.path.exists(file_path):
        # Load the existing sequence
        random_seq = np.load(file_path)
        print("Loaded existing sequence:"+ file_path)
    else:
        # Set a fixed seed
        np.random.seed(5)
        random_seq = np.random.normal(loc=0, scale=sigma_ot)
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
    x = deepcopy( xb_ens )
    xm = np.mean( x )
    hx = deepcopy( hxb_ens  )

    # Innovation
    y_hxm = yo - hxb_mean

    # Variance in modeled observation space
    var_hx = 0
    hx_pert = []
    for ihxb in hx:
        tmp = ihxb - np.mean(hx) # perturbation
        hx_pert.append( tmp )
        var_hx = var_hx + tmp*tmp

    # Total observation variance
    d = var_hx/(Nens-1) + sigma_ot**2

    # Alpha used in EnSRF
    alpha = 1/(1+np.sqrt(sigma_ot**2/d))

    # Compute Kalman gain matrix
    km = 0
    for i in range(Nens):
        km = km + xb_ens[i]*hx_pert[i] # covariance
    km = ( km/(Nens-1) )/d

    #---------------------------------
    # Update model
    #---------------------------------
    # update model perturbation
    for i in range(Nens):
        x[i] = x[i] - km * alpha * hx_pert[i]

    # update model mean
    xm = xm + km * y_hxm

    #---------------------------------
    # Update simulated observations
    #---------------------------------
    #for inob in range(len(yo)):   # inob: index of next obs
    hx_m = np.mean( hxb_ens )
    print('hxb_ens: ',hxb_ens)
    print('hx_m ', hx_m)

    cov_obs = 0
    for i in range(Nens):
        cov_obs = cov_obs + hx_pert[i]*hx_pert[i]   #( hx_pert[in] -  )

    # update perturbation
    for i in range(Nens):
        hx[i] = hx[i] - alpha * cov_obs / (Nens-1) * hx_pert[i]/d

    #????

    # update mean
    hx_m = hx_m + cov_obs / (Nens-1) * (yo-hx_m)/d  

    print('After: hx', hx)
    print('After: hx_m', hx_m)
    d_da = {'xa_ens':x+xm,'xa_mean':xm,'hxa_ens':hx,'hxa_mean':hx_m}

    return d_da

if __name__ == '__main__':

    #------ Default Setup -------
    # number of ensemble members
    Nens = 10
    # truth
    yt = 200
    xt = perfect_HxtoX( yt )
    # observations
    #xo = [25,]
    #yo = perfect_XtoHx( xo )
    # instrument observation error
    sigma_ot = 1

    # model error (also assume Gaussian)
    sigma_x_inobs = 1.5 # model error in observation space

    # prior ensemble mean in cloud amount
    ens_mean_xb = 4
    # standard deviation of prior ensemble in cloud amount
    stddev_ens_xb = 3

    #------------------------------------

    # Make the observation in observation space
    yo =  yt_to_yo( yt ) 

    # Generate the ensemble in cloud space
    xb_ens = generate_ensemble()
    xb_mean = np.mean( xb_ens )

    # Simulate observations 
    hxb_ens = imperfect_XtoHx( xb_ens )

    # Generate the ensemble mean in IR-Tbs space
    hxb_mean = np.mean( hxb_ens )

    # Perform the EnSRF
    d_da = EnSRF( )
    
    # Generate the covariance between the cloud and IR Tbs
    #cov = covxy( xb_ens, hxb_ens )
    # Generate the variance for ensemble in IR-Tbs space
    #var_hxb = calculate_HBH( hxb_ens, hxb_mean )
    # Generate the updated x
    #xa = ens_mean_xb + (cov/var_hxb)*(yo-hxb_mean)


    plot_test()





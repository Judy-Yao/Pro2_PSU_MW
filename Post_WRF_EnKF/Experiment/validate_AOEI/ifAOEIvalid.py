############################################################################################################
# Symbol definition: 
#    yo: observation, <H(xb)>: ensemble mean in observation space, H(xb)_i: ensemble member,
#    sigma_ot: true observation error, sigma_Hxb: ensemble-estimated standard deviation in observation space.
#    d: innovations/first guess departures of observations: yo - <H(xb)>
#    B: background error covariance, R: true observation error covariance
#    HBH^: background error covariance transformed to observation space with addition of error variance in forward operator
#-----------------------------------------------------------------------------------------------------------
# The AOEI (Minamide and Zhang, 2019) claims that the below relationship holds true:
# - for a single observation:
#    [ yo - <H(xb)> ]**2 = sigma_Hxb**2 + sigma_ot**2   
# - for a vector of observations:
#    dd^ = <HBH^> + R
#
# Note: with an member of 1000 members, <HBH^> is very close to stddev_ens_xb**2. 


# However, notice that in Geer and Bauer (2011) or in Andersson (2004), the relationship is
#    <dd^> = HBH^ + R
# Namely, the covariance of innnovation/FG departure is a sum of all the errors.
# In the ensemble-based DA, the covariance of innovation/FG departure can be estimated as:
# - for a single observation:
#      sum_var = 0
#      for i = 1: Nes:
#          sum_var = sum_var + [ yo - H(xb)_i ]**2
#      <dd^> = sum_var /(Nens - 1)
#
# Notice there is a clear inconsistency between the AOEI assumption and the Geer and Bauer (2011). This script is to test if AOEI is accurate or if it should be revised using a single example for one IR Tb observation (K) .


import os,fnmatch # functions for interacting with the operating system
import numpy as np
import math
import matplotlib.pyplot as plt

def plot_impact_limited_member(  ):

    # ------------ Figure 1
    fig = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    ax = fig.add_subplot(1,1,1)
    # Create a scatter plot
    ax.scatter(var_inno_single, var_inno_overEns, s=5, c='black', alpha=0.7, edgecolors='black', linewidth=0.5)
    
    # Axes
    min_axis = int(min( min(var_inno_single),min(var_inno_overEns) )) - 5
    max_axis = int(max( max(var_inno_single),max(var_inno_overEns) )) + 5
    # Set the same axis limits for both x and y
    ax.set_xlim(min_axis,max_axis)
    ax.set_ylim(min_axis,max_axis)
    # Set custom x and y ticks
    ax.set_xticks(range(min_axis, max_axis, 30))  # Tick marks from 0 to 12 with a step of 2
    ax.set_yticks(range(min_axis, max_axis, 30)) 
    # Add labels and legend
    ax.set_title("Variance of innovation: AOEI v.s. Geer and Bauer 2011",fontsize=15)
    ax.set_xlabel("AOEI: d,"+r"$d^T$",fontsize=15)
    ax.set_ylabel("Geer and Bauer 2011: <d,"+r"$d^T$"+'>',fontsize=15)

    # Add the line y = x
    line_x = [min_axis, max_axis]  # Ensure it matches the x-axis range
    line_y = [min_axis, max_axis]  # Ensure it matches the y-axis range
    ax.plot(line_x, line_y, color='gray', linestyle='-', linewidth=2,)

    # Optional: Add a grid
    ax.grid(True, linestyle='--', alpha=0.5)
    
    plt.savefig("scatter_plot.png", dpi=300, )


    # ------------ Figure 2
    fig = plt.figure( figsize=(6.5,6),dpi=200) # standard: 6.5,8.5
    ax = fig.add_subplot(1,1,1)
    # Create a scatter plot
    diff_Geer = [var_inno_overEns[i]-HBH_and_R[i] for i in range(len(HBH_and_R))]
    diff_aoei = [var_inno_single[i]-HBH_and_R[i] for i in range(len(HBH_and_R))]
    print('mean of diff_Geer:',np.mean( diff_Geer ))
    print('mean of diff_aoei:',np.mean( diff_aoei ))
    # axis range
    min_y_axis = int(min( min(diff_Geer),min(diff_aoei) )) - 5
    max_y_axis = int(max( max(diff_Geer),max(diff_aoei) )) + 5

    ax.scatter(HBH_and_R, diff_Geer, s=5, c='blue', alpha=0.7, edgecolors='black', linewidth=0.5)
    ax.set_xlabel("HB"+r"$H^T$"+"+R", fontsize=12)
    ax.set_ylabel("<d,"+r"$d^T$"+'> - '+"(HB"+r"$H^T$"+"+R)",fontsize=15)
    # Set the same axis limits for both x and y
    ax.set_ylim(min_y_axis,max_y_axis)
    ax.set_yticks(range(min_y_axis, max_y_axis, 20))

    # Create a secondary y-axis for C
    ax1 = ax.twinx()
    ax1.scatter(HBH_and_R, diff_aoei, s=5, color='red', )
    ax1.set_ylabel("d,"+r"$d^T$"+' - '+"(HB"+r"$H^T$"+"+R)", color='red', fontsize=15)
    ax1.set_ylim(min_y_axis,max_y_axis)
    ax1.set_yticks(range(min_y_axis, max_y_axis, 20))


    # Set the same axis limits for both x and y
    # Set custom x and y ticks
    #ax.set_xticks(range(min_axis, max_axis, 30))  # Tick marks from 0 to 12 with a step of 2
    #ax.set_yticks(range(min_axis, max_axis, 30))
    # Add labels and legend
    #ax.set_title("Variance of innovation: AOEI v.s. Geer and Bauer 2011",fontsize=15)
    #ax.set_xlabel("AOEI: d,"+r"$d^T$",fontsize=15)
    #ax.set_ylabel("Geer and Bauer 2011: <d,"+r"$d^T$"+'>',fontsize=15)

    # Add the line y = x
    #line_x = [min_axis, max_axis]  # Ensure it matches the x-axis range
    #line_y = [min_axis, max_axis]  # Ensure it matches the y-axis range
    #ax.plot(line_x, line_y, color='gray', linestyle='-', linewidth=2,)

    # Optional: Add a grid
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.savefig("scatter_plot1.png", dpi=300, )











# Calculate the dd^
def calculate_var_inno_single( y, mean_xb ):

    return ( y - mean_xb )**2

# Generate N samples
def generate_ensemble():
   
    return np.random.normal(loc=ens_mean_xb, scale=stddev_ens_xb, size=Nens)

# Calculate <dd^> 
def calculate_var_inno_overEns( y,ens ):

    sum_var = np.sum((y - ens) ** 2)
    return sum_var/(Nens-1)
    
# Calculate HBH^ 
def calculate_HBH( ens, y ):

    sum_var = np.sum((ens - y ) ** 2)
    return sum_var/(Nens-1)




if __name__ == '__main__':

    #------ Default Setup -------
    # number of ensemble members
    Nens = 60
    # truth
    yt = 220
    # observation
    yo = 219
    # prior ensemble mean in radiance space
    ens_mean_xb = 220
    # standard deviation of prior ensemble in radiance space
    stddev_ens_xb = 15
    # true observation error
    sigma_ot = 1
    #------------------------------------
    impact_limited_member = True
    
    # Same limited ensemble size; Repeat the experiment for numerous times
    if impact_limited_member:

        Nens = 1000
        Realization = 500
        yt = 220
        ens_mean_xb = 220
        stddev_ens_xb = 15
        sigma_ot = 0

        HBH_and_R = []
        var_inno_single = []
        var_inno_overEns = []
        for ir in range(Realization):
            # Generate the ensemble for the ir-th experiment
            ens = generate_ensemble()
            # Calculate the ensemble mean
            sample_mean = np.mean( ens )
            # HBH^: background error variance in observation space
            HBH_and_R.append( calculate_HBH( ens, sample_mean) + sigma_ot**2  )
            # Calculate the dd^ (innovation squared)
            var_inno_single.append( calculate_var_inno_single( yt,sample_mean ) )
            # Calculate the covariance of innovation squared evaluated over the ensemble
            var_inno_overEns.append( calculate_var_inno_overEns( yt,ens ) )
        # Plot
        plot_impact_limited_member()

    # HBH^
    #HBH = calculate_var_inno_overEns( yt )

    # Calculate the dd^ (innovation squared)
    #var_inno_single = var_inno_single()
    
    # Calcuate the variance of innovation
    #var_inno_overEns = calculate_var_inno_overEns( yo )

    #print('HBH^ + R is: ')
    #print( HBH + sigma_ot**2 )
    #print(stddev_ens_xb**2+sigma_ot**2)

    #print('dd^ is: ')
    #print( var_inno_single )

    #print('<dd^> is')
    #print(var_inno_overEns)












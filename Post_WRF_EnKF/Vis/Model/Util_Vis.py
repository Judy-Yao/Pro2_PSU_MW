#!/usr/bin/env python3

# This module consists of functions for visualizing post-EnKF products.
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap

# --------------------------------------------------------------
# ---------------------- COLOR MAP -----------------------------
# --------------------------------------------------------------
def HydroIncre( N_nega, N_posi ):
   
    print( N_nega)
    print( N_posi)
    N = N_nega+N_posi
    buff = 1                 
    # Define default colors         
    end_nega_color = np.array([0, 0, 1])
    end_posi_color = np.array([1,0,0])      
    white_color = np.array([1,1,1])
    black_color = np.array([0,0,0])             
                                                    
    # Initialize map
    Hydromap = np.zeros([N,3])                          
                                                            
    # place red, white, blue                                                
    for i in range( N ):
        if i <= ( N_nega-buff ): 
            near_nega = i
            near_center = (N_nega-buff)-i                                                                  
            tmp = (near_center/(N_nega-buff))*end_nega_color + (near_nega/(N_nega-buff))*white_color
            Hydromap[i,:] = tmp
        elif i >= ( N_nega ):
            near_posi = (N_nega+N_posi-buff-i) 
            near_center = i-( N_nega )
            tmp = (near_center/(N_posi-buff))*end_posi_color + (near_posi/(N_posi-buff))*white_color 
            Hydromap[i,:] = tmp 
        else:
            Hydromap[i,:] = white_color                                                                     
    arr = np.append(Hydromap, np.ones( [N,1] ), axis=1)     
    HydroMap = ListedColormap( arr )                                                                        
    return HydroMap 


# Plot full Tbs
def newJet(max_T, min_T, min_Jet):

    max_T = 300
    min_T = 80
    min_Jet = 150
    jetLength = max_T - min_Jet + 1
    notJetLength = (min_Jet - 1) - (min_T + 1) + 1 

    jet = cm.get_cmap('jet', max_T-min_T)
    Myjet = jet(np.linspace(0,1,jetLength))

    jet_red = Myjet[:,0]
    jet_green = Myjet[:,1]
    jet_blue = Myjet[:,2]

    jetAdd_red = np.linspace(1.0, jet_red[0], notJetLength)
    jetAdd_green = np.linspace(1.0, jet_green[0], notJetLength)
    jetAdd_blue  = np.linspace(1.0, jet_blue[0], notJetLength) 

    cm_red = np.concatenate([np.insert(jetAdd_red,0,0),np.append(jet_red,0)])
    cm_green = np.concatenate([np.insert(jetAdd_green,0,0),np.append(jet_green,0)])
    cm_blue = np.concatenate([np.insert(jetAdd_blue,0,0),np.append(jet_blue,0)])

    newJet_value = np.column_stack([cm_red,cm_green,cm_blue, np.ones(len(cm_red))])
    newJet = ListedColormap(newJet_value)
    
    return newJet

# Plot Tb differences
def newRWB(max_T, min_T, min_RWB):

    max_T = 10
    min_T = 10
    min_RWB = 8
    RWBLength = 18 # max_T - min_RWB + 1
    notRWBLength = 2 #(min_RWB - 1) - (min_T + 1) + 1

    RWB = cm.get_cmap('RdBu_r', max_T-min_T)
    MyRWB = RWB(np.linspace(0,1,RWBLength))

    RWB_red = MyRWB[:,0]
    RWB_green = MyRWB[:,1]
    RWB_blue = MyRWB[:,2]

    #RWBAdd_red = np.linspace(1.0, RWB_red[0], notRWBLength)
    #RWBAdd_green = np.linspace(1.0, RWB_green[0], notRWBLength)
    #RWBAdd_blue  = np.linspace(1.0, RWB_blue[0], notRWBLength)

    cm_red = np.append( np.insert(RWB_red,0,0), 0)
    cm_green = np.append( np.insert(RWB_green,0,0), 0)
    cm_blue = np.append( np.insert(RWB_blue,0,0), 0)

    newRWB_value = np.column_stack([cm_red,cm_green,cm_blue])
    #newRWB_value = np.column_stack([cm_red,cm_green,cm_blue, np.ones(len(cm_red))])
    newRWB = ListedColormap(newRWB_value)

    return newRWB



# colormap for plotting infrared BT
def IRcmap(c_int):
    
    cmap = np.zeros( (int(140/c_int),3) ) 

    cmap[0:int(5/c_int), 0] = np.append( np.arange(0.5, 1, c_int/(5-c_int)*0.5), 1.0)
    cmap[0:int(5/c_int), 1] = 0
    cmap[0:int(5/c_int), 2] = np.append( np.arange(0.5, 1, c_int/(5-c_int)*0.5), 1.0)
    

    cmap[int(5/c_int)-1:int(10/c_int):int(c_int/c_int),0] = 1
    cmap[int(5/c_int)-1:int(10/c_int):int(c_int/c_int),1] = np.append( np.arange(0, 0.5, c_int/5*0.5), 0.5)
    cmap[int(5/c_int)-1:int(10/c_int):int(c_int/c_int),2] = 1

    cmap[int(10/c_int)-1:int(20/c_int):int(c_int/c_int),0] = np.append( np.arange(0.8, 0.0, -c_int/10*0.8), 0.0)
    cmap[int(10/c_int)-1:int(20/c_int):int(c_int/c_int),1] = np.append( np.arange(0.8, 0.0, -c_int/10*0.8), 0.0)
    cmap[int(10/c_int)-1:int(20/c_int):int(c_int/c_int),2] = np.append( np.arange(0.8, 0.0, -c_int/10*0.8), 0.0)

    cmap[int(20/c_int)-1:int(30/c_int):int(c_int/c_int),0] = np.append( np.arange(0, 1.0, c_int/10), 1)
    cmap[int(20/c_int)-1:int(30/c_int):int(c_int/c_int),1] = 0
    cmap[int(20/c_int)-1:int(30/c_int):int(c_int/c_int),2] = 0

    cmap[int(30/c_int)-1:int(40/c_int):int(c_int/c_int),0] = 1
    cmap[int(30/c_int)-1:int(40/c_int):int(c_int/c_int),1] = np.append( np.arange(0, 1.0, c_int/10), 1)
    cmap[int(30/c_int)-1:int(40/c_int):int(c_int/c_int),2] = 0

    cmap[int(40/c_int)-1:int(50/c_int):int(c_int/c_int),0] = np.append( np.arange(1.0, 0, -c_int/10), 0) 
    cmap[int(40/c_int)-1:int(50/c_int):int(c_int/c_int),1] = 1
    cmap[int(40/c_int)-1:int(50/c_int):int(c_int/c_int),2] = 0

    cmap[int(50/c_int)-1:int(60/c_int):int(c_int/c_int),0] = 0
    cmap[int(50/c_int)-1:int(60/c_int):int(c_int/c_int),1] = np.append( np.arange(1.0, 0, -c_int/10), 0)  
    cmap[int(50/c_int)-1:int(60/c_int):int(c_int/c_int),2] = np.append( np.arange(0, 1.0, c_int/10), 1)
    #cmap[(50:c_int:60)/c_int,0] = 0
    #cmap[(50:c_int:60)/c_int,1] = 1:-c_int/10:0
    #cmap[(50:c_int:60)/c_int,2] = 0:c_int/10:1

    cmap[int(60/c_int)-1:int(70/c_int):int(c_int/c_int),0] = 0
    cmap[int(60/c_int)-1:int(70/c_int):int(c_int/c_int),1] = np.append( np.arange(0, 1.0, c_int/10), 1)
    cmap[int(60/c_int)-1:int(70/c_int):int(c_int/c_int),2] = 1
    #cmap[(60:c_int:70)/c_int,0] = 0
    #cmap[(60:c_int:70)/c_int,1] = 0:c_int/10:1
    #cmap[(60:c_int:70)/c_int,2] = 1
    
    cmap[int(70/c_int)-1:int(140/c_int):int(c_int/c_int),0] = np.append( np.arange(1.0, 0, -c_int/70), 0)
    cmap[int(70/c_int)-1:int(140/c_int):int(c_int/c_int),1] = np.append( np.arange(1.0, 0, -c_int/70), 0)
    cmap[int(70/c_int)-1:int(140/c_int):int(c_int/c_int),2] = np.append( np.arange(1.0, 0, -c_int/70), 0)
    #cmap[(70:c_int:140)/c_int,0] = 1:-c_int/70:0
    #cmap[(70:c_int:140)/c_int,1] = 1:-c_int/70:0
    #cmap[(70:c_int:140)/c_int,2] = 1:-c_int/70:0

    IRcmap = ListedColormap( cmap )
    
    return IRcmap

if __name__ == '__main__':
    
    #RWB = cm.get_cmap('RdBu_r', 256)    
    #MyRWB = RWB(np.linspace(0,1,255))
    #print(MyRWB)
    np.random.seed(19680801)
    data = np.random.randn(30, 30)
    cmap = HydroIncre(7,6)
    fig, ax = plt.subplots(1, 1, figsize=(6, 3), constrained_layout=True)
    psm = ax.pcolormesh(data, cmap=cmap, rasterized=True)
    fig.colorbar(psm, ax=ax)
    plt.show()


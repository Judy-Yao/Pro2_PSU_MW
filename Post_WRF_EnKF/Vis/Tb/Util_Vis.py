#!/usr/bin/env python3

# This module consists of functions for visualizing post-EnKF products.
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap

# --------------------------------------------------------------
# ---------------------- COLOR MAP -----------------------------
# --------------------------------------------------------------
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

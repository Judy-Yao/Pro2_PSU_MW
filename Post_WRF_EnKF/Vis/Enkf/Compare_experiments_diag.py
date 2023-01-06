#!/work2/06191/tg854905/stampede2/opt/anaconda3/lib/python3.7

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import pickle
import netCDF4 as nc
import Diagnostics as Diag
import math
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time






def Diagnose(big_dir, Storm, Exper, time, v_interest):
   
    # Read in diagnostics for each experiment
    Diag_Expers = {}
    for iExper in Exper:
        filename = big_dir+Storm+'/'+iExper+'/run/'+time+'/enkf/d03/fort.10000'
        Diag_Expers[iExper] = Diag.read_diag( filename, v_interest )

    # 




if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'

    # configuration
    Storm = 'HARVEY'
    Exper_name = ['JerryRun/MW_THO','J_DA+J_WRF+J_init']
    time = '201708221200'
    #filename = '/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/IR+MW-J_DA+J_WRF+J_init-SP-intel19/run/201708221200/enkf/d03/fort.10000'
    v_interest = ['obs_type','lat','lon','obs','prior_mean','posterior_mean']
    d_diag = read_diag( filename, v_interest )

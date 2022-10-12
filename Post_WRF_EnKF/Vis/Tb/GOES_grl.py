#!/usr/bin/env python                                                                                                                                      
# This script plots BTs observed by GOES-16 during the times chosen in the settings below.  It also plots simulated BTs as calculated by the CRTM.
# It works with the BTHV.py script in the DIR: 4 runs.

import functions
import netCDF4 as nc
import numpy as np
import datetime
import os
from math import pi
import sys 
import scipy as sp
import scipy.ndimage
from scipy import interpolate
from scipy.interpolate import griddata

def G16(num):    
        from datetime import datetime
        time = datetime.strptime(str(num), '%Y%m%d%H%M')
        storm_name = 'HARVEY'
        dom=3
        ch_list = [8] # channel from observations
        ch_fig = 8 #the channel you want to plot
        ch_list2 = [8,9,10,14] # list of channels you calculated with CRTM
        xmax=5424  # Obs max x value
        ymax=5424  # Obs max y value
        xmax2=297  # CRTM maximum x value
        ymax2=297  # CRTM maximum y value
        n_ch = len(ch_list)

        ### This subsection reads in the observations ###
        G16OBS_DIR = '/work/06191/tg854905/stampede2/data/GOES_16/G16_obs/'

        import datetime
        day = time - datetime.datetime(time.year-1,12,31,0,0) # This finds the day of the year. Julian Date
        Tb_obs = np.zeros([n_ch,ymax,xmax])
        for ch in ch_list:
            filename_part = 'OR_ABI-L2-CMIPF-M3C'+"{0:02}".format(ch)+'_G16_s'+time.strftime('%Y')+str(day.days)+time.strftime('%H%M')
            files = os.listdir(G16OBS_DIR)
            for filenc in files:
                if filename_part in filenc:
                    filename = G16OBS_DIR+'/'+filenc
                    filecheck = True
            Tb_tmp = functions.read_GOES16(filename)
            Tb_obs[ch_list.index(ch),:,:] = Tb_tmp['tb']
            del filename
        lons_obs = np.array(Tb_tmp['lons'], dtype=float)   # This forces lons_obs to be an array
        lats_obs = np.array(Tb_tmp['lats'], dtype=float)
        lons_obs[np.isnan(lons_obs)] = 1e9
        lats_obs[np.isnan(lats_obs)] = 1e9

         ### This subsection reads in the simulated BTs ###
        MW_a = '/work/06191/tg854905/stampede2/data/GOES_16/simu_new/MW_lmt19/Radiance_mean_output_d03_2017-08-24_00:00.bin'
        MW_f = '/work/06191/tg854905/stampede2/data/GOES_16/simu_new/Forecast/Radiance_MW_lmt19_mean_input_d03_2017-08-24_00:00.bin'
        IR_a = '/work/06191/tg854905/stampede2/data/GOES_16/simu_new/IR_CONV/Radiance_mean_output_d03_2017-08-24_00:00.bin' 
        IR_f = '/work/06191/tg854905/stampede2/data/GOES_16/simu_new/Forecast/Radiance_IR_mean_input_d03_2017-08-24_00:00.bin'

 
        data_MW_fore = functions.read_crtm(MW_f,xmax2,ymax2,ch_list2)
        MW_fore = data_MW_fore[ch_fig]
        
        data_MW_ana = functions.read_crtm(MW_a,xmax2,ymax2,ch_list2)
        MW_ana = data_MW_ana[ch_fig]

        data_IR_fore = functions.read_crtm(IR_f,xmax2,ymax2,ch_list2)
        IR_fore = data_IR_fore[ch_fig]

        data_IR_ana = functions.read_crtm(IR_a,xmax2,ymax2,ch_list2)
        IR_ana = data_IR_ana[ch_fig]

        lons_model = np.array(data_IR_ana['lons'], dtype=float)     # This forces lons_model to be an array
        lats_model = np.array(data_IR_ana['lats'], dtype=float)

        x_min = np.min(lons_model[0,:])
        x_max = np.max(lons_model[0,:])
        y_min = np.min(lats_model[:,0])
        y_max = np.max(lats_model[:,0])

        lons_obs_new = lons_obs[(lons_obs >= x_min) & (lons_obs <= x_max) & (lats_obs >= y_min) & (lats_obs <= y_max)]
        lats_obs_new = lats_obs[(lats_obs >= y_min) & (lats_obs <= y_max) & (lons_obs >= x_min) & (lons_obs <= x_max)]
        Tb_obs_new = Tb_obs[0, (lons_obs >= x_min) & (lons_obs <= x_max) & (lats_obs >= y_min) & (lats_obs <= y_max)]
        lons_obs = lons_obs_new
        lats_obs = lats_obs_new
        Tb_obs = Tb_obs_new
        #print('max_model_lon:', np.max(lons_model), 'min_model_lon:', np.min(lons_model))

        #plt.savefig('/work/06191/tg854905/stampede2/mwir/test.png', dpi=300)
        my_list = [lons_obs, lats_obs, Tb_obs, lons_model, lats_model, MW_fore, MW_ana, IR_fore, IR_ana]
        return(my_list)







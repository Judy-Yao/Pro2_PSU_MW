#!/usr/bin/env python                                                                                                                                      

from netCDF4 import Dataset
from wrf import getvar,interplevel
import numpy as np
import scipy as sp
import scipy.ndimage
import datetime
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange
from math import pi
from matplotlib.colors import LinearSegmentedColormap as LSC
import sys

Req = 6378137.
Rpol = 6356752.31414
H = 42164160.
ramda_o = -89.5*pi/180.
def read_GOES16(filename, lonlat=True):
    Tb_dict = {}
    print(filename)
    ncfile = Dataset(filename)
    # cloud and moisture imagery TB
    tb = ncfile.variables['CMI'][:]
    Tb_dict['Yo'] = tb[::-1,:]
    if lonlat:
        # GOES fixed grid projection x-coordinate (units: rad)
        x = ncfile.variables['x'][:]
        # GOES fixed grid projection y-coordinate (units: rad) 
        y = ncfile.variables['y'][:]
        xx, yy = np.meshgrid(x,y)
        a = np.sin(xx)**2 + np.cos(xx)**2 * (np.cos(yy)**2 + Req**2/Rpol**2 * np.sin(yy)**2)
        b = - 2.* H * np.cos(xx) * np.cos(yy)
        c = H**2 - Req**2
        Rs = (-b - np.sqrt(b**2 - 4.*a*c)) / (2.*a)
        Sx = Rs * np.cos(xx) * np.cos(yy)
        Sy = -Rs * np.sin(xx)
        Sz = Rs * np.cos(xx) * np.sin(yy)
        lats = np.arctan( Req**2/Rpol**2 * Sz / np.sqrt((H - Sx)**2 + Sy**2))
        lons = ramda_o - np.arctan(Sy / (H - Sx))
        Tb_dict['lon'] = lons[::-1,:]/pi*180.
        Tb_dict['lat'] = lats[::-1,:]/pi*180.

    ncfile.close()
    return Tb_dict

def read_crtm(filename,xmax,ymax,ch_list):
    print(filename)
    data = np.fromfile(filename,dtype='>f4')
    n_ch = len(data)/(xmax*ymax)-2
    n_ch = int(n_ch)
    #print('n_ch:',type(n_ch))
    #print('xmax', type(xmax))
    if n_ch != len(ch_list):
        print('Error!! # of channels in data is '+str(n_ch))
    sim = data[:].reshape(n_ch+2,ymax,xmax)
    Tb_dict = {}
    Tb_dict['lons'] = sim[0,:,:]
    Tb_dict['lats'] = sim[1,:,:]
    for rec in range(n_ch):
        Tb_dict[ch_list[rec]] = sim[rec+2,:,:]
    return Tb_dict

def ir_colormap(vmax=325): 
        ## for 173.15 K ~ 333.15 K
        vmin = 185  
        clevs = 0.5       
        interval = np.arange(vmin,vmax+clevs,clevs)
        clevs = 20        
        ticks = np.arange(180,int(vmax)+clevs,clevs)
        theresholds = [15., 20., 30., 40., 50., 55., 60., 70.]
        cdict = {'red':  ((0.0,    0.9297, 0.9297), # violet
                          (theresholds[0]/(float(vmax)-vmin), 0.1   , 0.1   ), # black red
                          (theresholds[1]/(float(vmax)-vmin), 0.4844, 0.4844), # Barn red
                          (theresholds[2]/(float(vmax)-vmin), 1.0   , 1.0   ), # red
                          (theresholds[3]/(float(vmax)-vmin), 1.0   , 1.0   ), # yellow
                          (theresholds[4]/(float(vmax)-vmin), 0.0   , 0.0   ), # green
                          (theresholds[5]/(float(vmax)-vmin), 0.0977, 0.0977), # midnight blue
                          (theresholds[6]/(float(vmax)-vmin), 0.0   , 0.0   ), # blue
                          (theresholds[7]/(float(vmax)-vmin), 1.0   , 1.0   ), # white
                          (1.0,    0.0   , 0.0   )),# black
                'green':((0.0000, 0.5078, 0.5078),
                         (theresholds[0]/(float(vmax)-vmin),0.0   , 0.0   ),
                         (theresholds[1]/(float(vmax)-vmin),0.0391, 0.0391),
                         (theresholds[2]/(float(vmax)-vmin),0.0   , 0.0   ),
                         (theresholds[3]/(float(vmax)-vmin),1.0   , 1.0   ),
                         (theresholds[4]/(float(vmax)-vmin),1.0   , 1.0   ),
                         (theresholds[5]/(float(vmax)-vmin),0.0977, 0.0977),
                         (theresholds[6]/(float(vmax)-vmin),0.0   , 0.0   ),
                         (theresholds[7]/(float(vmax)-vmin),1.0   , 1.0   ),
                         (1.0,    0.0   , 0.0   )),
                'blue': ((0.0000, 0.9297, 0.9297),
                         (theresholds[0]/(float(vmax)-vmin),0.0   , 0.0   ),
                         (theresholds[1]/(float(vmax)-vmin),0.0078, 0.0078),
                         (theresholds[2]/(float(vmax)-vmin),0.0   , 0.0   ),
                         (theresholds[3]/(float(vmax)-vmin),0.0   , 0.0   ),
                         (theresholds[4]/(float(vmax)-vmin),0.0   , 0.0   ),
                         (theresholds[5]/(float(vmax)-vmin),0.4375, 0.4375),
                         (theresholds[6]/(float(vmax)-vmin),1.0   , 1.0   ),
                         (theresholds[7]/(float(vmax)-vmin),1.0   , 1.0   ),
                         (1.0,    0.0   , 0.0   ))
                }
        ircmap = LSC('ircmap',cdict)
        ircmap.set_under('white')
        return ircmap, vmin, vmax, interval, ticks












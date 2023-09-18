#!/usr/bin/env python3

import os # functions for interacting with the operating system
import numpy as np
import glob
from netCDF4 import Dataset
from wrf import getvar
from datetime import datetime, timedelta
import pickle # To save any object in python to the disk (transportation)
import scipy as sp
import scipy.ndimage

import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap

matplotlib.rcParams['lines.linewidth'] = 0.5
matplotlib.rcParams['lines.markersize'] = 2
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['font.size'] = 6

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True

HARVEY_REPORT = {
        'time': [datetime(2017,8,23,12,00,00) + timedelta(hours=t) for t in [0,6,12,18,24,30,36,42,48,54,60,63,66,72,78,84] ],
        'max_ws': np.array([30, 35, 40, 50, 60, 70, 80, 90, 95, 105, 115, 115,  105, 65, 50, 45]) * 0.51444,
        'slp': np.array([1006, 1005, 1003, 997, 986, 978, 973, 966, 949, 943, 941,937, 948, 978, 991, 995]),
        'lat': np.array([21.4,21.6,22.0,22.8,23.7,24.4,25.0,25.6,26.3,27.1,27.8,28.0,28.2,28.7,29.0,29.2]),
        'lon': np.array([92.3,92.4,92.5,92.6,93.1,93.6,94.4,95.1,95.8,96.3,96.8,96.9,97.1,97.3,97.5,97.4]) * (-1),
        }

def read_wrfout(directory, filename_pickle=None, force_reload=False):
    # get all wrfout_d03_* files
    if filename_pickle is None:
        filename_pickle = directory + '/HPI.pickle'
    if not os.path.exists( filename_pickle ) or force_reload:
        filenames = sorted(glob.glob(directory + '/wrfout_d03_*') )
        nstep = len(filenames)
        HPI = {}
        HPI['time'] = [datetime(2017,1,1) for i in range(nstep)]
        HPI['max_ws'] = np.zeros( nstep )
        HPI['slp'] = np.zeros( nstep )
        HPI['lat'] = np.zeros( nstep )
        HPI['lon'] = np.zeros( nstep )
        for ifile in range(nstep):
            filename = filenames[ifile]
            print(filename)
            with Dataset(filename) as ncid:
                start_time = datetime.strptime( ncid.SIMULATION_START_DATE, '%Y-%m-%d_%H:%M:%S')
                dtime = timedelta( minutes= float(ncid.variables['XTIME'][:][0]) )
                HPI['time'][ifile] = start_time + dtime
                ws = np.sqrt( ncid.variables['U10'][:]**2 + ncid.variables['V10'][:]**2 )
                HPI['max_ws'][ifile] = np.max( ws )
                slp = getvar(ncid, 'slp')
                HPI['slp'][ifile] = np.min( slp )
                
                # smooth to find min slp location
                slp_smooth = sp.ndimage.filters.gaussian_filter(slp, [11, 11] )
                idx = np.nanargmin( slp_smooth )
                HPI['lat'][ifile] = ncid.variables['XLAT'][:].flatten()[idx]
                HPI['lon'][ifile] = ncid.variables['XLONG'][:].flatten()[idx]

        # write to pickle
        with open(filename_pickle, 'wb') as f:
            pickle.dump(HPI,f)
    else:
        with open(filename_pickle, 'rb') as f:
            HPI = pickle.load(f)
    return HPI

def plot_one( ax1, ax2, ax3, map, fcst, linestyle, label, step=1):
    dates = matplotlib.dates.date2num( fcst['time'] )
    hours = np.array([date.hour for date in fcst['time']])
    idx = np.mod(hours.astype(int), step) == 0
    ax1.plot_date(dates[idx], fcst['slp'][idx], linestyle, label=label)
    ax2.plot_date(dates[idx], fcst['max_ws'][idx], linestyle, label=label)
#    ax3.plot( fcst['lon'], fcst['lat'], linestyle, label=label)
    x, y = map( fcst['lon'][idx], fcst['lat'][idx] )
    map.plot(x,y,linestyle)


def test():
    directory = '/scratch/05012/tg843115/EnKF_crtm_new_run/harvey_mwir5/run/201708240000/wrf_vf_201708241200_201708270000_scott_namelist'
    #directory = '/scratch/05012/tg843115/EnKF_crtm_new_run/harvey_mwir5/run/201708240000/wrf_vf_201708241200_201708270000'
    fcst = read_wrfout(directory)
    dates = matplotlib.dates.date2num( fcst['time'] )
    plt.subplot(1,2,1)
    plt.plot_date(dates, fcst['slp'])
    plt.subplot(1,2,2)
    plt.plot_date(dates, fcst['max_ws'])
    plt.show()

def plot_hpi(output_dir=None):
    reload_data = False
    step = 1
    
    f, ax=plt.subplots(1, 3, subplot_kw=dict(adjustable='box'), sharex='none', sharey='none', figsize=(12,4), dpi=150)
    
    map = Basemap(projection='merc',llcrnrlat=20,urcrnrlat=31,llcrnrlon=-100.5, urcrnrlon=-89.5, resolution='l', ax=ax[2])
    map.drawcoastlines(linewidth=0.8)
    map.fillcontinents(color='coral',lake_color='aqua')
    map.drawparallels( np.arange(20,33,2), labels=[True, False, False, False], color='w', linewidth=0.5)
    map.drawmeridians( np.arange(-100,-88,2), labels=[False, False, False, True], color='w', linewidth=0.5)

    
    plot_one( ax[0], ax[1], ax[2], map, read_wrfout('/scratch/05012/tg843115/EnKF_crtm_new_run/harvey_mwir5/output/201708240000', force_reload=reload_data), 'b-', 'MW/IR new namelist', step=step)
#    plot_one( ax[0], ax[1], ax[2], map, read_wrfout('/scratch/05012/tg843115/EnKF_crtm_new_run/harvey_mwir5/run/201708240000/wrf_vf_201708241200_201708270000_scott_namelist', force_reload=reload_data), 'b--', 'MW/IR old namelist', step=step)
    plot_one( ax[0], ax[1], ax[2], map, read_wrfout('/scratch/05012/tg843115/EnKF_forecast/Harvey_conv_GPM-03x36x24-13x36x24_limitedIR_2312_abei_Green/output/201708240000', force_reload=reload_data), 'g-', 'MW Scott', step=step)
#    plot_one( ax[0], ax[1], ax[2], map, read_wrfout('/scratch/05012/tg843115/EnKF_crtm_upgrade1_run/test_scott_ranch/201708240000/SSMIS_12x36x24/run', force_reload=reload_data), 'y--', 'MW ranch backup', step=step)
    plot_one( ax[0], ax[1], ax[2], map, read_wrfout('/scratch/05012/tg843115/EnKF_forecast/201708231200_GOES16_conv/output/201708240000', force_reload=reload_data), 'r-', 'GOES+conv', step=step)
    plot_one( ax[0], ax[1], ax[2], map, read_wrfout('/scratch/05012/tg843115/EnKF_crtm_new_run/harvey_mw_allhf/output/201708240000', force_reload=reload_data), 'c-', 'MW HF all', step=step)
    plot_one( ax[0], ax[1], ax[2], map, HARVEY_REPORT, 'k-', 'Best track', step=3)
    

    # horizontal lines
    ws_line = [18, 33, 42, 49, 58]
    for ws in ws_line:
        ax[1].axhline(ws)



    ax[0].set_xlim([datetime(2017, 8, 23, 12, 0, 0), datetime(2017, 8, 27)])
    ax[1].set_xlim([datetime(2017, 8, 23, 12, 0, 0), datetime(2017, 8, 27)])
    ax[0].tick_params(axis='x', labelrotation=45)
    ax[1].tick_params(axis='x', labelrotation=45)
    ax[0].legend(frameon=False)
    ax[0].set_title( 'minimum SLP' )
    ax[1].set_title( 'maximum wind speed' )
    ax[2].set_title( 'track' )
    
    if output_dir is not None:                             
        plt.savefig('%s/hpi_fcst_%s_%s.png' % (output_dir, '201708240000', '201708270000') )                                                
    else:
        plt.show()
    plt.clf()


if __name__ == '__main__':
    plot_hpi( output_dir = './out')

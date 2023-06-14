#!/usr/bin/env python3

import os # functions for interacting with the operating system
import numpy as np
from datetime import datetime, timedelta
import glob
import pickle
from netCDF4 import Dataset
import math
import matplotlib
matplotlib.use("agg")
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from Track_intensity import read_bestrack


USAF_idx_var = {'GMT':0,'GPSA':8,'LAT':12,'LON':13,'WSpd':30,'TD':31}


class Info_USAF_clt:
    '''This object is designed to collect information for each wrf day'''
    def __init__(self, filedir, beg_time, end_time):
        self.filedir = filedir
        self.beg_time = beg_time
        self.end_time = end_time

# Check the beginning and ending time for a mission
def Time_range_USAF_mission( mission_path ):
    
    Day_start = os.path.split( mission_path )[1][0:8] # e.g., 20170823
    Day_start_dt = datetime.strptime( Day_start, '%Y%m%d' )
    Day_end_dt = Day_start_dt

    with open(mission_path, 'r') as f:
        all_lines = f.readlines()
    # hhmmss_start and hhmmss_end
    hhmmss_start = all_lines[6].split()[0]
    hhmmss_end = all_lines[-1].split()[0]

    # Figure out the Day_end
    for line in all_lines[6:]:
        split_line = line.split()
        hhmmss = split_line[0] # GMT  
        if hhmmss == '00:00:00':
            Day_end_dt = Day_end_dt + timedelta( days=1 )

    # Update time_start and time_end
    time_start = Day_start_dt + timedelta(hours=int(hhmmss_start[0:2]),minutes=int(hhmmss_start[3:5]),seconds=int(hhmmss_start[6:8]))
    time_end = Day_end_dt + timedelta(hours=int(hhmmss_end[0:2]),minutes=int(hhmmss_end[3:5]),seconds=int(hhmmss_end[6:8]))

    return time_start, time_end

def read_USAF_mission( mission_path, attr_itt ):

    with open(mission_path, 'r') as f:
        all_lines = f.readlines()

    # find the indices of the attributes of interest
    idx_itt = [] 
    for it in attr_itt:
        idx_itt.append( USAF_idx_var[it] )
    
    # create a list of several empty sub lists
    lists = [ [] for i in range(len(attr_itt)) ]    

    # put variables of interest to the sub lists of the large list
    time_mission = os.path.split( mission_path )[1][0:8] # e.g., 20170823
    time_mission_dt = datetime.strptime(time_mission,"%Y%m%d")
    dayAddone = False
    for line in all_lines[6:]:
        split_line = line.split()
        # combine day and hour together
        hhmmss = split_line[idx_itt[0]]
        if hhmmss == '00:00:00':
            dayAddone = True
        if dayAddone == False:
            time_dt = time_mission_dt + timedelta(hours=int(hhmmss[0:2]),minutes=int(hhmmss[3:5]),seconds=int(hhmmss[6:8]))
        if dayAddone == True:
            time_dt = time_mission_dt + timedelta(days=1,hours=int(hhmmss[0:2]),minutes=int(hhmmss[3:5]),seconds=int(hhmmss[6:8]))
        
        lists[0].append( time_dt )
        for it in range( 1,len(attr_itt) ):
            #print(split_line[31])
            lists[it].append( float(split_line[idx_itt[it]])  )
    
    # gather these values into a dictionary
    dict_AF_mission = { }
    for it in range( len(attr_itt) ):
        dict_AF_mission[ attr_itt[it] ] = lists[it]
    
    # sanity check to make sure the correct variables are read
    print( 'Sanity check values of a record...' )
    for ik,iv in dict_AF_mission.items():
        print( ik, iv[15000] )
    
    # identify the unique hour of days
    #all_times_str = [ datetime.strftime(it, "%Y%m%d%H%M%S")[0:10] for it in dict_AF_mission['GMT'] ]
    #uni_hhdd = list(sorted(set( dict_AF_mission['GMT'] )))
    #print('Unique hours of days are: ', uni_hhdd)
    
    return dict_AF_mission


# Comine records of several missions related to a certain day and read them into the memory
def read_USAF(wrf_time, Info_USAF_clt_list, attr_interest):

    print(' ----------- Read USAF obs for the day: ', wrf_time, '--------------')
    time_all = []
    gpsa_all = []
    lat_all = []
    lon_all = []    
    
    for idx in range(0,len(Info_USAF_clt_list)):
        data_one_time = read_USAF_period(wrf_time, Info_USAF_clt_list[idx], attr_interest)
        time_all = time_all + data_one_time['times']
        gpsa_all = gpsa_all + data_one_time['gpsa']
        lat_all = lat_all + data_one_time['lat']
        lon_all = lon_all + data_one_time['lon']

    dict_AF_all_today = {'times':time_all, 'gpsa': gpsa_all, 'lat': lat_all, 'lon':lon_all}

    # Get unique hours where observation exists
    dd_all = []
    hh_all = []
    for time in dict_AF_all_today['times']:
        dd_all.append( time.day )
        hh_all.append( time.hour )
    uni_dd = list(set(dd_all))
    uni_hh = list(set(hh_all))
    print('All of unique hours are: ')
    print(uni_hh)

    # Raise an error in case it happens
    if len(uni_dd) != 1:
        raise ValueError("Current algorithm does not handle. Modify it as needed.")

    return dict_AF_all_today, uni_hh

# Read the obs of a file/Mission only realted to a certain day to the memory
def read_USAF_period(wrf_time, Info_USAF_clt, attr_interest):
    # read the obs only on the day
    filename = Info_USAF_clt.filedir
    beg_time = Info_USAF_clt.beg_time
    end_time = Info_USAF_clt.end_time

    print('Read the obs from file: ', filename)
 
    with open(filename, 'r') as f:
        all_lines = f.readlines()
    
    # process header and filter attributes (line0 - line6)
    print(all_lines[0]) # e.g., Export of 1-second Data for 0215A MARIA at 17:49:38 on 09/18/2017
    str_day = wrf_time

    idx_itt = [] # only collect the attributes of interest
    all_attr = all_lines[5].split() # all of attribute names
    for it in attr_interest:
        idx = [it == it_all for it_all in all_attr]        
        idx_itt.append( int(np.where(idx)[0]) )        

    time = []
    gpsa = []
    lat = []
    lon = []
    # set the time bounds    
    date_format = '%Y%m%d %H:%M:%S'
    begT = str_day +' '+ beg_time
    endT = str_day +' '+ end_time
    print('Time bounds are: ', begT, endT)
    begT_dt = datetime.strptime(begT, date_format)
    endT_dt = datetime.strptime(endT, date_format)    

    for line in all_lines[6:]:
        split_line = line.split()
        # combine day and hour together
        hhmmss = split_line[idx_itt[0]]
        str_time = str_day +' '+ hhmmss
        date_format = '%Y%m%d %H:%M:%S'
        time_dt = datetime.strptime(str_time, date_format)
        if begT_dt <= time_dt <= endT_dt:
            time.append( time_dt )
            gpsa.append( float(split_line[idx_itt[1]-1]) )
            lat.append( float(split_line[idx_itt[2]-1]) ) # GMT Time is supposed to be the same entry
            lon.append( float(split_line[idx_itt[3]-1]) )
        else:
            continue

    dict_AF_one_time = {'times':time, 'gpsa': gpsa, 'lat': lat, 'lon':lon}
    
    # Get unique hours where observation exists
    dd_all = []
    hh_all = []
    for time in dict_AF_one_time['times']:
        dd_all.append( time.day )
        hh_all.append( time.hour )
    uni_dd = list(set(dd_all))
    uni_hh = list(set(hh_all))
    print('Unique hours are: ')
    print(uni_hh)

    # Raise an error in case it happens
    if len(uni_dd) != 1:
        raise ValueError("Current algorithm does not handle. Modify it as needed.") 
    
    return dict_AF_one_time

# check if the USAF file spans two days
def check_timespan( filename, nowtime, nexttime ):
    print('check the file: ', filename)
    with open(filename, 'r') as f:
        all_lines = f.readlines()

    span_days = False
    time_all = []
    for line in all_lines[6:]:
        split_line = line.split()
        hhmmss = split_line[0] 
        time_all.append( hhmmss )
        if hhmmss == '00:00:00':
            span_days = True
    
    if span_days:
       print('The USAF file spans two days...')
       span_days_info =  {nowtime: [time_all[0], '23:59:59'], nexttime: ['00:00:00',time_all[-1]]}
       print(span_days_info)
    else:
       print('The USAF file does not span days...')
       span_days_info =  {nowtime: [time_all[0], time_all[-1]], nexttime: None}
       print(span_days_info)

    return span_days,span_days_info        

# Select records within the range around a certain hour (input)
def mask_time( dict_AF_all, wrf_hh, Nminutes ):
    
    # Get data within the valid time range
    year = dict_AF_all['GMT'][0].year
    month = dict_AF_all['GMT'][0].month
    day= dict_AF_all['GMT'][0].day
    hh = wrf_hh

    time_start = datetime(year,month,day,hh,0) - timedelta(minutes=Nminutes) 
    time_end = datetime(year,month,day,hh,0) + timedelta(minutes=Nminutes)
    
    dict_AF_masked = {}
    # check if the time range exists
    if time_start not in dict_AF_all['GMT']:
        print('AF obs does not exist within the time range at this hour!')
        dict_AF_masked = None
    elif time_end not in dict_AF_all['GMT']:
        print('AF obs does not exist within the time range at this hour!')
        dict_AF_masked = None
    else:
        print('Time starts: '+ datetime.strftime(time_start,"%Y%m%d %H:%M"))
        print('Time ends: '+ datetime.strftime(time_end,"%Y%m%d %H:%M"))
        idx = [time_start <= time_possible <= time_end for time_possible in dict_AF_all['GMT']]
        idx_masktime = np.where(idx)[0]
        for ikey, ivalue in dict_AF_all.items():
            dict_AF_masked[ikey] = np.array(ivalue)[idx_masktime]

    return dict_AF_masked



# Select records within the range around a specific time (input)
def mask_time_of_t( dict_AF_all, uni_dt, Nminutes ):
  
    time_start = uni_dt - timedelta(minutes=Nminutes)
    time_end = uni_dt + timedelta(minutes=Nminutes)
   
    dict_AF_masked = {}
    # check if the time range exists
    if time_start not in dict_AF_all['GMT']:
        print('AF obs does not exist within the time span!')
        dict_AF_masked = None
    elif time_end not in dict_AF_all['GMT']:
        print('AF obs does not exist within the time span!')
        dict_AF_masked = None
    else:
        print('Time starts: '+ datetime.strftime(time_start,"%Y%m%d %H:%M"))
        print('Time ends: '+ datetime.strftime(time_end,"%Y%m%d %H:%M"))
        idx = [time_start <= time_possible <= time_end for time_possible in dict_AF_all['GMT']]
        idx_masktime = np.where(idx)[0]
        for ikey, ivalue in dict_AF_all.items():
            dict_AF_masked[ikey] = np.array(ivalue)[idx_masktime]

    return dict_AF_masked


def Plot_hourly_track( Storm, dict_AF_mission, dict_btk, uni_hhdd, attrs_itt, plot_dir):
     
    lon_min = -100 #-75 #-100 #-71 #np.amin( dict_AF_mission['lon'] )
    lon_max = -85 #-51 #-85 #-57 #np.amax( dict_AF_mission['lon'] )
    lat_min = 15 #11 #15 #10 #np.amin( dict_AF_mission['lat'] )
    lat_max = 31 #25 #31 #20 #np.amax( dict_AF_mission['lat'] )
    gpsa_min = 0
    gpsa_max = np.amax( dict_AF_mission['GPSA'] )

    btk_dt = [datetime.strptime(it_str,"%Y%m%d%H%M") for it_str in dict_btk['time']]

    # Plot figure for each unique hours of days
    Nminutes = 30
    for it_str in uni_hhdd:
        it_dt = datetime.strptime(it_str,"%Y%m%d%H")
        print('Dealing with time: ', it_dt)
        dict_AF_hour = mask_time_of_t( dict_AF_mission, it_dt, Nminutes )
    
        if dict_AF_hour == None:
            print('No AF obs is available at this time span!')
            continue

        # Find the best-track position
        bool_match = [it_dt == it for it in btk_dt]
        if True in bool_match:
            if_btk_exist = True
            idx_btk = np.where( bool_match )[0][0] # the second[0] is due to the possibility of multiple records at the same time
        else:
            if_btk_exist = False

        # figure
        fig = plt.figure()

        ax = plt.subplot(1,1,1,projection=ccrs.PlateCarree())
        ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        ax.coastlines (resolution='10m', color='black', linewidth=1)
        #  Plot the track for the whole mission
        gpsa_min = 0
        gpsa_max = 10000
        ax.scatter(dict_AF_mission['LON'], dict_AF_mission['LAT'], 1, 'grey', cmap='jet', vmin=gpsa_min, vmax=gpsa_max,transform=ccrs.PlateCarree())
        # Plot the hourly track
        AF = ax.scatter(dict_AF_hour['LON'], dict_AF_hour['LAT'], 2, dict_AF_hour['GPSA'], cmap='jet', vmin=gpsa_min, vmax=gpsa_max, transform=ccrs.PlateCarree())
        AF_bar = fig.colorbar(AF,fraction=0.046, pad=0.04) 
        AF_bar.ax.set_ylabel('Flight Height (m)') 
        AF_bar.ax.tick_params(labelsize=7)
        # Mark the best track
        if if_btk_exist:
            ax.scatter(dict_btk['lon'][idx_btk],dict_btk['lat'][idx_btk], 5, 'red', marker='*',transform=ccrs.PlateCarree())
        # Title    
        ax.set_title( datetime.strftime(it_dt, '%Y-%m-%d %H'),fontweight='bold',fontsize=10)
        # Axis labels
        lon_ticks = list(range(math.ceil(lon_min)-2, math.ceil(lon_max)+2,2))
        lat_ticks = list(range(math.ceil(lat_min)-2, math.ceil(lat_max)+2,2))
        gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=False,linewidth=0.1, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.xlabels_bottom = True
        gl.ylabels_left = True
        gl.ylabels_right = False
        gl.ylocator = mticker.FixedLocator(lat_ticks)
        gl.xlocator = mticker.FixedLocator(lon_ticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 6}
        gl.ylabel_style = {'size': 6}

        plt.savefig( plot_dir+Storm+'_'+datetime.strftime(it_dt, '%Y-%m-%d_%H')+'_gpsa.png', dpi=300 )
        print('Saving the figure: ', plot_dir+Storm+'_'+datetime.strftime(it_dt, '%Y-%m-%d_%H')+'.png')
        plt.close()

def track_height( Storm, big_dir, small_dir, plot_dir):

    USAF_list = sorted(glob.glob( small_dir + Storm + '/USAF/2017*.txt' ))

    # define attributes of interest to read
    attrs_itt = ['GMT','GPSA','LAT','LON'] # GMT time / GPS Altimeter (height of the air plane)

    # Plot: track of a mission and supersede hourly obs over the track
    for imission in USAF_list:
        print( '--------- The mission is ', imission, '--------------' )
        # Read the obs from the mission
        mission_head_tail = os.path.split( imission )
        time_mission = mission_head_tail[1][0:8] # e.g., '20170918'
        dict_AF_mission,uni_hhdd = read_USAF_mission( time_mission, imission, attrs_itt )

        # Best-track
        dict_btk = read_bestrack(Storm)

        # Plot
        Plot_hourly_track( Storm, dict_AF_mission, dict_btk, uni_hhdd, attrs_itt, plot_dir )


if __name__ == '__main__':

    Storm = 'HARVEY'
    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'
    plot_dir = small_dir + Storm + '/USAF/Vis/'

    track_height( Storm, big_dir, small_dir, plot_dir) 








































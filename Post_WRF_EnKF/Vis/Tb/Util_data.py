#!/usr/bin/env python3

import numpy as np
import time
import glob

# ------------------------------------------------------------------------------------------------------
#           Object: TCvital
# ------------------------------------------------------------------------------------------------------
# Read the location of the storm from TCvital file
def read_TCvitals(tc_file, DAtime):

    with open(tc_file) as tmp:
        tc_all = tmp.readlines()

    tc_lat = []
    tc_lon = []
    for line in tc_all:
        line_split = line.split()
        tc_time = line_split[3]+line_split[4]

        if tc_time == DAtime:
            print('Time from TCvitals:', tc_time)
            # Read latitude
            if 'N' in line_split[5]:
                tc_lat.append(float(line_split[5].replace('N',''))/10)
            else:
                tc_lat.append( 0-float(line_split[5].replace('S',''))/10)
            # Read longitude
            if 'W' in line_split[6]:
                tc_lon.append(0-float(line_split[6].replace('W',''))/10)
            else:
                tc_lon.append(float(line_split[6].replace('E',''))/10)

            break

    return tc_lon, tc_lat

# ------------------------------------------------------------------------------------------------------
#           Object: Model-IR Tb
# ------------------------------------------------------------------------------------------------------

# Read crtm calcuated IR data from one binary file
def read_simu_Tb_single(Hxb_file, Hxa_file, ch_list):

    xmax = 297
    ymax = 297
    print('Hxb_file:' + Hxb_file)
    print('Hxa_file:' + Hxa_file)
    print('xmax, ymax: '+str(xmax)+' '+str(ymax))

    Hxb_data = np.fromfile(Hxb_file,dtype='<f4') # <: little endian; f: float; 4: 4 bytes
    n_ch = len(Hxb_data)/(xmax*ymax) - 2
    n_ch = int(n_ch)
    if n_ch != len(ch_list):
        print('Error!! # of channels in data is '+str(n_ch))
    Hxb_sim = Hxb_data[:].reshape(n_ch+2,ymax,xmax)

    dict_simu_Tb = {}
    dict_simu_Tb['Lon_x'] = Hxb_sim[0,:,:]
    dict_simu_Tb['Lat_x'] = Hxb_sim[1,:,:]
    dict_simu_Tb['Ch_x'] = ch_list[0]
    dict_simu_Tb['Yb_x'] = Hxb_sim[2,:,:]

    Hxa_data = np.fromfile(Hxa_file,dtype='<f4')
    n_ch = len(Hxa_data)/(xmax*ymax) - 2
    n_ch = int(n_ch)
    if n_ch != len(ch_list):
        print('Error!! # of channels in data is '+str(n_ch))
    Hxa_sim = Hxa_data[:].reshape(n_ch+2,ymax,xmax)
    dict_simu_Tb['Ya_x'] = Hxa_sim[2,:,:]
    #for rec in range(n_ch):
    #    dict_simu_Tb[ch_list[rec]] = sim[rec+2,:,:]

    return dict_simu_Tb

# Read the mean of crtm calcuated IR data from an ensemble of binary files
def read_simu_mean_Hx(Hx_dir, ch_list):

    print("Calculate the mean of the ensemble Tb...")
    start_time=time.process_time()

    # Number of ensemble members
    num_ens = 60
    # Dimension of the domain
    xmax = 297
    ymax = 297

    # List the Yb and Ya files
    file_yb = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_input_mem0*.bin') )
    file_ya = sorted( glob.glob(Hx_dir + '/TB_GOES_CRTM_output_mem0*.bin') )
    # Sanity check the number of calculated ens
    if np.size(file_yb) != num_ens:
        raise ValueError('Wait until all of the ens is calculated!')
    if np.size(file_ya) != num_ens:
        raise ValueError('Wait until all of the ens is calculated!')

    # Read attributes from a member
    tmp_control = np.fromfile( file_yb[0],dtype='<f4') # <: little endian; f: float; 4: 4 bytes
    n_ch = len(tmp_control)/(xmax*ymax) - 2
    n_ch = int(n_ch)
    if n_ch != len(ch_list):
        print('Error!! # of channels in data is '+str(n_ch))
    tmp_data = tmp_control[:].reshape(n_ch+2,ymax,xmax)
    Lon_all =  tmp_data[0,:,:].flatten()  #longitude
    Lat_all = tmp_data[1,:,:].flatten()  #latitude

    # ---- Read prior Tb from the ens ----
    sum_yb = np.zeros( shape=np.shape(Lon_all) )
    # Iterate thru input ens
    for ifile in file_yb:
        print('Reading the file: ' + ifile)
        tmp = np.fromfile( ifile,dtype='<f4')

        # Sanity check
        n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
        tmp_data = tmp[:].reshape(n_ch+2,ymax,xmax)
        sum_yb = sum_yb + tmp_data[2,:,:].flatten()

    Yb_all_mean = sum_yb / num_ens

    # ---- Read posterior Tb from the ens -----
    sum_ya = np.zeros( shape=np.shape(Lon_all) )
    # Iterate thru input ens
    for ifile in file_ya:
        print('Reading the file: ' + ifile)
        tmp = np.fromfile( ifile,dtype='<f4')

        # Sanity check
        n_ch = int( len(tmp_control)/(xmax*ymax) - 2 )
        tmp_data = tmp[:].reshape(n_ch+2,ymax,xmax)
        sum_ya = sum_ya + tmp_data[2,:,:].flatten()

    Ya_all_mean = sum_ya / num_ens
     
    d_mean_Hx = {'Lon_x':Lon_all,'Lat_x':Lat_all,'mean_Hxb':Yb_all_mean,'mean_Hxa':Ya_all_mean}
    return d_mean_Hx


# ------------------------------------------------------------------------------------------------------
#           Object: Model 
# ------------------------------------------------------------------------------------------------------

def read_wrf_domain( wrf_file ):

    print('Read domain info from: ' + wrf_file)
    ncdir = nc.Dataset(wrf_file, 'r')

    Lat_x = ncdir.variables['XLAT'][0,:,:] #latitude: XLAT(time, y, x)
    Lon_x = ncdir.variables['XLONG'][0,:,:] #longitude: XLONG(time, y, x)

    lat_min = np.min( Lat_x.flatten() )
    lat_max = np.max( Lat_x.flatten() )
    lon_min = np.min( Lon_x.flatten() )
    lon_max = np.max( Lon_x.flatten() )

    d_wrf_d = {'lat_min':lat_min, 'lat_max':lat_max, 'lon_min':lon_min, 'lon_max':lon_max}
    return d_wrf_d





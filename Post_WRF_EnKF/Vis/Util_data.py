
from numba import njit, prange
import os,fnmatch # functions for interacting with the operating system
import numpy as np
import time
import glob
import netCDF4 as nc
import math
from math import radians, sin, cos, sqrt, atan2
from datetime import datetime, timedelta
import pandas as pd
from geographiclib.geodesic import Geodesic

# constants
P1000MB = 100000.0
G = 9.81
R_D = 287.0
R = 287.04
CP = 7.0*R_D/2.0
GAMMA = 0.0065
NX = 297
NY = 297
NZ = 42


# ------------------------------------------------------------------------------------------------------
#           Operation: Basic Geography Modelling
# ------------------------------------------------------------------------------------------------------
# Function to convert km to degrees of lat or lon
def km_to_deg(km):
    earth_radius_km = 6371.0
    return km/earth_radius_km * (180.0/np.pi)

# Function to convert degrees to radians
def deg_to_rad(deg):
    return math.radians(deg)

# Function to calculate the distance between two points on Earth(KM)
def mercator_distance(lon1, lat1, lon2, lat2):
    # Convert latitude and longitude from degrees to radians
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    
    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    distance = 6371 * c  # Radius of Earth in kilometers
    return distance

# Function to calculate azimuth & distance from point A to point B
def twoP_inverse(lat1, lon1, lat2, lon2):
    geod = Geodesic.WGS84
    g = geod.Inverse(lat1, lon1, lat2, lon2)
    # s12: distance from the first point to the second in meters
    # azi1: azimuth
    #       ^
    #       | 
    #       |   
    #       A ~ azimuth at the first point 
    #      /
    #    / 
    #  B
    return g['s12']/1000,g['azi1']

# Function to convert from geo azimuth to math definition in Cartesian coordinate
# geo zimuth: North clockwise to South: 0 -> 180; North counter-clockwise to South: 0 -> -180
# math: East counter-clockwise:  0->90->180->270->360 (0)
def azi_geo_to_math(azimuth):
    math_angle = 90 - azimuth
    if math_angle < 0:
        math_angle += 360
    return math_angle

# ------------------------------------------------------------------------------------------------------
#           Operation: Perform Operations
# ------------------------------------------------------------------------------------------------------

# Generate time series constraint by start_date, end_date, and interval
def generate_times( start_time_str, end_time_str, interval ):

    time_diff = datetime.strptime(end_time_str,"%Y%m%d%H%M") - datetime.strptime(start_time_str,"%Y%m%d%H%M")
    time_diff_hour = time_diff.total_seconds() / 3600
    time_interest_dt = [datetime.strptime(start_time_str,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+interval, interval))]
    time_interest_str = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return time_interest_str

# generate linearly-interpolated locations between synoptic scale times
def interpolate_locations( DAtimes, bestrack ):

    # group pairs of synoptic times 
    syno_ts = []
    syno_pair = [[] for i in range(2)]
    for it in bestrack['time']:
        syno_ts.append( datetime.strptime(it, '%Y%m%d%H%M') )
    for i in range(len(syno_ts)-1):
        syno_pair[0].append( syno_ts[i] )
        syno_pair[1].append( syno_ts[i+1] )

    # Linearly interpolate the locations between two synoptic times
    time_hr = []
    lat_hr = []
    lon_hr = []
    for target in DAtimes:
        
        time_hr.append( target )
        if any( hh in target[8:10] for hh in ['00','06','12','18']):
            boolean_compare = [target == itime for itime in bestrack['time']]
            if np.shape(np.where(boolean_compare))[1] == 1:
                idx = int( np.where(boolean_compare)[0] )
            else:
                idx = int( np.where(boolean_compare)[0][0] )
            lat_hr.append( bestrack['lat'][idx] )
            lon_hr.append( bestrack['lon'][idx] )
        else:
            target_t = datetime.strptime(target, '%Y%m%d%H%M')
            boolean_compare = [ syno_pair[0][ip] < target_t < syno_pair[1][ip] for ip in range(len(syno_pair[0]))]
            idx = int( np.where(boolean_compare)[0] )

            # Calculate time differences
            time_diff_total = (syno_pair[1][idx] - syno_pair[0][idx]).total_seconds()
            time_diff_target = (target_t - syno_pair[0][idx]).total_seconds()

            # Calculate interpolation factor
            interp_factor = time_diff_target / time_diff_total

            # Find indices for time list
            boolean_time1 = [ datetime.strftime(syno_pair[0][idx], '%Y%m%d%H%M') == itime for itime in bestrack['time']]
            boolean_time2 = [ datetime.strftime(syno_pair[1][idx], '%Y%m%d%H%M') == itime for itime in bestrack['time']]

            if np.shape(np.where(boolean_time1))[1] == 1:
                idx1 = int( np.where(boolean_time1)[0] )
            else:
                idx1 = int( np.where(boolean_time1)[0][0] )
            if np.shape(np.where(boolean_time2))[1] == 1:
                idx2 = int( np.where(boolean_time2)[0] )
            else:
                idx2 = int( np.where(boolean_time2)[0][0] )

			# Perform linear interpolation for latitude and longitude
            lat_hr.append( bestrack['lat'][idx1] + interp_factor * (bestrack['lat'][idx2] - bestrack['lat'][idx1]) )
            lon_hr.append( bestrack['lon'][idx1] + interp_factor * (bestrack['lon'][idx2] - bestrack['lon'][idx1]) )
    
    dict_bestrack = {'time': time_hr, 'lat': lat_hr, 'lon': lon_hr}

    return dict_bestrack


# Automate the experiment name for one storm
def generate_one_name( Storm,DA,MP ):

    if MP == 'THO':
        if DA == 'IR+MW':
            if Storm == 'HARVEY':
                return 'JerryRun/MW_THO'
            elif Storm == 'IRMA' or Storm == 'JOSE':
                return 'IR+MW-J_DA+J_WRF+J_init-SP-intel17-THO-30hr-hroi900'
            elif Storm == 'MARIA':
                return 'IR+MW-J_DA+J_WRF+J_init-SP-intel17-THO-24hr-hroi900'
            else:
                raise ValueError('No corresponding storm!')
        elif DA == 'IR':
            if Storm == 'HARVEY':
                return 'JerryRun/IR_THO'
            elif Storm == 'IRMA' or Storm == 'JOSE':
                return 'IR-J_DA+J_WRF+J_init-SP-intel17-THO-30hr-hroi900'
            elif Storm == 'MARIA':
                return 'IR-J_DA+J_WRF+J_init-SP-intel17-THO-24hr-hroi900'
            else:
                raise ValueError('No corresponding storm!')
        elif DA == 'CONV':
            if Storm == 'HARVEY':
                return 'J_DA+J_WRF+J_init-Expanse-THO-24hr-hroi300'
            elif Storm == 'IRMA' or Storm == 'JOSE':
                return 'J_DA+J_WRF+J_init-SP-intel17-THO-30hr-hroi900'
            elif Storm == 'MARIA':
                return 'J_DA+J_WRF+J_init-SP-intel17-THO-24hr-hroi900'
            else:
                raise ValueError('No corresponding storm!')

    elif MP == 'WSM6':
        if DA == 'IR+MW':
            if Storm == 'HARVEY':
                return 'JerryRun/MW_WSM6'
            elif Storm == 'IRMA' or Storm == 'JOSE':
                return 'IR+MW-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'
            elif Storm == 'MARIA':
                return 'IR+MW-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900'
            else:
                raise ValueError('No corresponding name!')
        elif DA == 'IR':
            if Storm == 'HARVEY':
                return 'JerryRun/IR_WSM6'
            elif Storm == 'IRMA' or Storm == 'JOSE':
                return 'IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'
            elif Storm == 'MARIA':
                return 'IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900'
            else:
                raise ValueError('No corresponding name!')
        elif DA == 'CONV':
            if Storm == 'HARVEY':
                return 'J_DA+J_WRF+J_init-Expanse-WSM6-24hr-hroi300'
            elif Storm == 'IRMA' or Storm == 'JOSE':
                return 'J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'
            elif Storm == 'MARIA':
                return 'J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900'
            else:
                raise ValueError('No corresponding name!')
    elif MP == 'TuneWSM6':
        #if DA == 'IR+MW':
        #    if Storm == 'HARVEY':
        #        return 'JerryRun/MW_WSM6'
        #    elif Storm == 'IRMA' or Storm == 'JOSE':
        #        return 'IR+MW-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'
        #    elif Storm == 'MARIA':
        #        return 'IR+MW-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900'
        #    else:
        #        raise ValueError('No corresponding name!')
        if DA == 'IR':
            if Storm == 'HARVEY':
                return 'IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi300'
            elif Storm == 'IRMA' or Storm == 'JOSE':
                return 'IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900'
            elif Storm == 'MARIA':
                return 'IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900'
            else:
                raise ValueError('No corresponding name!')
    else:
        raise ValueError('No corresponding MP!')



# Written by wrf-python
def destagger(var, stagger_dim, meta=False):
    """Return the variable on the unstaggered grid.

    This function destaggers the variable by taking the average of the
    values located on either side of the grid box.

    Args:

        var (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A variable
            on a staggered grid.

        stagger_dim (:obj:`int`): The dimension index to destagger.
            Negative values can be used to choose dimensions referenced
            from the right hand side (-1 is the rightmost dimension).

        meta (:obj:`bool`, optional): Set to False to disable metadata and
            return :class:`numpy.ndarray` instead of
            :class:`xarray.DataArray`.  Default is False.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`:
        The destaggered variable.  If xarray is enabled and
        the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.

    """
    var_shape = var.shape
    num_dims = var.ndim
    stagger_dim_size = var_shape[stagger_dim]

    # Dynamically building the range slices to create the appropriate
    # number of ':'s in the array accessor lists.
    # For example, for a 3D array, the calculation would be
    # result = .5 * (var[:,:,0:stagger_dim_size-2]
    #                    + var[:,:,1:stagger_dim_size-1])
    # for stagger_dim=2.  So, full slices would be used for dims 0 and 1, but
    # dim 2 needs the special slice.
    full_slice = slice(None)
    slice1 = slice(0, stagger_dim_size - 1, 1)
    slice2 = slice(1, stagger_dim_size, 1)

    # default to full slices
    dim_ranges_1 = [full_slice] * num_dims
    dim_ranges_2 = [full_slice] * num_dims

    # for the stagger dim, insert the appropriate slice range
    dim_ranges_1[stagger_dim] = slice1
    dim_ranges_2[stagger_dim] = slice2

    result = .5*(var[tuple(dim_ranges_1)] + var[tuple(dim_ranges_2)])
    return result

# ------------------------------------------------------------------------------------------------------
#           Object: Post-storm track
# ------------------------------------------------------------------------------------------------------
# Read the post-storm best-track file of the storm at all times
def read_bestrack( small_dir,Storm ):

    Best_track_path = sorted(fnmatch.filter(os.listdir(small_dir+Storm+'/TC_Guidance/'),'bal*'))
    Best_track_file = small_dir+Storm+'/TC_Guidance/'+Best_track_path[0]
    with open(Best_track_file, 'r') as f:
        all_lines = f.readlines()

    # Process all of records to our format/unit 
    time_all = []
    lat_all = []
    lon_all = []
    maxV_all = []
    minP_all = []
    for line in all_lines:
        # split one record into different parts
        split_line = line.split()
        # Read time
        time_all.append(split_line[2].replace(',','') + '00')
        # Read latitude
        lat_line = split_line[6].replace(',','')
        if 'N' in lat_line:
            lat_all.append(float(lat_line.replace('N',''))/10)
        else:
            lat_all.append(0-float(lat_line.replace('S',''))/10)
        # Read longitute
        lon_line = split_line[7].replace(',','')
        if 'W' in lon_line:
            lon_all.append(0-float(lon_line.replace('W',''))/10)
        else:
            lon_all.append(float(lon_line.replace('E',''))/10)
        # Read max wind
        maxV_all.append(float(split_line[8].replace(',',''))*0.51444) # knots to m/s
        # Read min sea level pressure
        minP_all.append(float(split_line[9].replace(',',''))) # mb

    dict_bestrack = {'time': time_all, 'lat': lat_all, 'lon': lon_all, 'max_ws': maxV_all, 'min_slp': minP_all}
    return dict_bestrack

# Return records only at times of interest from the complete best-track file 
def btk_in_duration(Storm, Btk_start, Btk_end, hour_step):

    # Calculate the duration from the start to the end of deterministic forecast in hour
    Btk_diff = datetime.strptime(Btk_end,"%Y%m%d%H%M") - datetime.strptime(Btk_start,"%Y%m%d%H%M")
    Btk_diff_hour = Btk_diff.total_seconds() / 3600
    # List the times (in string format) every 6 hours in the duration 
    time_interest_dt = [datetime.strptime(Btk_start,"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(Btk_diff_hour+1), hour_step))]
    time_interest_str = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    # Read the complete best-track file
    dict_all_btk = read_bestrack(small_dir,Storm)
    # Get the indices in the best-track file corresponded to the times of interest
    idx_in_btk = []

    lat_btk = []
    lon_btk = []
    max_ws_btk = []
    min_slp_btk = []
    for time_str in time_interest_str:
        boolean_compare = [t_btk == time_str for t_btk in dict_all_btk['time'][:]]
        if any(boolean_compare):
            idx = int(np.where(boolean_compare)[0][0])
            idx_in_btk.append(idx)
            lat_btk.append(dict_all_btk['lat'][idx])
            lon_btk.append(dict_all_btk['lon'][idx])
            max_ws_btk.append(dict_all_btk['max_ws'][idx])
            min_slp_btk.append(dict_all_btk['min_slp'][idx])

    dict_btk = {'time': time_interest_str, 'lat': lat_btk, 'lon': lon_btk, 'max_ws': max_ws_btk, 'min_slp': min_slp_btk}

    return dict_btk


# ------------------------------------------------------------------------------------------------------
#           Object: Real-Time TC Guidance Project
# ------------------------------------------------------------------------------------------------------
# read model simulations (a-deck texts)
#def read_adeck( Storm,model_name ):
 
# subset forecasts that were generated for a specific model
def subset_adeck( Storm,model_name ):   
    
    # specify input and output names    
    input_path = sorted(fnmatch.filter(os.listdir(small_dir+Storm+'/TC_Guidance/'),'aal*'))
    input_name = small_dir+Storm+'/TC_Guidance/'+input_path[0]
    output_name = small_dir+Storm+'/TC_Guidance/'+model_name+'_'+input_path[0]
    # subset
    with open(input_name, 'r') as file, open(output_name, 'w') as output_file:
        for line in file:
            if model_name in line:
                output_file.write(line)

# read forecasts from a model
def read_adeck_model( Storm,model_name ):
    
    # Read all lines
    aal_name = sorted(fnmatch.filter(os.listdir(small_dir+Storm+'/TC_Guidance/'),'aal*'))
    file_name = small_dir+Storm+'/TC_Guidance/'+model_name+'_'+aal_name[0]
    with open(file_name,'r') as f:
        all_lines = f.readlines()

    # Read the file into a pandas DataFrame
    df = pd.read_csv( file_name,delimiter=',',header=None, dtype=str )
    # Find out all possible values
    all_fc_init_times = df[2].values # get all records of forecast initialization time
    fc_init_times = list(set(all_fc_init_times))
    fc_init_times = [it.replace(" ","") for it in fc_init_times if " " in it ]#remove space
    all_fc_lens =  df[5].values # get all records of forecast lengths (in hours)
    fc_lens = list(set(all_fc_lens))
    fc_lens = [int(it) for it in fc_lens]

    # Set up dictionaries
    fc = {}
    for init in fc_init_times:
        fc[init] = {}
        # set up times
        times_dt = [datetime.strptime(init,"%Y%m%d%H") + timedelta(hours=t) for t in fc_lens] 
        times_str = [it.strftime("%Y%m%d%H") for it in times_dt]
        fc[init]['time'] = times_str
        fc[init]['lat'] = np.full(np.shape(times_str),np.nan)
        fc[init]['lon'] = np.full(np.shape(times_str),np.nan)
        fc[init]['MSLP'] = np.full(np.shape(times_str),np.nan)
        fc[init]['Vmax'] = np.full(np.shape(times_str),np.nan)
    
    # Gather forecasted attributes 
    for init in fc_init_times:
        lat_all = []
        lon_all = []
        mslp_all = []
        vmax_all = []
        # Read forecasted attributes for a certain forecast
        for line in all_lines:
            if init in line:    
                # split one record into different parts
                split_line = line.split()
                # Read latitude
                lat_line = split_line[6].replace(',','')
                if 'N' in lat_line:
                    lat_all.append(float(lat_line.replace('N',''))/10)
                else:
                    lat_all.append(0-float(lat_line.replace('S',''))/10)
                # Read longitute
                lon_line = split_line[7].replace(',','')
                if 'W' in lon_line:
                    lon_all.append(0-float(lon_line.replace('W',''))/10)
                else:
                    lon_all.append(float(lon_line.replace('E',''))/10)
                # Read max wind
                vmax_all.append(float(split_line[8].replace(',',''))*0.51444) # knots to m/s
                # Read min sea level pressure
                mslp_all.append(float(split_line[9].replace(',',''))) # mb
            else:
                continue 
        # fill in
        fc[init]['lat'] = np.array(lat_all)
        fc[init]['lon'] = np.array(lon_all)
        fc[init]['mslp'] = np.array(mslp_all)
        fc[init]['vmax'] = np.array(vmax_all)

    # Sanity check
    #print(fc['2020082718']['time'])
    #print(fc['2020082718']['vmax'])

    return fc


# ------------------------------------------------------------------------------------------------------
#           Object: TCvital
# ------------------------------------------------------------------------------------------------------
# Read the location of the storm from TCvital file
def read_TCvitals(small_dir, Storm, DAtime):
 
    tc_file = small_dir+Storm+'/TCvitals/'+Storm+'_tcvitals'
    with open(tc_file) as tmp:
        tc_all = tmp.readlines()

    for line in tc_all:
        line_split = line.split()
        tc_time = line_split[3]+line_split[4]

        if tc_time == DAtime:
            # Read latitude
            if 'N' in line_split[5]:
                tc_lat = float(line_split[5].replace('N',''))/10
            else:
                tc_lat =  0-float(line_split[5].replace('S',''))/10
            # Read longitude
            if 'W' in line_split[6]:
                tc_lon = 0-float(line_split[6].replace('W',''))/10
            else:
                tc_lon = float(line_split[6].replace('E',''))/10
            # Read slp
            tc_slp = float(line_split[9]) 
            break
    return tc_lon, tc_lat, tc_slp

# ------------------------------------------------------------------------------------------------------
#           Object: Model-IR Tb
# ------------------------------------------------------------------------------------------------------

# For ungridded points, find ones that are inside an circle area with a center
def find_circle_area_ungrid( ct_lon, ct_lat, lon_location, lat_location, radius_th ):

    lg_location = []
    for im in range(len(lon_location)):
        distance = mercator_distance(ct_lon, ct_lat, lon_location[im], lat_location[im]) 
        lg_location.append( distance <= radius_th )
    # Get the indices where the elements are True
    index_area = [index for index, value in enumerate( lg_location ) if value]
    return index_area

# Read crtm calcuated IR data from one binary file
def read_simu_IR_single(Hxb_file, Hxa_file, ch_list):

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
def read_IRsimu_mean_Hx(Hx_dir, ch_list):

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

@njit(parallel=True)
# compute temperature in kelvin
def compute_tk( full_p, theta ):

    p_flat = full_p.flatten()
    theta_flat = theta.flatten()
    res = np.zeros( (len(theta_flat)) )
    res[:] = np.nan
    for im in prange( len(theta_flat) ):
        pi = (p_flat[im]/P1000MB)**(R_D/CP)
        res[im] = pi*theta_flat[im]
    # make sure all values are reasonable
    assert res.any() != np.nan
    return res

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

def calculate_convergence(u, v, dx, dy):
    """
    Calculate horizontal wind convergence.
    
    Parameters:
    - u: 3d array of east-west wind components (m/s)
    - v: 3d array of north-south wind components (m/s)
    - dx: grid spacing in the east-west direction (m)
    - dy: grid spacing in the north-south direction (m)
    
    Returns:
    - convergence: array of horizontal wind convergence (1/s)
    """
    # Calculate the spatial derivatives of u and v
    dudx, dudy = np.gradient(u, [dx, dy])#, axis=(2, 1, 0))
    dvdx, dvdy = np.gradient(v, [dx, dy])#, axis=(2, 1, 0))
    
    # Calculate horizontal wind convergence
    convergence = dudx + dvdy
    print(np.shape(convergence))

    return convergence


# Find the grids that are inside an circle area with a center
def find_circle_area_model_ij( wrf_file, ct_i, ct_j, radius_threshold, dx):

    # ct_i, ct_j: in grid space with the index starting from 0
    # radius_threshold, dx: in km

    ncdir = nc.Dataset( wrf_file )
    xlat = ncdir.variables['XLAT'][0,:,:]
    nx = xlat.shape[1]
    ny = xlat.shape[0]

    i_inArea = []
    j_inArea = []
    for i in range( nx ):
        for j in range( ny ):
            radius = (((i-ct_i)**2+(j-ct_j)**2)**0.5)*dx
            if radius > radius_threshold:
                continue
            else:
                i_inArea.append( i )
                j_inArea.append( j )

    Idx_inArea = []
    for io in range(len(i_inArea)):
        Idx_inArea.append( int(j_inArea[io]*nx+i_inArea[io]) )
  
    #         
    #(0,ny) ____________
    #   ^  |____________
    #   |  |____________
    #   |  |____________
    #      |____________
    #   (0,0)  -->     (nx,0)

    # check 
    #print('Location of the center: i-'+ str(ct_i) +' j-'+str(ct_j))
    #print('Location of the first selected grid: i-'+ str(i_inArea[0]) +' j-'+str(j_inArea[0]))
    #print('Location of the last selected grid: i-'+ str(i_inArea[-1]) +' j-'+str(j_inArea[-1]))
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # netCDF4 package inherits the indices starting rule from WRF. 
    # Indices start from lower-left corner:
    # lower-left -> lower-right == 0 -> nx in python == increasing longitude 
    # lower-left -> higher-left == 0 -> ny in python == increasing latitude
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    return Idx_inArea

# *************************************************************
#             Compute sea level pressure
# *************************************************************
@njit(parallel=True)
def least_zeta_level( full_p,PCONST ):

    level = np.zeros( (NX*NY) )
    level[:] = np.nan
    p_flat = full_p.reshape( full_p.shape[0],-1)
    for im in prange( len(level) ):
        k = 0
        found = False
        while (found is False) and (k<NZ):
            if p_flat[k,im] < (p_flat[0,im]-PCONST):
                level[im] = k
                found = True
            k = k+1

    if level.any() == np.nan:
        raise ValueError('Error in finding 100 hPa up!')
    return level    

@njit(parallel=True)
def t_surf_slp( level_flat,p,t,z,q,PCONST ):

    p_flat = p.reshape( p.shape[0],-1)
    t_flat = t.reshape( t.shape[0],-1)
    z_flat = z.reshape( z.shape[0],-1)
    q_flat = q.reshape( q.shape[0],-1)

    t_surf = np.zeros( (NX*NY) )
    t_surf[:] = np.nan
    t_sea_level = np.zeros( (NX*NY) )
    t_sea_level[:] = np.nan

    for im in prange( NX*NY):

        klo = max( level_flat[im]-1,0 )
        khi = min( klo+1,NZ-2 )
        if klo == khi:
            raise ValueError('Error trapping levels')

        plo = p_flat[klo,im]
        phi = p_flat[khi,im]
        tlo = t_flat[klo,im]*(1+0.608*q_flat[klo,im])
        thi = t_flat[khi,im]*(1+0.608*q_flat[khi,im])
        zlo = z_flat[klo,im]
        zhi = z_flat[khi,im]

        p_at_pconst = p_flat[0,im] - PCONST
        t_at_pconst = thi - (thi-tlo)*np.log(p_at_pconst/phi)*np.log(plo/phi)
        z_at_pconst = zhi - (zhi-zlo)*np.log(p_at_pconst/phi)*np.log(plo/phi)

        t_surf[im] = t_at_pconst * (p_flat[0,im]/p_at_pconst)**(GAMMA*R/G)
        t_sea_level[im] = t_at_pconst + GAMMA*z_at_pconst

    assert t_surf.any() != np.nan
    assert t_sea_level.any() != np.nan
    return t_surf, t_sea_level

# Rewrite the f_computeslp function from WRF in python
def compute_slp( ncdir ):

    ridiculous_mm5_test = True
    # specific constants for assumptions made in this routine:
    TC = 273.16+17.5
    PCONST=10000    

    # Read data
    # full pressure
    p = ncdir.variables['P'][0,:,:,:] # perturbation
    pb = ncdir.variables['PB'][0,:,:,:]
    full_p = p + pb 
    # full T in kelvin
    theta = ncdir.variables['T'][0,:,:,:] # theta perturbation
    full_theta = theta + 300 
    full_t = compute_tk( full_p, full_theta )
    full_t = full_t.reshape( (theta.shape[0],theta.shape[1],theta.shape[2]) )
    # geopotential height
    ph = ncdir.variables['PH'][0,:,:,:] # perturbation
    phb = ncdir.variables['PHB'][0,:,:,:]
    full_ph = (ph+phb)/G # in meter
    destag_ph = destagger(full_ph, -3) # destagger ph from nLevel=43 to nLevel=42
    # qvapor
    qvapor = ncdir.variables['QVAPOR'][0,:,:,:] 
    
    # prepare
    z = destag_ph
    t = full_t
    p = full_p
    q = qvapor

    # Find least zeta level that is PCONST Pa above the surface.  We
    # later use this level to extrapolate a surface pressure and
    # temperature, which is supposed to reduce the effect of the diurnal
    # heating cycle in the pressure field.
    level_flat = least_zeta_level( full_p, PCONST )
    level_flat = level_flat.astype(int) # convert from float array to int array

    # Get temperature PCONST Pa above surface.  Use this to extrapolate
    # the temperature at the surface and down to sea level.
    t_surf, t_sea_level = t_surf_slp( level_flat,p,t,z,q,PCONST )

    # If we follow a traditional computation, there is a correction to the
    # sea level temperature if both the surface and sea level
    # temperatures are *too* hot.
    if ridiculous_mm5_test:
        #print('Correction when both surface and sea level temperature are too hot!')
        for im in range(NX*NY):        
            l1 = t_sea_level[im] < TC
            l2 = t_surf[im] <= TC
            l3 = not l1
            if (l2 and l3):
                t_sea_level[im] = TC
        else:
                t_sea_level[im] = TC - 0.005*(t_surf[im]-TC)**2

    #     The grand finale: ta da!
    sea_level_pressure = np.zeros( (NX*NY) )
    sea_level_pressure[:] = np.nan
    z = z.reshape( z.shape[0],-1)
    p = p.reshape( p.shape[0],-1)
    for im in range(NX*NY):
        z_half_lowest = z[0,im]
        sea_level_pressure[im] = 0.01 * (p[0,im]*math.exp((2.0*G*z_half_lowest)/(R*(t_sea_level[im]+t_surf[im])))) #hPa
    assert sea_level_pressure.any() != np.nan

    #index = np.where( sea_level_pressure < 900 )[0]
    #lala = np.zeros( (NX*NY) )
    #for im in range(NX*NY):
    #    z_half_lowest = z[0,im]
    #    lala[im] = 0.01 * p[0,im]*math.exp((2.0*G*z_half_lowest)/(R*(t_sea_level[im]+t_surf[im])))
    sea_level_pressure = sea_level_pressure.reshape( (NX,NY) )
    return sea_level_pressure



if __name__ == '__main__':
   
    small_dir = '/expanse/lustre/projects/pen116/zuy121/Pro5_BP_OSSE/' 
    start_time=time.process_time()
    storm = 'Laura'
    read_adeck_model( storm,'HWRF' )
    end_time = time.process_time()
    print ('time needed: ', end_time-start_time, ' seconds')


































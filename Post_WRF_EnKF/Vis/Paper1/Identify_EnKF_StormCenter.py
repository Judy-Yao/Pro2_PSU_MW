
import numpy as np
import EnKF_minSLP_track as SC
import netCDF4 as nc
from datetime import datetime, timedelta

import Util_data as UD




# Generate time series
def generate_times( Storms, start_time_str, end_time_str, interval ):

    dict_times = {}
    for istorm in Storms:
        time_diff = datetime.strptime(end_time_str[istorm],"%Y%m%d%H%M") - datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M")
        time_diff_hour = time_diff.total_seconds() / 3600
        time_interest_dt = [datetime.strptime(start_time_str[istorm],"%Y%m%d%H%M") + timedelta(hours=t) for t in list(range(0, int(time_diff_hour)+interval, interval))]
        dict_times[istorm] = [time_dt.strftime("%Y%m%d%H%M") for time_dt in time_interest_dt]
    return dict_times

#----------------------------------------------------------------------------
# Maximum 10-meter wind speed
#----------------------------------------------------------------------------
def find_Vmax( wrfout ):

    ncdir = nc.Dataset( wrfout )
    # Geographical 
    lat = ncdir.variables['XLAT'][0,:,:]
    lon = ncdir.variables['XLONG'][0,:,:]
    # 10-meter wind
    U10 = ncdir.variables['U10'][0,:,:]
    V10 = ncdir.variables['V10'][0,:,:]
    windspeed = (U10 ** 2 + V10 ** 2) ** 0.5
    vmax_value = np.amax(windspeed)
    # Find the index of the maximum element
    max_index = np.unravel_index(np.argmax(windspeed), windspeed.shape) # e.g., (124, 143)
    # Identify the location of vmax
    vmax_lon = lon[max_index]
    vmax_lat = lat[max_index]

    return [vmax_lat,vmax_lon,vmax_value]

#----------------------------------------------------------------------------
# Read and process modelled data
#----------------------------------------------------------------------------

# Obtain min slp from WRF output
def Get_model_HPI( big_dir,istorm,iExper,DAtimes ):

    d_HPI = {}
    minslp_lat = []
    minslp_lon = []
    minslp_value = []
    vmax_lat = []
    vmax_lon = []
    vmax_value = []

    for DAtime in DAtimes:
        if analysis:
            x_dir = big_dir+istorm+'/'+iExper+'/fc/'+DAtime+'/wrf_enkf_output_d03_mean'
        else:
            x_dir = big_dir+istorm+'/'+iExper+'/fc/'+DAtime+'/wrf_enkf_input_d03_mean'
        print('Reading the EnKF mean from ', x_dir)
        # Read min slp
        list_slp = SC.find_minSLP( istorm,x_dir,DAtime )
        minslp_lat.append( list_slp[0] )
        minslp_lon.append( list_slp[1] )
        minslp_value.append( list_slp[2] )
        # Read max 10-meter wind
        list_vmax = find_Vmax( x_dir )
        vmax_lat.append( list_vmax[0] )
        vmax_lon.append( list_vmax[1] )
        vmax_value.append( list_vmax[2] )

    d_HPI = {'times':DAtimes,'mslp_lon':minslp_lon,'mslp_lat':minslp_lat,
             'mslp':minslp_value,'vmax_lon':vmax_lon,'vmax_lat':vmax_lat,'vmax':vmax_value}

    return d_HPI

#----------------------------------------------------------------------------
# Write into txt files
#----------------------------------------------------------------------------
def write_HPI( ist,d_HPI,save_dir ):

    DAtimes = d_HPI['times']
    mslp_lon =  [ "{0:.3f}".format(item) for item in d_HPI['mslp_lon'] ]
    mslp_lat = [ "{0:.3f}".format(item) for item in d_HPI['mslp_lat'] ]
    mslp = [ "{0:.3f}".format(item) for item in d_HPI['mslp'] ]
    vmax_lon =  [ "{0:.3f}".format(item) for item in d_HPI['vmax_lon'] ] #longitude
    vmax_lat = [ "{0:.3f}".format(item) for item in d_HPI['vmax_lat'] ] #latitude
    vmax = [ "{0:.3f}".format(item) for item in d_HPI['vmax'] ] 

    # Stack each list into an array
    all_attrs = np.column_stack( (DAtimes,mslp_lon,mslp_lat,mslp,vmax_lon,vmax_lat,vmax) )

    # ---- Write to file and save it to the disk ----
    header = ['Time','mslp_lon','mslp_lat','mslp(hpa)','vmax_lon','vmax_lat','vmax(ms-1)']
    if analysis:
        file_name = save_dir+'HPI_wrf_enkf_output_d03_mean.'+start_time_str[ist]+'_'+end_time_str[ist]+'.txt'
    else:
        file_name = save_dir+'HPI_wrf_enkf_input_d03_mean.'+start_time_str[ist]+'_'+end_time_str[ist]+'.txt'
    # Define the fixed width for each field
    fixed_widths = [12+3,8+3,8+3,9+3,8+3,8+3,10+3]
    # Function to pad strings to fixed width
    def pad_string(string, width):
        return string.ljust(width)

    with open(file_name,'w') as file:
        # Write the header
        padded_header = [pad_string(field, width) for field, width in zip(header, fixed_widths)]
        file.write("".join(padded_header) + "\n")

        # Write the records
        for record in all_attrs :
            padded_record = [pad_string(field, width) for field, width in zip(record, fixed_widths)]
            file.write("".join(padded_record) + "\n")

        # Add a success flag
        file.write('Success\n')

    print('Saved '+file_name)


if __name__ == '__main__':

    big_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/'
    small_dir = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'

    # ------------Configuration------------------
    Storms = ['HARVEY','IRMA','JOSE','MARIA']
    MP = ['WSM6','THO']
    DA = ['CONV','IR','IR+MW']

    start_time_str = {'HARVEY':'201708231200','IRMA':'201709040000','JOSE':'201709060000','MARIA':'201709170000'}
    end_time_str = {'HARVEY':'201708241200','IRMA':'201709050000','JOSE':'201709070000','MARIA':'201709180000'}

    # action
    analysis = False

    # -------------------------------------------

    # Create experiment names
    Exper_names = {}
    for istorm in Storms:
        Exper_names[istorm] = {}
        for imp in MP:
            Exper_names[istorm][imp] = {}
            for ida in DA:
                Exper_names[istorm][imp][ida] = UD.generate_one_name( istorm,ida,imp )

    # Identify DA times in the period of interest
    d_hrs = generate_times( Storms, start_time_str, end_time_str, 1 )
    d_6hrs = generate_times( Storms, start_time_str, end_time_str, 6 )

    # Obtain modeled HPI
    for istorm in Storms:
        for imp in MP:
            for ida in DA:
                d_HPI = Get_model_HPI( big_dir,istorm,Exper_names[istorm][imp][ida],d_hrs[istorm] )
                save_dir = small_dir+istorm+'/'+Exper_names[istorm][imp][ida]+'/Data_analyze/'
                write_HPI( istorm,d_HPI,save_dir )












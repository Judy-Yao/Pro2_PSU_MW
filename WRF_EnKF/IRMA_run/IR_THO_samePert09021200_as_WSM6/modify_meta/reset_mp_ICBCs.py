import glob
from netCDF4 import Dataset

# modify the metadata "MP_PHYSICS" for single WRF file
def modify_wrf_meta( file_path ):

    # Open the file
    nc_file = Dataset(file_path, "r+")

    # Check if attribute exists and modify it
    if "MP_PHYSICS" in nc_file.ncattrs():
        nc_file.setncattr("MP_PHYSICS", 8)  # Set the new value: THO
        print("MP_PHYSICS changed to 6")
    else:
        print("MP_PHYSICS not found. Adding it...")
        nc_file.setncattr("MP_PHYSICS", 8)

    # Close the file
    nc_file.close()

# modify files under fc
def modify_fc():

    fc_path = data_dir+'fc/201709021200/'
    fc_files = sorted( glob.glob(fc_path + '/wrf*') )

    for file in fc_files:
        print(file+'...')
        modify_wrf_meta( file )

# modify files under rc
def modify_rc():

    rc_path = data_dir+'rc/201709021200/'
    rc_files = sorted( glob.glob(rc_path + '/wrf*') )

    for file in rc_files:
        print(file+'...')
        modify_wrf_meta( file )


if __name__ == '__main__':

    data_dir = '/scratch/06191/tg854905/Pro2_PSU_MW/IRMA/CONV_THO_samePert09021200_as_WSM6/'
    
    # modify files under /fc
    modify_fc()

    # modify files under /rc
    modify_rc()

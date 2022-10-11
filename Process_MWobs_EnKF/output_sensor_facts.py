#!/usr/bin/env python3

# ================================================================================================================================
# Output necessary sensor parameters to a HDF5 file. 
# ================================================================================================================================
# The microwave-observation-preprocessing system (MOPS) requires a lookp table sort of thing for parameters related to sensors.
# This process, unfortunately, is impossible to automate at this point (June, 2022). Thus, it is the responsibility of the user who  
# uses this system to eyeball the technical details of the sensor and revise this program before running the system.
# ---------------------------------------------------------------------------------------------------------------------------------

# !!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# This file may need to be tweaked a bit by adding
# or removing the sensors/channels of the user's interest.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ----------------- Important Details -----------------------------------------
# Sensor name should be consistent with main_process_MWobs.m 
# Channel index/number should be consistent with the ones used in CRTM package
# Other parameters are gathered from their official sources
#------------------------------------------------------------------------------


import h5py 
import numpy as np

f = h5py.File('sensor_database.HDF5','w')

# Group: AMSR2
amsr2 = f.create_group("AMSR2")
amsr2.attrs['ScanType'] = 'Conical'
amsr2.attrs['FovType'] = 'IFOV'
# --- low f
amsr2_ch = amsr2.create_group("18.7GHzV-Pol") 

amsr2_ch_dset = amsr2_ch.create_dataset("Channel_num", data=7)

amsr2_ch_dset = amsr2_ch.create_dataset("fovs_alongTrack", data=22)
amsr2_ch_dset.attrs['Units'] = 'km'

amsr2_ch_dset = amsr2_ch.create_dataset("fovs_crossTrack", data=14)
amsr2_ch_dset.attrs['Units'] = 'km'

amsr2_lf_nsample = 243 #Number of Samples for low frequency in one scan.
amsr2_lf_scan = np.full((1,amsr2_lf_nsample),47.57) #Assume the earth is spherical; Use laws of cosine
amsr2_ch_dset = amsr2_ch.create_dataset("scan_angles", data=amsr2_lf_scan)
amsr2_ch_dset.attrs['Units'] = 'degree'
amsr2_ch_dset.attrs['Dimension'] = '1,nsamples'

amsr2_ch_dset = amsr2_ch.create_dataset("max_scan_angle", data=47.57)
amsr2_ch_dset.attrs['Units'] = 'degree'

# --- high f
amsr2_ch = amsr2.create_group("89GHzV-PolA-Scan")

amsr2_ch_dset = amsr2_ch.create_dataset("Channel_num", data=13)

amsr2_ch_dset = amsr2_ch.create_dataset("fovs_alongTrack", data=5)
amsr2_ch_dset.attrs['Units'] = 'km'

amsr2_ch_dset = amsr2_ch.create_dataset("fovs_crossTrack", data=3)
amsr2_ch_dset.attrs['Units'] = 'km'

amsr2_hf_nsample = 972 #combine A scan and B scan
amsr2_hf_scan = np.full((1,amsr2_hf_nsample),47.57)
amsr2_ch_dset = amsr2_ch.create_dataset("scan_angles", data=amsr2_hf_scan)
amsr2_ch_dset.attrs['Units'] = 'degree'
amsr2_ch_dset.attrs['Dimension'] = '1,nsamples'

amsr2_ch_dset = amsr2_ch.create_dataset("max_scan_angle", data=47.57)
amsr2_ch_dset.attrs['Units'] = 'degree'

# ----
amsr2_ch = amsr2.create_group("89GHzV-PolB-Scan")

amsr2_ch_dset = amsr2_ch.create_dataset("Channel_num", data=14)

amsr2_ch_dset = amsr2_ch.create_dataset("fovs_alongTrack", data=5)
amsr2_ch_dset.attrs['Units'] = 'km'

amsr2_ch_dset = amsr2_ch.create_dataset("fovs_crossTrack", data=3)
amsr2_ch_dset.attrs['Units'] = 'km'

amsr2_hf_nsample = 972 #combine A scan and B scan
amsr2_hf_scan = np.full((1,amsr2_hf_nsample),47.57)
amsr2_ch_dset = amsr2_ch.create_dataset("scan_angles", data=amsr2_hf_scan)
amsr2_ch_dset.attrs['Units'] = 'degree'
amsr2_ch_dset.attrs['Dimension'] = '1,nsamples'

amsr2_ch_dset = amsr2_ch.create_dataset("max_scan_angle", data=47.57)
amsr2_ch_dset.attrs['Units'] = 'degree'


# Group: ATMS
atms = f.create_group("ATMS")
atms.attrs['ScanType'] = 'Cross-track'
atms.attrs['FovType'] = 'EFOV'

atms_ch = atms.create_group("183.31+-7GHzQH-Pol")

atms_ch_dset = atms_ch.create_dataset("Channel_num", data=18)

atms_ch_dset = atms_ch.create_dataset('fovs_alongTrack',data=1.1)
atms_ch_dset.attrs['Units'] = 'degree'

atms_ch_dset = atms_ch.create_dataset('fovs_crossTrack',data=2.2)
atms_ch_dset.attrs['Units'] = 'degree'

atms_nsample = 96
atms_scan = np.linspace(-52.725,52.725,atms_nsample)
atms_ch_dset = atms_ch.create_dataset("scan_angles", data=atms_scan)
atms_ch_dset.attrs['Units'] = 'degree'
atms_ch_dset.attrs['Dimension'] = '1,nsamples'

atms_ch_dset = atms_ch.create_dataset("max_scan_angle", data=52.725)
atms_ch_dset.attrs['Units'] = 'degree'


# Group: GMI
gmi = f.create_group("GMI")
gmi.attrs['ScanType'] = 'Conical'
gmi.attrs['FovType'] = 'EFOV'

# --- low frequency
gmi_ch = gmi.create_group("18.7GHzV-Pol")

gmi_ch_dset = gmi_ch.create_dataset("Channel_num", data=3)

gmi_ch_dset = gmi_ch.create_dataset("fovs_alongTrack", data=10.9)
gmi_ch_dset.attrs["Units"] = 'km'

gmi_ch_dset = gmi_ch.create_dataset("fovs_crossTrack", data=18.1)
gmi_ch_dset.attrs["Units"] = 'km'

gmi_lf_nsample = 221
gmi_lf_scan = np.full((1,gmi_lf_nsample),48.5)
gmi_ch_dset = gmi_ch.create_dataset("scan_angles", data=gmi_lf_scan)
gmi_ch_dset.attrs["Units"] = 'degree'
gmi_ch_dset.attrs['Dimension'] = '1,nsamples'

gmi_ch_dset = gmi_ch.create_dataset("max_scan_angle", data=48.5)
gmi_ch_dset.attrs["Units"] = 'degree'

# --- high frequency
gmi_ch = gmi.create_group("183.31+-7GHzV-Pol")

gmi_ch_dset = gmi_ch.create_dataset("Channel_num", data=13)

gmi_ch_dset = gmi_ch.create_dataset("fovs_alongTrack", data=3.8)
gmi_ch_dset.attrs["Units"] = 'km'

gmi_ch_dset = gmi_ch.create_dataset("fovs_crossTrack", data=5.8)
gmi_ch_dset.attrs["Units"] = 'km'

gmi_hf_nsample = 221
gmi_hf_scan = np.full((1,gmi_hf_nsample),48.5)
gmi_ch_dset = gmi_ch.create_dataset("scan_angles", data=gmi_hf_scan)
gmi_ch_dset.attrs["Units"] = 'degree'
gmi_ch_dset.attrs['Dimension'] = '1,nsamples'

gmi_ch_dset = gmi_ch.create_dataset("max_scan_angle", data=48.5)
gmi_ch_dset.attrs["Units"] = 'degree'


# Group: MHS
mhs = f.create_group("MHS")
mhs.attrs['ScanType'] = 'Cross-track'
mhs.attrs['FovType'] = 'IFOV'

mhs_ch = mhs.create_group("190.31GHzV-Pol")

mhs_ch_dset = mhs_ch.create_dataset("Channel_num", data=5)

mhs_ch_dset = mhs_ch.create_dataset('fovs_alongTrack', data=1.1)
mhs_ch_dset.attrs['Units'] = 'degree'

mhs_ch_dset = mhs_ch.create_dataset('fovs_crossTrack', data=1.1)
mhs_ch_dset.attrs['Units'] = 'degree'

mhs_nsample = 90
mhs_scan = np.linspace(-49.5,49.5,mhs_nsample)
mhs_ch_dset = mhs_ch.create_dataset("scan_angles", data=mhs_scan)
mhs_ch_dset.attrs['Units'] = 'degree'
mhs_ch_dset.attrs['Dimension'] = '1,nsamples'

mhs_ch_dset = mhs_ch.create_dataset("max_scan_angle", data=49.5)
mhs_ch_dset.attrs['Units'] = 'degree'


# Group: SAPHIR
saphir = f.create_group("SAPHIR")
saphir.attrs['ScanType'] = 'Cross-track'
saphir.attrs['FovType'] = 'IFOV'

saphir_ch = saphir.create_group("183.31+-6.8GHz")

saphir_ch_dset = saphir_ch.create_dataset("Channel_num", data=5)

saphir_ch_dset = saphir_ch.create_dataset('fovs_alongTrack', data=0.6616)
saphir_ch_dset.attrs['Units'] = 'degree'

saphir_ch_dset = saphir_ch.create_dataset('fovs_crossTrack', data=0.6616)
saphir_ch_dset.attrs['Units'] = 'degree'

saphir_nsample = 182
saphir_scan = np.linspace(-42.96,42.96,saphir_nsample)
saphir_ch_dset = saphir_ch.create_dataset("scan_angles", data=saphir_scan)
saphir_ch_dset.attrs['Units'] = 'degree'
saphir_ch_dset.attrs['Dimension'] = '1,nsamples'

saphir_ch_dset = saphir_ch.create_dataset("max_scan_angle", data=42.96)
saphir_ch_dset.attrs['Units'] = 'degree'


# Group: SSMI
ssmi = f.create_group("SSMI")
ssmi.attrs['ScanType'] = 'Conical'
ssmi.attrs['FovType'] = 'EFOV'

# --- low frequency
ssmi_ch = ssmi.create_group("fcdr_tb19v")

ssmi_ch_dset = ssmi_ch.create_dataset("Channel_num", data=1)

ssmi_ch_dset = ssmi_ch.create_dataset("fovs_alongTrack", data=69)
ssmi_ch_dset.attrs['Units'] = 'km'

ssmi_ch_dset = ssmi_ch.create_dataset("fovs_crossTrack", data=43)
ssmi_ch_dset.attrs['Units'] = 'km'

ssmi_lf_nsample = 64
ssmi_lf_scan = np.full((1,ssmi_lf_nsample), 45)
ssmi_ch_dset = ssmi_ch.create_dataset("scan_angles", data=ssmi_lf_scan)
ssmi_ch_dset.attrs["Units"] = 'degree'
ssmi_ch_dset.attrs['Dimension'] = '1,nsamples'

ssmi_ch_dset = ssmi_ch.create_dataset('max_scan_angle', 45)
ssmi_ch_dset.attrs['Units'] = 'degree'

# --- high frequency
ssmi_ch = ssmi.create_group("fcdr_tb85v")

ssmi_ch_dset = ssmi_ch.create_dataset("Channel_num", data=6)

ssmi_ch_dset = ssmi_ch.create_dataset("fovs_alongTrack", data=15)
ssmi_ch_dset.attrs['Units'] = 'km'

ssmi_ch_dset = ssmi_ch.create_dataset("fovs_crossTrack", data=13)
ssmi_ch_dset.attrs['Units'] = 'km'

ssmi_hf_nsample = 128
ssmi_hf_scan = np.full((1,ssmi_hf_nsample), 45)
ssmi_ch_dset = ssmi_ch.create_dataset("scan_angles", data=ssmi_hf_scan)
ssmi_ch_dset.attrs["Units"] = 'degree'
ssmi_ch_dset.attrs['Dimension'] = '1,nsamples'

ssmi_ch_dset = ssmi_ch.create_dataset('max_scan_angle', 45)
ssmi_ch_dset.attrs['Units'] = 'degree'

# Group: SSMIS
ssmis = f.create_group("SSMIS")
ssmis.attrs['ScanType'] = 'Conical'
ssmis.attrs['FovType'] = 'IFOV'

# --- low frequency
ssmis_ch = ssmis.create_group("19.35GHzV-Pol")

ssmis_ch_dset = ssmis_ch.create_dataset("Channel_num", data=13)

ssmis_ch_dset = ssmis_ch.create_dataset("fovs_alongTrack", data=73)
ssmis_ch_dset.attrs['Units'] = 'km'

ssmis_ch_dset = ssmis_ch.create_dataset("fovs_crossTrack", data=47)
ssmis_ch_dset.attrs['Units'] = 'km'

ssmis_lf_nsample = 90
ssmis_lf_scan = np.full((1,ssmis_lf_nsample), 45)
ssmis_ch_dset = ssmis_ch.create_dataset("scan_angles", data=ssmis_lf_scan)
ssmis_ch_dset.attrs["Units"] = 'degree'
ssmis_ch_dset.attrs['Dimension'] = '1,nsamples'

ssmis_ch_dset = ssmis_ch.create_dataset('max_scan_angle', 45)
ssmis_ch_dset.attrs['Units'] = 'degree'

# --- high frequency
ssmis_ch = ssmis.create_group("183.31+-6.6GHzH-Pol")

ssmis_ch_dset = ssmis_ch.create_dataset("Channel_num", data=9)

ssmis_ch_dset = ssmis_ch.create_dataset("fovs_alongTrack", data=14)
ssmis_ch_dset.attrs['Units'] = 'km'

ssmis_ch_dset = ssmis_ch.create_dataset("fovs_crossTrack", data=13)
ssmis_ch_dset.attrs['Units'] = 'km'

ssmis_hf_nsample = 180
ssmis_hf_scan = np.full((1,ssmis_hf_nsample), 45)
ssmis_ch_dset = ssmis_ch.create_dataset("scan_angles", data=ssmis_hf_scan)
ssmis_ch_dset.attrs["Units"] = 'degree'
ssmis_ch_dset.attrs['Dimension'] = '1,nsamples'

ssmis_ch_dset = ssmis_ch.create_dataset('max_scan_angle', 45)
ssmis_ch_dset.attrs['Units'] = 'degree'



f.close()

# Header


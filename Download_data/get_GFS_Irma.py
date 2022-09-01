#!/usr/bin/env python
#################################################################
# Python Script to retrieve 3 online Data files of 'ds084.1',
# total 617.28M. This script uses 'requests' to download data.
#
# Highlight this script by Select All, Copy and Paste it into a file;
# make the file executable and run it on command line.
#
# You need pass in your password as a parameter to execute
# this script; or you can set an environment variable RDAPSWD
# if your Operating System supports it.
#
# Contact rdahelp@ucar.edu (RDA help desk) for further assistance.
#################################################################


import sys, os
import requests
import datetime as dt


def check_file_status(filepath, filesize):
    sys.stdout.write('\r')
    sys.stdout.flush()
    size = int(os.stat(filepath).st_size)
    percent_complete = (size/filesize)*100
    sys.stdout.write('%.3f %s' % (percent_complete, '% Completed'))
    sys.stdout.flush()

# Try to get password
if len(sys.argv) < 2 and not 'RDAPSWD' in os.environ:
    try:
        import getpass
        input = getpass.getpass
    except:
        try:
            input = raw_input
        except:
            pass
    pswd = input('Password: ')
else:
    try:
        pswd = sys.argv[1]
    except:
        pswd = os.environ['RDAPSWD']

url = 'https://rda.ucar.edu/cgi-bin/login'
values = {'email' : 'yao.zhu.91@gmail.com', 'passwd' : pswd, 'action' : 'login'}
# Authenticate
ret = requests.post(url,data=values)
if ret.status_code != 200:
    print('Bad Authentication')
    print(ret.text)
    exit(1)
    
dspath = 'https://rda.ucar.edu/data/ds084.1/'

# --------- Download data for the spin up --------- 
spinup_start_str = ['2017090218',] #!!!!!!!!!!!!!!!!!!!!! 
file_id = ['f' + str(id).zfill(3) for id in list(range(0, 13, 6)) ]
spinup_list = []

for date in spinup_start_str:
    yyyy = date[0:4] 
    yyyymmdd = date[0:8]
    for id in file_id:
        spinup_list.append(yyyy+'/'+yyyymmdd+'/gfs.0p25.'+ date + '.' + id + '.grib2')

for file in spinup_list:
    filename=dspath+file
    file_base = os.path.basename(file)
    print('Downloading',file_base)
    req = requests.get(filename, cookies = ret.cookies, allow_redirects=True, stream=True)
    filesize = int(req.headers['Content-length'])
    with open(file_base, 'wb') as outfile:
        chunk_size=1048576
        for chunk in req.iter_content(chunk_size=chunk_size):
            outfile.write(chunk)
            if chunk_size < filesize:
                check_file_status(file_base, filesize)
    check_file_status(file_base, filesize)
    print()
    
# --------- Download data for WRF-DA cyle --------
# Set up times using GFS files
t_startDA_str = '2017090306' #!!!!!!!!!!!!!!!!!!!!!
t_endDA_str = '2017090318'  #!!!!!!!!!!!!!!!!!!!!!

DA_diff = dt.datetime.strptime(t_endDA_str, '%Y%m%d%H') - dt.datetime.strptime(t_startDA_str, '%Y%m%d%H') 
DA_diff_hour = DA_diff.total_seconds() / 3600

DAtimes_dt = [dt.datetime.strptime(t_startDA_str, '%Y%m%d%H')+dt.timedelta(hours=t) for t in list(range(0, int(DA_diff_hour)+1,6))]
DAtimes_str = [time_dt.strftime('%Y%m%d%H') for time_dt in DAtimes_dt]

# Set up file ids
DAfile_id = ['f' + str(id).zfill(3) for id in list(range(0, 145, 6)) ]

DAcycle_list = []

for date in DAtimes_str:
    yyyy = date[0:4] 
    yyyymmdd = date[0:8]
    for id in DAfile_id:
        DAcycle_list.append(yyyy+'/'+yyyymmdd+'/gfs.0p25.'+ date + '.' + id + '.grib2')

for file in DAcycle_list:
    filename=dspath+file
    file_base = os.path.basename(file)
    print('Downloading',file_base)
    req = requests.get(filename, cookies = ret.cookies, allow_redirects=True, stream=True)
    filesize = int(req.headers['Content-length'])
    with open(file_base, 'wb') as outfile:
        chunk_size=1048576
        for chunk in req.iter_content(chunk_size=chunk_size):
            outfile.write(chunk)
            if chunk_size < filesize:
                check_file_status(file_base, filesize)
    check_file_status(file_base, filesize)
    print()

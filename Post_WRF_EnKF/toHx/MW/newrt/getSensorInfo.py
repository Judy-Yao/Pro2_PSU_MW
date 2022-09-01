#!/usr/bin/env python3

import sys

def main(MW_File):
    # Read the content inside the microwave obs to a list
    with open(MW_File) as f:
        all_lines = f.readlines() 

    # Declare an empty list (Datatype list in python is very similar to datatype cell in matlab)
    sensorCh_rep = [] 
    sensor_rep = []

    # For each string (elements in this list), separate it into words
    for line in all_lines:
        split_all_lines = line.split() 
        sensor_rep.append(split_all_lines[1])
        # Combine sensor name and channel number together 
        sensorCh_rep.append(split_all_lines[1] + ' ' + split_all_lines[2])

    # Find the uqiue sensor / sensor-and-channel combination
    sensor_uni_set = set(sensor_rep)
    sensor_uni = list(sensor_uni_set)

    sensorCh_uni_set = set(sensorCh_rep)
    sensorCh_uni = list(sensorCh_uni_set)

    # For each combination, separate it into sensor name and channel number
    sensor_Ch = []
    for record in sensorCh_uni: 
        sensor_Ch.append(record.split())

    Ch_perSS = []
    # Loop through each unique sensor
    for iss in sensor_uni:
        for ir in sensor_Ch:      
            Ch_perSS.append([]) 
            if sensor_Ch[sensor_Ch.index(ir)][0] == iss:
                Ch_perSS[sensor_uni.index(iss)].append(sensor_Ch[sensor_Ch.index(ir)][1])
            else:
                continue

    # Build a dictionary: sensor = channel1, channel2, ...
    dict_sensor = {}
    for iss in sensor_uni:
        dict_sensor[sensor_uni[sensor_uni.index(iss)]] = Ch_perSS[sensor_uni.index(iss)]

    # Convert dictionary to list in order to easily pass values to bash variables
    List_SS_Ch = [ [] for _ in range(len(sensor_uni))]

    for iss in sensor_uni:
        List_SS_Ch[sensor_uni.index(iss)].append(iss)
        for ich in Ch_perSS[sensor_uni.index(iss)]:
            List_SS_Ch[sensor_uni.index(iss)].append(ich)
        

    return List_SS_Ch


if __name__ == '__main__':
    Storm = sys.argv[1]
    Exper_name = sys.argv[2]
    MW_time = sys.argv[3]
    MW_File = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/'+Storm+'/Obs_y/MW/microwave_d03_'+MW_time+'_so';
    List_SS_Ch = main(MW_File)
    with open('/scratch/06191/tg854905/Pro2_PSU_MW/'+Storm+'/'+Exper_name+'/fc/'+MW_time+'/'+MW_time+'_sensorCh','w') as f:
        for isCh in List_SS_Ch:
            f.write(" ".join(isCh))
            f.write('\n')

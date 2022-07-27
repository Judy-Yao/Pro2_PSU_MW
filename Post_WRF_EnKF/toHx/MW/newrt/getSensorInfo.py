#!/usr/bin/env python3

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

    return dict_sensor

if __name__ == '__main__':
    MW_File = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/HARVEY/Obs_y/MW/microwave_d03_201708221200_so';
    dict_sensor = main(MW_File)
    print(len(dict_sensor))
    print(dict_sensor)

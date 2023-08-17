% ==========================================================================================
% Find the sensor factual parameters from facts_'sensor'.HDF5 database
% Use sensor and channel name as identifiers
% ==========================================================================================

function [scantype,Ch_num,fov_alongTrack,fov_crossTrack,max_scan_angle,scan_angles] = SensorInfo_read(sensor,ChName)

    % Process channel names which consist of "/"
    % Some channel names read from L1C MW files consist of "/" while channel names used in facts_'sensor'.HDF5 don't contain "/" because it is sort of "not allowed" in HDF5 file writing.
    switch(ChName)
        case '183.31+/-7GHzV-Pol'
            ChName = '183.31+-7GHzV-Pol';
        case '183.31+/-6.6GHzH-Pol'
            ChName = '183.31+-6.6GHzH-Pol';
        case '183.31+/-6.8GHz'
            ChName = '183.31+-6.8GHz';
    end
        
    % Read sensor information from the HDF5 file
    scantype = h5readatt('sensor_database.HDF5',"/"+sensor,'ScanType');
    search_name = "/" + sensor + "/" + ChName;
    
    Ch_num = h5read("sensor_database.HDF5",search_name + '/Channel_num'); % Channel number/index that is consistent with CRTM coefficients
    max_scan_angle = h5read("sensor_database.HDF5",search_name + '/max_scan_angle');
    scan_angles = h5read("sensor_database.HDF5",search_name + '/scan_angles');
    % conical scan: FOV parameters to be read is applicable to all pixels; 
    % cross-track scan: FOV parameters to be read is only applicable to pixel at nadir
    fov_alongTrack = h5read("sensor_database.HDF5",search_name + '/fovs_alongTrack');
    fov_crossTrack = h5read("sensor_database.HDF5",search_name + '/fovs_crossTrack');

end

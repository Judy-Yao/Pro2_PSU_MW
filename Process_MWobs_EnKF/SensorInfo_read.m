function [Ch_num,fov_alongTrack,fov_crossTrack,max_scan_angle,scan_angles] = SensorInfo_read(sensor,ChName)

    % Process channel names which consist of slash
    switch(ChName)
        case '183.31+/-7GHzV-Pol'
            ChName = '183.31+-7GHzV-Pol';
        case '183.31+/-6.6GHzH-Pol'
            ChName = '183.31+-6.6GHzH-Pol';
        case '183.31+/-6.8GHz'
            ChName = '183.31+-6.8GHz';
    end
        
    search_name = '/'+ sensor +'/'+ ChName;
    
    Ch_num = h5read("sensor_database.HDF5",search_name + '/Channel_num'); % Channel number/index that is consistent with CRTM coefficients
    fov_alongTrack = h5read("sensor_database.HDF5",search_name + '/fovs_alongTrack');
    fov_crossTrack = h5read("sensor_database.HDF5",search_name + '/fovs_crossTrack');
    max_scan_angle = h5read("sensor_database.HDF5",search_name + '/max_scan_angle');
    scan_angles = h5read("sensor_database.HDF5",search_name + '/scan_angles');

end
function [myTimes,mySat_name,myChNum,myLat,myLon,myTb,myROI_hydro,myROI_other,myObsErr,mySat_alt] = ProduceForEnKF(DAtime, Ch, min_WRF_lat, max_WRF_lat, min_WRF_lon, max_WRF_lon, Tb_file, control)

    % Read GOES obs 
    [lat_col,lon_col,Tb_col] = Read_GOES(Tb_file);

    % Narrow obs area to the domain of interest
    idx_d3 = (lat_col >= min_WRF_lat) & (lat_col <= max_WRF_lat) & (lon_col >= min_WRF_lon) & (lon_col <= max_WRF_lon);
    lat_d3 = lat_col(idx_d3); lon_d3 = lon_col(idx_d3); Tb_d3 = Tb_col(idx_d3);
    
    % Locate the start point for each ROI plan
    filter_grid_step = control.filter_reso / control.dx;
    grid_start(2) = floor(filter_grid_step(2) / 2);
    grid_start(1) = grid_start(2) + .5*(2*filter_grid_step(2) - filter_grid_step(1));
    
    % ---------- Process obs for each ROI plan --------------
    % Preallocating memory
    myTimes = cell(size(control.roi_oh));
    mySat_name = cell(size(control.roi_oh));
    myChNum = cell(size(control.roi_oh));
    myLat = cell(size(control.roi_oh));
    myLon = cell(size(control.roi_oh));
    myTb  = cell(size(control.roi_oh));
    myROI_hydro = cell(size(control.roi_oh));
    myROI_other = cell(size(control.roi_oh));
    myObsErr = cell(size(control.roi_oh));
    mySat_alt = cell(size(control.roi_oh));

    for iroi = 1:length(control.roi_oh) % ROI first references [200,0] second [30,30]
        obs_idx = grid_start(iroi):filter_grid_step(iroi):length(lat_d3);
        lat = lat_d3(obs_idx);
        lon = lon_d3(obs_idx);
        Tb = Tb_d3(obs_idx);
        % randomize IR records for this ROI
        randOrder = randperm(length(cat(1,Tb(:))));
        
        myTimes{iroi} = repmat([DAtime+"00"], length(Tb), 1);
        mySat_name{iroi} = repmat("abi_gr", length(Tb), 1);
        myChNum{iroi} = repmat(Ch, length(Tb), 1);  
        myLat{iroi} = lat(randOrder); myLon{iroi} = lon(randOrder); myTb{iroi} = Tb(randOrder);
        myROI_hydro{iroi} = repmat(control.roi_oh{iroi}(2), length(Tb), 1); 
        myROI_other{iroi} = repmat(control.roi_oh{iroi}(1), length(Tb), 1);
        myObsErr{iroi} =  repmat(control.obsError, length(Tb), 1);   
        mySat_alt{iroi} = repmat(control.Sat_alt, length(Tb), 1);    
    end
    











end

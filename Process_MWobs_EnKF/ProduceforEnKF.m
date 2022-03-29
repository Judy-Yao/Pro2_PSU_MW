function [sat_name,myLat,myLon,myTb,mySat_lat,mySat_lon,mySat_alt,mySat_azimuth,myScan_angle,myZenith_angle,myFov_crossTrack,myFov_alongTrack,myTimes,myChNum,myRoi_hydro,myRoi_otherVars,myObsErr] = ProduceforEnKF(iTb,Swath_used,ChIdx_ps,ChName_ps,if_swath_good,DAtime_ps,loc_DAtime_ps,Tb_file,control) 
    % ---------------------------------------------------------------------
    % ---- For each raw observation file
    % ---- Loop through each possibly good channel AND
    % ---- Read characteristics of sensor into memory
    % ---------------------------------------------------------------------
    
    % Preallocating memory
    lat = cell(size(Swath_used{iTb}));
    lon = cell(size(Swath_used{iTb}));
    Tb = cell(size(Swath_used{iTb}));
    zenith_angle = cell(size(Swath_used{iTb}));
    sat_lat = cell(size(Swath_used{iTb}));
    sat_lon = cell(size(Swath_used{iTb}));
    sat_alt = cell(size(Swath_used{iTb}));
    sat_azimuth = cell(size(Swath_used{iTb}));
    outime = cell(size(Swath_used{iTb}));

    % Get platform and sensor information from the new name 
    [filepath,filename,filext] = fileparts(Tb_file);
    ss_info = split(filename,'.'); platform = ss_info(2); sensor = ss_info(3);

    % modify satellite-and-sensor name so that it is consistent with what is used in the CRTM package 
    if (contains(platform,'GCOMW1'))
        sat_name = "amsr2_gcom-w1";
    elseif (contains(platform,'NPP'))
        sat_name = "atms_npp";
    elseif (contains(platform,'GPM'))
        sat_name = "gmi_gpm";
    elseif (contains(platform,'METOPA'))
        sat_name = "mhs_metop-a";
    elseif (contains(platform,'METOPB'))
        sat_name = "mhs_metop-b";
    elseif (contains(platform,'NOAA18'))
        sat_name = "mhs_n18";
    elseif (contains(platform,'NOAA19'))
        sat_name = "mhs_n19";
    elseif (contains(platform,'MT1'))
        sat_name = "saphir_meghat";
    elseif (contains(platform,'F16'))
        sat_name = "ssmis_f16";
    elseif (contains(platform,'F17'))
        sat_name = "ssmis_f17";
    elseif (contains(platform,'F18'))
        sat_name = "ssmis_f18";
    end

    % read data from L1C file
    for it = 1:length(Swath_used{iTb}) % it: item
        % quality control
        if if_swath_good{iTb}(it) == 0
           continue;
        end
        
        % **latitude**
        lat{it} = h5read(Tb_file,Swath_used{iTb}(it) + '/Latitude'); % npixel, nscan

        % **longitude**
        lon{it} = h5read(Tb_file,Swath_used{iTb}(it) + '/Longitude'); % npixel, nscan

        % **Tb**
        Tb_allCh = h5read(Tb_file,Swath_used{iTb}(it) + '/Tc');
        Tb{it} = squeeze(Tb_allCh(ChIdx_ps{iTb}(it),:,:)); % npixel, nscan

        % **zenith angle**: the angle of the satellite from the local zenith as seen at the pixel locatio on the earth
        % Interpretation on nChUIA
        % For example, GMI sensor has two swaths. Channel 1-9 are under
        % swath 1 and channel 10-14 are under swath 2, which respectively
        % corresponds to nChUIA1 and nChUIA2. If nChUIA1 = 1, at a specific (a pixel on a scan)
        % location, only one value of incidence angle exist for all 9 channels; 
        % If nChUIA1 = n (n >= 2), at a specific location, each channel
        % could in theory reference n values.
        incidenceAngles = h5read(Tb_file,Swath_used{iTb}(it) + '/incidenceAngle'); % nChUIA1, npixel, nscan
        iAIndices = h5read(Tb_file,Swath_used{iTb}(it) + '/incidenceAngleIndex'); % nchannel, nscan
        if size(incidenceAngles,1) == 1
            zenith_angle{it}(:,:) = squeeze(incidenceAngles(1,:,:));
        else
            num_Ch_perSW = size(iAIndices,1);
            num_unique_perCh = zeros(1,num_Ch_perSW); % size(iAIndices,1): number of channels under this swath
            for ich = 1:size(iAIndices,1)
                num_unique_perCh(ich) = length(unique(iAIndices(ich,:)));
            end
            if sum(num_unique_perCh) == num_Ch_perSW
                iAIndices = iAIndices(:,1); % For a channel, all scans are references to the same set of zenith angles
                zenith_angle{it}(:,:) = squeeze(incidenceAngles(iAIndices(ChIdx_ps{iTb}(it)),:,:)); % npixel, nscan
            else
                disp("Error: current algorithm does not work!! Please modify it.");
            end
        end

        % **latitude of satellite scan**
        sat_lat{it} = h5read(Tb_file,Swath_used{iTb}(it)+ '/SCstatus/SClatitude'); % nscan

        % **longitude of satellite scan**
        sat_lon{it} = h5read(Tb_file,Swath_used{iTb}(it)+ '/SCstatus/SClongitude'); % nscan

        % **altitude of satellite scan**
        sat_alt{it} = h5read(Tb_file,Swath_used{iTb}(it)+'/SCstatus/SCaltitude'); % nscan (unit: km)

        % **azimuth of satellite scan (used in CRTM)**: the angle subtended
        % by the horizontal projection of a direct line from the satellite
        % to the FOV and the North-South axis measured cloclwise from North
        % (0-> 360 degrees)
        %                     North
        %                      |AZ /
        %                      |~ FOV
        %                      | /
        %                      |/
        % West ---------- Satellite ------------- East
        perscan_length{it} = size(Tb{it},1); % number of pixels per scan
        sat_azimuth{it} = geodetic2aer(lat{it},                                       lon{it},                                       0, ...
                                             repmat(sat_lat{it}',[perscan_length{it} 1]), repmat(sat_lon{it}',[perscan_length{it} 1]), repmat(sat_alt{it}'*1000,[perscan_length{it} 1]), ...
                                             referenceEllipsoid('WGS 84'));

        % **Time** 
        year   = double(h5read(Tb_file, Swath_used{iTb}(it) + '/ScanTime/Year')); % nscan
        month  = double(h5read(Tb_file, Swath_used{iTb}(it) + '/ScanTime/Month')); % nscan
        day    = double(h5read(Tb_file, Swath_used{iTb}(it) + '/ScanTime/DayOfMonth')); % nscan
        hour   = double(h5read(Tb_file, Swath_used{iTb}(it) +  '/ScanTime/Hour')); % nscan
        minute = double(h5read(Tb_file, Swath_used{iTb}(it) +  '/ScanTime/Minute')); % nscan
        second = double(h5read(Tb_file, Swath_used{iTb}(it) +  '/ScanTime/Second')); % nscan
        
        num_scans= size(Tb{it},2); % number of scans

        for i_time = 1:num_scans
            datetime_temp = datetime(year(i_time), month(i_time), day(i_time), hour(i_time), minute(i_time), second(i_time));
            if ( second(i_time) >= 30)
                datetime_temp = datetime_temp + minutes(1);
            end
            outime{it}(i_time,1) = string(datestr(datetime_temp,'yyyymmddHHMM')); % nscan
        end

    end

    % Special treatment to AMSR2 89GHz
    if (sensor == "AMSR2") & control.comnine_AMSR89GHz
        [Swath_used,ChIdx_ps,ChName_ps,lat,lon,Tb,zenith_angle,sat_lat,sat_lon,sat_alt,sat_azimuth,outime] = Combine_AMSR2(iTb,Swath_used,ChIdx_ps,ChName_ps,DAtime_ps,lat,lon,Tb,zenith_angle,sat_lat,sat_lon,sat_alt,sat_azimuth,outime);
    elseif (sensor == "AMSR2") & (~control.comnine_AMSR89GHz)
        disp(['Two 89 GHz scans on AMSR2 are not combined!']);
    end


    % Time cost to get obs_index for one Channel: 200 ~ 300 seconds 
    % ---------------------------------------------------------------------
    % ---- Define area of interest for simulation AND
    % ---- Select the raw obs for every grid point for EnKF assimilation
    % ---------------------------------------------------------------------
    % Prepare area: best-track location followed
    nx = control.nx*control.domain_buffer; % zoom out
    ny = control.ny*control.domain_buffer; % zoom out
    min_XLAT = loc_DAtime_ps{iTb}{1}(1) - (ny/2*control.dx)/(cos(loc_DAtime_ps{iTb}{1}(1)*(pi/180))*111);
    max_XLAT = loc_DAtime_ps{iTb}{1}(1) + (ny/2*control.dx)/(cos(loc_DAtime_ps{iTb}{1}(1)*(pi/180))*111);
    min_XLONG = loc_DAtime_ps{iTb}{1}(2) - (nx/2*control.dx)/111;
    max_XLONG = loc_DAtime_ps{iTb}{1}(2) + (nx/2*control.dx)/111;
    disp(['min of xlong: ',num2str(min_XLONG), ', max of xlong: ',num2str(max_XLONG)]);
    disp(['min of xlat: ',num2str(min_XLAT), ', max of xlat: ',num2str(max_XLAT)]);
    latitudes  = linspace(min_XLAT,max_XLAT,ny);
    longitudes = linspace(min_XLONG,max_XLONG,nx);
    [XLAT, XLONG] = meshgrid(latitudes,longitudes);
    
    % Separate gird points of low frequency from of high frequency
    % Original grid points:         Filtered grid points:
    %       * * * * *                *   *   * 
    %       * * * * *                          
    %       * * * * *                *   *   * 
    %       * * * * *                          
    %       * * * * *                *   *   * 
    slots_x = cell(size(Swath_used{iTb}));
    slots_y = cell(size(Swath_used{iTb}));
    obs_index = cell(size(Swath_used{iTb}));
        
    filter_ratio_grid = control.filter_ratio / control.dx; 
    grid_start_hf =  floor(filter_ratio_grid(2) / 2); % high frequency
    grid_start_lf = grid_start_hf + .5*(2*filter_ratio_grid(2) - filter_ratio_grid(1)); % low frequency

    for it = 1:length(Swath_used{iTb})
        if (contains(ChName_ps{iTb}(it), "18.7GHzV-Pol")) || (contains(ChName_ps{iTb}(it), "19.35GHzV-Pol"))
            slots_x{it} = grid_start_lf:filter_ratio_grid(1):nx;
            slots_y{it} = grid_start_lf:filter_ratio_grid(1):ny;
            obs_index{it} = PickRawforCRTM(lat{it},lon{it},Tb{it},min_XLONG,max_XLONG,min_XLAT,max_XLAT,latitudes,longitudes,slots_x{it},slots_y{it},control);
        else
            slots_x{it} = grid_start_hf:filter_ratio_grid(2):nx;
            slots_y{it} = grid_start_hf:filter_ratio_grid(2):ny;
            obs_index{it} = PickRawforCRTM(lat{it},lon{it},Tb{it},min_XLONG,max_XLONG,min_XLAT,max_XLAT,latitudes,longitudes,slots_x{it},slots_y{it},control);
        end
    end
    
    % ---------------------------------------------------------------------
    % ---- Select MW records for EnKF assimilation
    % ---------------------------------------------------------------------
    
    % Preallocating memory
    DA_lat = cell(size(Swath_used{iTb}));
    DA_lon = cell(size(Swath_used{iTb}));
    DA_Tb = cell(size(Swath_used{iTb}));

    DA_sat_lat = cell(size(Swath_used{iTb}));
    DA_sat_lon = cell(size(Swath_used{iTb}));
    DA_sat_alt = cell(size(Swath_used{iTb}));
    DA_sat_azimuth = cell(size(Swath_used{iTb})); 
    DA_scan_angle = cell(size(Swath_used{iTb}));
    DA_zenith_angle = cell(size(Swath_used{iTb}));
    DA_fov_crossTrack = cell(size(Swath_used{iTb}));
    DA_fov_alongTrack = cell(size(Swath_used{iTb}));

    DA_times = cell(size(Swath_used{iTb}));
    DA_chNum = cell(size(Swath_used{iTb}));
    DA_obsError = cell(size(Swath_used{iTb}));

    % Assign values
    for it = 1:length(Swath_used{iTb}) % it: item
        obs_index_array = obs_index{it};
        obs_index_1d = obs_index_array(obs_index_array(:) == obs_index_array(:)); % get rid of obs_index with value NaN

        allTb = reshape(Tb{it},[],1);
        all_lat = reshape(lat{it},[],1);
        all_lon = reshape(lon{it},[],1);
        DA_Tb{it}          = allTb(obs_index_1d);
        DA_lat{it}         = all_lat(obs_index_1d);
        DA_lon{it}         = all_lon(obs_index_1d);

        [scan_position, scan_num] = ind2sub(size(Tb{it}),obs_index_1d); % size(Tb{it}): npixel, nscan; output from ind2sub: row, column ??

        DA_sat_lat{it} = sat_lat{it}(scan_num);
        DA_sat_lon{it} = sat_lon{it}(scan_num);
        DA_sat_alt{it} = sat_alt{it}(scan_num);
        all_sat_azimuth = reshape(sat_azimuth{it},[],1);
        DA_sat_azimuth{it} = all_sat_azimuth(obs_index_1d);
        all_zenith_angles = reshape(zenith_angle{it},[],1);
        DA_zenith_angle{it} = all_zenith_angles(obs_index_1d);

        [scantype,ch_num,fov_alongTrack,fov_crossTrack,max_scan_angle,scan_angles] = SensorInfo_read(sensor,ChName_ps{iTb}{it});
        DA_scan_angle{it} = scan_angles(scan_position)';
        [DA_fov_crossTrack{it}, DA_fov_alongTrack{it}] = Get_pixel_resolution(scantype,ch_num,fov_alongTrack,fov_crossTrack,DA_lat{it},DA_lon{it},DA_zenith_angle{it},DA_sat_lat{it},DA_sat_lon{it},DA_sat_alt{it});
    
         for my_scan_num_idx = 1:length(scan_num)
             DA_times{it}(my_scan_num_idx,1) = outime{it}(scan_num(my_scan_num_idx));
         end
        DA_chNum{it} = ones(numel(obs_index_1d),1,'int64')*ch_num;
        DA_obsError{it} = ones(numel(obs_index_1d),1)*control.obsError(it);
    end
    clear Tb lat lon sat_lat sat_lon sat_alt sat_azimuth zenith_angle scan_angles outime 
    clear allTb all_lat all_lon all_sat_azimuth all_zenith_angles 

    % ---------------------------------------------------------------------
    % ---- Produce MW records for EnKF assimilation
    % ---------------------------------------------------------------------

    myLat = cell(size(control.roi_oh));
    myLon = cell(size(control.roi_oh));
    myTb = cell(size(control.roi_oh));

    mySat_lat = cell(size(control.roi_oh));
    mySat_lon = cell(size(control.roi_oh));
    mySat_alt = cell(size(control.roi_oh));
    mySat_azimuth = cell(size(control.roi_oh)); 
    myScan_angle = cell(size(control.roi_oh));
    myZenith_angle = cell(size(control.roi_oh));
    myFov_crossTrack = cell(size(control.roi_oh));
    myFov_alongTrack = cell(size(control.roi_oh));

    myTimes = cell(size(control.roi_oh));
    myChNum = cell(size(control.roi_oh));
    myRoi_hydro = cell(size(control.roi_oh));
    myRoi_otherVars = cell(size(control.roi_oh));
    myObsErr = cell(size(control.roi_oh));

    for ir = 1:length(control.roi_oh)
        % randomize MW records for this ROI 
        randOrder = randperm(length(cat(1,DA_Tb{:}))); % DA_Tb may contain 1 channel (low or high) or contain both low and high channels

        tem_lat = cat(1,DA_lat{:});                 myLat{ir} = tem_lat(randOrder); %clear tem_lat
        tem_lon = cat(1,DA_lon{:});                 myLon{ir} = tem_lon(randOrder); %clear tem_lon
        tem_Tb = cat(1,DA_Tb{:});                   myTb{ir} = tem_Tb(randOrder); %clear tem_Tb
        
        tem_sat_lat = cat(1,DA_sat_lat{:}); mySat_lat{ir} = tem_sat_lat(randOrder); %clear tem_sat_lat
        tem_sat_lon = cat(1,DA_sat_lon{:}); mySat_lon{ir} = tem_sat_lon(randOrder); %clear tem_sat_lon
        tem_sat_alt = cat(1,DA_sat_alt{:}); mySat_alt{ir} = tem_sat_alt(randOrder); %clear tem_sat_alt
        tem_sat_azimuth = cat(1,DA_sat_azimuth{:}); mySat_azimuth{ir} = tem_sat_azimuth(randOrder); %clear tem_sat_azimuth
        tem_scan_angle = cat(2,DA_scan_angle{:}); myScan_angle{ir} = tem_scan_angle(randOrder)'; %clear tem_scan_angle
        tem_zenith_angle = cat(1,DA_zenith_angle{:}); myZenith_angle{ir} = tem_zenith_angle(randOrder); %clear  tem_zenith_angle
        tem_fov_crossTrack = cat(1,DA_fov_crossTrack{:}); myFov_crossTrack{ir} = tem_fov_crossTrack(randOrder); %clear  tem_fov_crossTrack
        tem_fov_alongTrack = cat(1,DA_fov_alongTrack{:}); myFov_alongTrack{ir} = tem_fov_alongTrack(randOrder); %clear tem_fov_alongTrack

        tem_times = cat(1,DA_times{:}); myTimes{ir} = tem_times(randOrder); %clear tem_times
        tem_chNum = cat(1,DA_chNum{:}); myChNum{ir} = tem_chNum(randOrder); %clear myTb_chNum

        tem_obsErr = cat(1, DA_obsError{:}); myObsErr{ir} = tem_obsErr(randOrder); %clear tem_obsErr
 
        % deal with ROI
        DA_ROI_other = cell(size(Swath_used{iTb}));
        DA_ROI_hydro = cell(size(Swath_used{iTb}));
        for it = 1:length(Swath_used{iTb}) 
            obs_index_array = obs_index{it};
            obs_index_1d = obs_index_array(obs_index_array(:) == obs_index_array(:));
            DA_ROI_other{it} = ones(numel(obs_index_1d),1)*control.roi_oh{ir}(1);
            DA_ROI_hydro{it} = ones(numel(obs_index_1d),1)*control.roi_oh{ir}(2);
        end
        tem_ROI_other = cat(1,DA_ROI_other{:}); myRoi_otherVars{ir} = tem_ROI_other(randOrder); %clear tem_ROI_other
        tem_ROI_hydro = cat(1,DA_ROI_hydro{:}); myRoi_hydro{ir} = tem_ROI_hydro(randOrder); %clear tem_ROI_hydro
    end
    clear DA_lat DA_lon DA_sat_lat DA_sat_lon DA_sat_alt DA_sat_azimuth DA_scan_angle DA_zenith_angle DA_fov_crossTrack DA_fov_alongTrack DA_times DA_chNum DA_obsError DA_ROI_other DA_ROI_hydro
	

end

% Ideas behind the algorithm
%   If there is only one Tb file at a time
%  -------------------------------
% | Channel 1 ~ ROI combination 1 |  randomize    ---------
% |                               | -----------> |         | ----
% | Channel 2 ~ ROI combination 1 |               ---------      |             ---------
%  -------------------------------                               |concatenate |         | 
%                                                                |----------->|         |                                       
%                                                                |             ---------                                                               |
%  -------------------------------                               |
% | Channel 1 ~ ROI combination 2 |  randomize    ---------      |
% |                               | -----------> |         |-----
% | Channel 2 ~ ROI combination 2 |               ---------
%  -------------------------------


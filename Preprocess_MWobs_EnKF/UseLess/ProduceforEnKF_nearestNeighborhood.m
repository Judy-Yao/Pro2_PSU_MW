% ========================================================================================================================
% This function reads L1C data from MW obs files into the memory &
% generates a grid mesh centered at the location of the storm at the DA_time using the same setup used for WRF simulation (nx,ny,dx,dy) &
% separates the mesh into different parts for different frequencies of interest &
% for each grid point in a mesh part, selects a L1C MW obs and return its characteristics such as Tb value, location, angles, scan time...   
% ========================================================================================================================

function [sat_name,myLat,myLon,myTb,mySat_lat,mySat_lon,mySat_alt,myAzimuth,myScan_angle,myZenith,myFov_crossTrack,myFov_alongTrack,myTimes,myChNum,myRoi_hydro,myRoi_otherVars,myObsErr] = ProduceforEnKF(istorm,iTb,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,Tb_file,control) 

% ================================================================================================
%                                         Step 1
% ================================================================================================
% For each L1C MW observation file, loop through each channel/frequency of interest and
% read characteristics of sensor into memory
% ------------------------------------------------------------------------------------------------
    
    % Preallocating memory
    lat = cell(size(Swath_used{iTb}));
    lon = cell(size(Swath_used{iTb}));
    Tb = cell(size(Swath_used{iTb}));
    zenith = cell(size(Swath_used{iTb}));
    sat_lat = cell(size(Swath_used{iTb}));
    sat_lon = cell(size(Swath_used{iTb}));
    sat_alt = cell(size(Swath_used{iTb}));
    azimuth = cell(size(Swath_used{iTb}));
    outime = cell(size(Swath_used{iTb}));

    % Get platform and sensor names from the symbolically-linked Tb file name 
    [filepath,filename,filext] = fileparts(Tb_file);
    ss_info = split(filename,'.'); platform = ss_info(2); sensor = ss_info(3);

    % Modify satellite-and-sensor name so that it is consistent with what is used in the CRTM package 
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
    elseif (contains(platform,'F15'))
        sat_name = "ssmi_f15";
    elseif (contains(platform,'F16'))
        sat_name = "ssmis_f16";
    elseif (contains(platform,'F17'))
        sat_name = "ssmis_f17";
    elseif (contains(platform,'F18'))
        sat_name = "ssmis_f18";
    end

    % Read values from MW obs file
    for it = 1:length(Swath_used{iTb}) % it: item
        
		% HDF5 file
        if contains(filext,"HDF5")  
        	% **latitude**
        	lat{it} = h5read(Tb_file,Swath_used{iTb}(it) + '/Latitude'); % npixel, nscan

        	% **longitude**
        	lon{it} = h5read(Tb_file,Swath_used{iTb}(it) + '/Longitude'); % npixel, nscan

        	% **Tb**
	        Tb_allCh = h5read(Tb_file,Swath_used{iTb}(it) + '/Tc');
    	    Tb{it} = squeeze(Tb_allCh(ChIdx_all{iTb}(it),:,:)); % npixel, nscan

			% **zenith angle**: the angle of the satellite from the local zenith as seen at the pixel location on the earth
			%     satellite   local zenith
			%         \       /
			%          \zenith
			%           \~ /
			%            \/
			%            FOV
			% ----------------------- earth surface
            %
			% Note: Interpretation on nChUIA (number of unique incidence angles for this swath)
            % For a sensor, there might be n sets of scanning setup/incidence angles. 
            % For example, GMI sensor has two swaths, each of which uses different scanning setup; 
            % the dimension name for IncidenceAngle for swath 1 or 2 is (nscan1, npixel1, nchUIA1) or (nscan2, npixel2, nchUIA2)
            % However, it is possible that even for one swath there might be more than one set of scanning setup such as nchUIA1 might be 2.
            % In such a case, for example, there are 2 channels under one swath and each channel points at different set of geometry value.

        	incidenceAngles = h5read(Tb_file,Swath_used{iTb}(it) + '/incidenceAngle'); % nChUIA, npixel, nscan
        	iAIndices = h5read(Tb_file,Swath_used{iTb}(it) + '/incidenceAngleIndex'); % nchannel, nscan
        	if size(incidenceAngles,1) == 1 % only one scanning setup for the swath
            	zenith{it}(:,:) = squeeze(incidenceAngles(1,:,:));
        	else
            	num_Ch_perSW = size(iAIndices,1); % number of channels under this swath 
            	num_unique_perCh = zeros(1,num_Ch_perSW); 
            	for ich = 1:num_Ch_perSW
                	num_unique_perCh(ich) = length(unique(iAIndices(ich,:)));
            	end
            	if sum(num_unique_perCh) == num_Ch_perSW % i.e., all of scans of a channel only points at one set of geometry values
                	iAIndices = iAIndices(:,1); % For a channel, all scans point to the same set of zenith angles
                	zenith{it}(:,:) = squeeze(incidenceAngles(iAIndices(ChIdx_all{iTb}(it)),:,:)); % npixel, nscan
            	else % more complicated: some scans points to first set of geometry while others point to other sets of geometry
                	disp(['Error: current algorithm does not work!! Please modify it.']); 
            	end
        	end

        	% **latitude of satellite scan**
        	sat_lat{it} = h5read(Tb_file,Swath_used{iTb}(it)+ '/SCstatus/SClatitude'); % nscan

        	% **longitude of satellite scan**
        	sat_lon{it} = h5read(Tb_file,Swath_used{iTb}(it)+ '/SCstatus/SClongitude'); % nscan

        	% **altitude of satellite scan**
        	sat_alt{it} = h5read(Tb_file,Swath_used{iTb}(it)+'/SCstatus/SCaltitude'); % nscan (unit: km)

  	      	% **azimuth of FOV (used in CRTM)**: the angle subtended
    	 	% by the horizontal projection of a direct line from the satellite
	      	% to the FOV and the North-South axis of the FOV measured clockwise from North
    	  	% (0-> 360 degrees)
  	      	%                     North
  	      	%                      |AZ /
      	  	%                      |~ Satellite
          	%                      | /
          	%                      |/
          	% West -------------- FOV ------------- East
        	Npixel_perscan{it} = size(Tb{it},1); % number of pixels per scan
            azimuth{it} = geodetic2aer( repmat(sat_lat{it}',[Npixel_perscan{it} 1]),  repmat(sat_lon{it}',[Npixel_perscan{it} 1]), repmat(sat_alt{it}'*1000,[Npixel_perscan{it} 1]),...
                                        lat{it},                                       lon{it},                                       0, referenceEllipsoid('WGS 84')); % npixel, nscan
            % geodetic2aer: transforms the geodetic coordinates specified by lat, lon, and h to the local azimuth-elevation-range (AER) spherical coordinates specified by az, elev, and slantRange.
            % The local AET system should be specified by the FOV, consistent with CRTM definition.

        	% ** (Scan) Time** 
        	year_scans   = double(h5read(Tb_file, Swath_used{iTb}(it) + '/ScanTime/Year')); % nscan
        	month_scans  = double(h5read(Tb_file, Swath_used{iTb}(it) + '/ScanTime/Month')); % nscan
        	day_scans    = double(h5read(Tb_file, Swath_used{iTb}(it) + '/ScanTime/DayOfMonth')); % nscan
        	hour_scans   = double(h5read(Tb_file, Swath_used{iTb}(it) +  '/ScanTime/Hour')); % nscan
        	minute_scans = double(h5read(Tb_file, Swath_used{iTb}(it) +  '/ScanTime/Minute')); % nscan
        	second_scans = double(h5read(Tb_file, Swath_used{iTb}(it) +  '/ScanTime/Second')); % nscan
        
        	num_scans= size(Tb{it},2); % number of scans
        	for i_time = 1:num_scans
            	datetime_temp = datetime(year_scans(i_time), month_scans(i_time), day_scans(i_time), hour_scans(i_time), minute_scans(i_time), second_scans(i_time));
            	if ( second_scans(i_time) >= 30)
                	datetime_temp = datetime_temp + minutes(1);
            	end
            	outime{it}(i_time,1) = string(datestr(datetime_temp,'yyyymmddHHMM')); % nscan
        	end

		elseif contains(filext,"nc")
       
		    % **latitude**
            lat{it} = ncread(Tb_file, ['lat_' + Swath_used{iTb}(it)]); % npixel, nscan

            % **longitude**
			lon{it} = ncread(Tb_file, ['lon_' + Swath_used{iTb}(it)]); % npixel, nscan

            % **Tb**
			Tb{it} = ncread(Tb_file, ChName_all{iTb}(it)); % npixel, nscan

            % **zenith angle**: the angle of the satellite from the local zenith as seen at the pixel location on the earth
            zenith{it} = ncread(Tb_file, ['eia_' + Swath_used{iTb}(it)]); % npixel, nscan

            % **latitude of satellite scan**
			sat_lat{it} = ncread(Tb_file, ['spacecraft_lat_' + Swath_used{iTb}(it)]); % nscan

            % **longitude of satellite scan**
            sat_lon{it} = ncread(Tb_file, ['spacecraft_lon_' + Swath_used{iTb}(it)]); % nscan

            % **altitude of satellite scan**
			sat_alt{it} = ncread(Tb_file, ['spacecraft_alt_' + Swath_used{iTb}(it)]); % nscan

            % **azimuth of satellite scan (used in CRTM)**: 
            Npixel_perscan{it} = size(Tb{it},1); % number of pixels per scan
            azimuth{it} = geodetic2aer( repmat(sat_lat{it}',[Npixel_perscan{it} 1]),  repmat(sat_lon{it}',[Npixel_perscan{it} 1]), repmat(sat_alt{it}'*1000,[Npixel_perscan{it} 1]),...
                                        lat{it},                                       lon{it},                                       0, referenceEllipsoid('WGS 84')); % npixel, nscan

            % ** (Scan) Time** 
            scan_datetime_name = ['scan_datetime_' + Swath_used{iTb}(it)]; 
            scan_datetime = ncread(Tb_file,scan_datetime_name)'; % nscan, numchar

            year_scans = str2num(scan_datetime(:,1:4));
            month_scans = str2num(scan_datetime(:,6:7));
            day_scans = str2num(scan_datetime(:,9:10));
            hour_scans = str2num(scan_datetime(:,12:13));
            minute_scans = str2num(scan_datetime(:,15:16));        
			second_scans = str2num(scan_datetime(:,18:19));

            num_scans= size(Tb{it},2); % number of scans

            for i_time = 1:num_scans
                datetime_temp = datetime(year_scans(i_time), month_scans(i_time), day_scans(i_time), hour_scans(i_time), minute_scans(i_time), second_scans(i_time));
                if ( second_scans(i_time) >= 30) 
                    datetime_temp = datetime_temp + minutes(1);
                end
                outime{it}(i_time,1) = string(datestr(datetime_temp,'yyyymmddHHMM')); % nscan
            end
		
		end
    end

    % Special treatment to  8xGHz
    if (sensor == "AMSR2") | (sensor == "SSMI")
        [Swath_used,ChIdx_all,ChName_all,lat,lon,Tb,zenith,sat_lat,sat_lon,sat_alt,azimuth,outime] = Handle_8xGHz(iTb,sensor,Swath_used,ChIdx_all,ChName_all,DAtime_all,lat,lon,Tb,zenith,sat_lat,sat_lon,sat_alt,azimuth,outime,control);
    end


% ================================================================================================
%                                         Step 2
% ================================================================================================
%  Read domain of interest from geo_em.d03.nc generated by WPS geogrid.exe 
%  Select the raw obs for staggered grid point for EnKF assimilation (indicated by obs_index)
% -----------------------------------------------------------------------------------------------

    geo_file =  strcat(control.geogrid_dir,control.storm_phase{istorm},'/',DAtime_all{iTb(1)}(1),'/geo_em.',control.domain,'.nc');
    disp(strcat('Reading ', geo_file, '......'));
    xlat_m = ncread(geo_file,'XLAT_M');
    xlon_m = ncread(geo_file,'XLONG_M');
    min_xlat = min(xlat_m,[],'all')-1;
    max_xlat = max(xlat_m,[],'all')+1;
    min_xlon = min(xlon_m,[],'all')-1;
    max_xlon = max(xlon_m,[],'all')+1;
   
    % initialize an empty container for indices of obs to be selected for ROI plans
    obs_index = cell(length(control.roi_oh),length(Swath_used{iTb}));
    % locate the start point of model grid for each ROI plan
    filter_grid_step = control.filter_reso / control.dx;
    grid_start(2) = floor(filter_grid_step(2) / 2); % start point for ROI plan 1 
    grid_start(1) = grid_start(2) + .5*(2*filter_grid_step(2) - filter_grid_step(1)); % start point for ROI plan 2 
    
    t_Start_filter = tic;
    %parpool(48); % use 48 cores in a node
    % Note: for MW data, the (i,j) obs of different channels might have different location. Therefore, loop over channels is needed.
    % For each ROI plan, find the nearest obs for each WRF grid
    for it = 1:length(Swath_used{iTb})
        disp(["  "+ChName_all{iTb}(it)+'...']);

        % reduce workload by find the approximate area of interest
        lat_col = reshape(lat{it},[],1);
        lon_col = reshape(lon{it},[],1);
        lg_Areaidx_inAll = (lat_col >= min_xlat) & (lat_col <= max_xlat) & (lon_col >= min_xlon) & (lon_col <= max_xlon);
        lat_obs = lat_col(lg_Areaidx_inAll);
        lon_obs = lon_col(lg_Areaidx_inAll);
        Areaidx_inAll = find( lg_Areaidx_inAll );

        for iroi = 1:length(control.roi_oh)
            disp(['  Resolution to be filtered for this ROI is : ', num2str(control.filter_reso(iroi)), 'km ......']);
            if iroi == 1
                lat_obs_toPick = lat_obs; lon_obs_toPick = lon_obs;
            end
            [lat_obs_toPick, lon_obs_toPick, idx_PickObs] = PickRawforCRTM(lat_obs_toPick,lon_obs_toPick,xlat_m,xlon_m,grid_start(iroi),filter_grid_step(iroi),control);
            idx_PickObs_inAll = Areaidx_inAll( idx_PickObs );
            obs_index{iroi,it} =  idx_PickObs_inAll;
            [v, w] = unique(obs_index{iroi,it},'stable' );
            duplicate_idx = setdiff(1:numel(obs_index{iroi,it}),w );
            if sum(duplicate_idx) ~= 0
                warning(['  Number of repeated obs selected for this ROI selection: ', num2str(length(duplicate_idx))]);
            end 
        end
        % find the possible repeated selected obs for all ROI plans
        [Cdata] = intersect(obs_index{1,it},obs_index{2,it});
        if sum(Cdata) ~= 0
            disp(['Number of repeated obs across all ROI selections: ', num2str(length(Cdata))]);
        end
    end
    delete(gcp('nocreate')); % shut down current parallel pool
    t_End_filter = toc(t_Start_filter);
    disp(['Picking the nearest obs for model grids for all ROI plans over all channels takes ', num2str(t_End_filter), ' seconds.']);
    % Original grid points:         Filtered grid points:
    %       * * * * *                *   *   * 
    %       * * * * *                          
    %       * * * * *                *   *   * 
    %       * * * * *                          

% ================================================================================================
%                                         Step 3
% ================================================================================================ 
%  Select MW records for EnKF assimilation
% -------------------------------------------------------------------------------------------------
    
    % Preallocating memory
    DA_lat = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_lat_perROI = cell(length(control.roi_oh));
    DA_lon = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_lon_perROI = cell(length(control.roi_oh));
    DA_Tb = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_Tb_perROI = cell(length(control.roi_oh)); 

    DA_sat_lat = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_satLat_perROI = cell(length(control.roi_oh));
    DA_sat_lon = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_satLon_perROI = cell(length(control.roi_oh));
    DA_sat_alt = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_satAlt_perROI = cell(length(control.roi_oh));
    DA_azimuth = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_azimuth_perROI = cell(length(control.roi_oh));
    DA_scan_angle = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_scan_perROI = cell(length(control.roi_oh));
    DA_zenith = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_zenith_perROI = cell(length(control.roi_oh));
    DA_fov_crossTrack = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_fovCross_perROI = cell(length(control.roi_oh));
    DA_fov_alongTrack = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_fovAlong_perROI = cell(length(control.roi_oh));

    DA_times = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_times_perROI = cell(length(control.roi_oh));
    DA_chNum = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_chNum_perROI = cell(length(control.roi_oh));
    DA_obsError = cell(length(control.roi_oh),length(Swath_used{iTb})); DA_obsError_perROI = cell(length(control.roi_oh));
    % Assign values
	for iroi = 1:length(control.roi_oh)
		for it = 1:length(Swath_used{iTb})
            obs_index_array = obs_index{iroi,it};
			obs_index_1d = obs_index_array(obs_index_array(:) == obs_index_array(:)); % get rid of obs_index with value NaN
            % **Tb,lat,lon**
			allTb = reshape(Tb{it},[],1); all_lat = reshape(lat{it},[],1); all_lon = reshape(lon{it},[],1);
			DA_Tb{iroi,it}          = allTb(obs_index_1d);
			DA_lat{iroi,it}         = all_lat(obs_index_1d);
			DA_lon{iroi,it}         = all_lon(obs_index_1d);
			% **sat_lat,sat_lon,sat_alt,zenith angle**
			[scan_position, scan_num] = ind2sub(size(Tb{it}),obs_index_1d); % size(Tb{it}): npixel, nscan; output from ind2sub: row, column
			DA_sat_lat{iroi,it} = sat_lat{it}(scan_num);
			DA_sat_lon{iroi,it} = sat_lon{it}(scan_num);
			DA_sat_alt{iroi,it} = sat_alt{it}(scan_num);
			all_azimuth = reshape(azimuth{it},[],1);
			DA_azimuth{iroi,it} = all_azimuth(obs_index_1d);
			all_zenith = reshape(zenith{it},[],1);
			DA_zenith{iroi,it} = all_zenith(obs_index_1d);
			% **scan angle** 
			[scantype,ch_num,fov_alongTrack,fov_crossTrack,max_scan_angle,scan_angles] = SensorInfo_read(sensor,ChName_all{iTb}{it});
            DA_scan_angle{iroi,it} = scan_angles(scan_position);
            % **cross-track and along-track FOV**
			[DA_fov_crossTrack{iroi,it}, DA_fov_alongTrack{iroi,it}] = Get_pixel_resolution(scantype,ch_num,fov_alongTrack,fov_crossTrack,DA_lat{iroi,it},DA_lon{iroi,it},DA_zenith{iroi,it},DA_sat_lat{iroi,it},DA_sat_lon{iroi,it},DA_sat_alt{iroi,it});
			% **scan time**
			for my_scan_num_idx = 1:length(scan_num)
				DA_times{iroi,it}(my_scan_num_idx,1) = outime{it}(scan_num(my_scan_num_idx));
			end
			% **Channel number, obs error**
			DA_chNum{iroi,it} = ones(numel(obs_index_1d),1,'int64')*ch_num;
			DA_obsError{iroi,it} = ones(numel(obs_index_1d),1)*control.obsError(it);
		end
        % For each ROI plan, combine low-frequency and high-frequency data into one column if there are at least 2 channels 
		DA_lat_perROI{iroi} = cat(1,DA_lat{iroi,:}); DA_lon_perROI{iroi} = cat(1,DA_lon{iroi,:}); DA_Tb_perROI{iroi} = cat(1,DA_Tb{iroi,:});
        DA_satLat_perROI{iroi} = cat(1,DA_sat_lat{iroi,:}); DA_satLon_perROI{iroi} = cat(1,DA_sat_lon{iroi,:}); DA_satAlt_perROI{iroi} = cat(1,DA_sat_alt{iroi,:});

        DA_azimuth_perROI{iroi} = cat(1,DA_azimuth{iroi,:}); DA_scan_perROI{iroi} = cat(1,DA_scan_angle{iroi,:}); DA_zenith_perROI{iroi} = cat(1,DA_zenith{iroi,:});    
		DA_fovCross_perROI{iroi} = cat(1,DA_fov_crossTrack{iroi,:}); DA_fovAlong_perROI{iroi} = cat(1,DA_fov_alongTrack{iroi,:}); 
		DA_times_perROI{iroi} = cat(1,DA_times{iroi,:});  DA_chNum_perROI{iroi} = cat(1,DA_chNum{iroi,:});  DA_obsError_perROI{iroi} = cat(1,DA_obsError{iroi,:});
	end
    clear Tb lat lon sat_lat sat_lon sat_alt azimuth zenith scan_angles outime 
    clear allTb all_lat all_lon all_azimuth all_zenith 
	clear DA_lat DA_lon DA_Tb DA_sat_lat DA_sat_lon DA_sat_alt DA_azimuth DA_azimuth DA_zenith DA_fov_crossTrack DA_fov_alongTrack 
	clear DA_times DA_chNum DA_obsError DA_ROI_hydro DA_ROI_other 

% ================================================================================================
%                                         Step 4
% ================================================================================================
% Randomize selected MW records
% ------------------------------------------------------------------------------------------------

    myLat = cell(size(control.roi_oh));
    myLon = cell(size(control.roi_oh));
    myTb = cell(size(control.roi_oh));

    mySat_lat = cell(size(control.roi_oh));
    mySat_lon = cell(size(control.roi_oh));
    mySat_alt = cell(size(control.roi_oh));
    myAzimuth = cell(size(control.roi_oh)); 
    myScan_angle = cell(size(control.roi_oh));
    myZenith = cell(size(control.roi_oh));
    myFov_crossTrack = cell(size(control.roi_oh));
    myFov_alongTrack = cell(size(control.roi_oh));

    myTimes = cell(size(control.roi_oh));
    myChNum = cell(size(control.roi_oh));
    myRoi_hydro = cell(size(control.roi_oh));
    myRoi_otherVars = cell(size(control.roi_oh));
    myObsErr = cell(size(control.roi_oh));

    for ir = 1:length(control.roi_oh) % ROI first references [200,0] second [60,60]
        % randomize MW records for this ROI
        if (control.random)
            randOrder = randperm(length(DA_Tb_perROI{ir}));
        else
            randOrder = 1:length(DA_Tb_perROI{ir});
        end
            
        tem_lat = DA_lat_perROI{ir}; myLat{ir} = tem_lat(randOrder); %clear tem_lat
        tem_lon = DA_lon_perROI{ir}; myLon{ir} = tem_lon(randOrder);
        tem_Tb = DA_Tb_perROI{ir};   myTb{ir} = tem_Tb(randOrder);
           
        tem_sat_lat = DA_satLat_perROI{ir}; mySat_lat{ir} = tem_sat_lat(randOrder);
        tem_sat_lon = DA_satLon_perROI{ir}; mySat_lon{ir} = tem_sat_lon(randOrder);
        tem_sat_alt = DA_satAlt_perROI{ir}; mySat_alt{ir} = tem_sat_alt(randOrder);
        tem_azimuth = DA_azimuth_perROI{ir}; myAzimuth{ir} = tem_azimuth(randOrder);
        tem_scan_angle = DA_scan_perROI{ir}; myScan_angle{ir} = tem_scan_angle(randOrder);
        tem_zenith = DA_zenith_perROI{ir}; myZenith{ir} = tem_zenith(randOrder);
        tem_fov_crossTrack = DA_fovCross_perROI{ir}; myFov_crossTrack{ir} = tem_fov_crossTrack(randOrder);
        tem_fov_alongTrack = DA_fovAlong_perROI{ir}; myFov_alongTrack{ir} = tem_fov_alongTrack(randOrder);
        tem_times = DA_times_perROI{ir}; myTimes{ir} = tem_times(randOrder);
        tem_chNum = DA_chNum_perROI{ir}; myChNum{ir} = tem_chNum(randOrder);

        tem_obsErr = DA_obsError_perROI{ir}; myObsErr{ir} = tem_obsErr(randOrder);

        % deal with ROI
        DA_ROI_other = cell(size(Swath_used{iTb}));
        DA_ROI_hydro = cell(size(Swath_used{iTb}));
        for it = 1:length(Swath_used{iTb}) 
            obs_index_array = obs_index{ir,it};
            obs_index_1d = obs_index_array(obs_index_array(:) == obs_index_array(:));
            DA_ROI_other{it} = ones(numel(obs_index_1d),1)*control.roi_oh{ir}(1);
            DA_ROI_hydro{it} = ones(numel(obs_index_1d),1)*control.roi_oh{ir}(2);
        end
        tem_ROI_other = cat(1,DA_ROI_other{:}); myRoi_otherVars{ir} = tem_ROI_other(randOrder); %clear tem_ROI_other
        tem_ROI_hydro = cat(1,DA_ROI_hydro{:}); myRoi_hydro{ir} = tem_ROI_hydro(randOrder); %clear tem_ROI_hydro 
    end
	

    % -------------------------------------------------------------------------
    % Function to find the nearest obs for filtered WRF grids    
    % -------------------------------------------------------------------------
    function [ obs_lat, obs_lon, idx_getObs ] = PickRawforCRTM(obs_lat,obs_lon,m_xlat,m_xlon,start_grid,step_grid,control)

        % Filter lat and lon of WRF domain for each ROI
        idx_step = start_grid:step_grid:size(m_xlat,1); % xlat_m/xlon_m has same dimension value along x and y axis 
        xlat = m_xlat(idx_step, idx_step);
        xlon = m_xlon(idx_step, idx_step);
        xlat_col = xlat(:);
        xlon_col = xlon(:);

        % For each WRF grid of interest, select the nearest obs 
        idx_getObs = nan(length(xlon_col),1);

        for id =1:length(xlon_col)
            dis_ig = distance(xlat_col(id),xlon_col(id),obs_lat,obs_lon);
            idx = find( dis_ig == min(dis_ig) );
            % Get the obs index 
            idx_getObs(id,1) = idx;
            % Mark the selected obs and gurantee it will be never used again
            obs_lat(idx) = nan; obs_lon(idx) = nan;
        end
    end

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


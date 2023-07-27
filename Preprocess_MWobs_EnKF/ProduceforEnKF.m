% ========================================================================================================================
% This function reads L1C data from MW obs files into the memory &
% generates a grid mesh centered at the location of the storm at the DA_time using the same setup used for WRF simulation (nx,ny,dx,dy) &
% separates the mesh into different parts for different frequencies of interest &
% for each grid point in a mesh part, selects a L1C MW obs and return its characteristics such as Tb value, location, angles, scan time...   
% ========================================================================================================================

function [sat_name,myLat,myLon,myTb,mySat_lat,mySat_lon,mySat_alt,myAzimuth,myScan_angle,myZenith,myFov_crossTrack,myFov_alongTrack,myTimes,myChNum,myRoi_hydro,myRoi_otherVars,myObsErr] = ProduceforEnKF(iTb,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,Tb_file,control) 

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

  	      	% **azimuth of satellite scan (used in CRTM)**: the angle subtended
    	 	% by the horizontal projection of a direct line from the satellite
	      	% to the FOV and the North-South axis measured clockwise from North
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
             

        	% ** (Scan) Time** 
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

            % **azimuth of satellite scan (used in CRTM)**: the angle subtended
            Npixel_perscan{it} = size(Tb{it},1); % number of pixels per scan
            azimuth{it} = geodetic2aer(lat{it},                                       lon{it},                                       0, ...
                                             repmat(sat_lat{it}',[Npixel_perscan{it} 1]), repmat(sat_lon{it}',[Npixel_perscan{it} 1]), repmat(sat_alt{it}'*1000,[Npixel_perscan{it} 1]), ...
                                             referenceEllipsoid('WGS 84')); % npixel, nscan

            % ** (Scan) Time** 
            scan_datetime_name = ['scan_datetime_' + Swath_used{iTb}(it)];                                                                 
            scan_datetime = ncread(Tb_file,scan_datetime_name)'; % nscan, numchar

            year = str2num(scan_datetime(:,1:4));
            month = str2num(scan_datetime(:,6:7));
            day = str2num(scan_datetime(:,9:10));
            hour = str2num(scan_datetime(:,12:13));
            minute = str2num(scan_datetime(:,15:16));        
			second = str2num(scan_datetime(:,18:19));

            num_scans= size(Tb{it},2); % number of scans

            for i_time = 1:num_scans
                datetime_temp = datetime(year(i_time), month(i_time), day(i_time), hour(i_time), minute(i_time), second(i_time));
                if ( second(i_time) >= 30)
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

    geo_file =  [control.geogrid_dir,control.storm_phase{istorm},'/',DAtime,'/geo_em.',control.domain,'.nc'];
    disp(['Reading ', geo_file, '......']);
    xlat_m = ncread(geo_file,'XLAT_M');
    xlon_m = ncread(geo_file,'XLONG_M');
    min_xlat = min(xlat_m,[],'all')-0.5;
    max_xlat = max(xlat_m,[],'all')+0.5;
    min_xlon = min(xlon_m,[],'all')-0.5;
    max_xlon = max(xlon_m,[],'all')+0.5;
    %% Old method: when there is no geo_em.d03.nc 
    %% Prepare area: best-track location followed; a slightly larger area that WRF simulation
    %nx = control.nx*control.domain_buffer; % zoom out
    %ny = control.ny*control.domain_buffer; % zoom out
    %% Below algorithm works if only for all frequencies of interest, the DA_time are the same
	%min_XLAT = loc_DAtime_all{iTb}(1) - (ny/2*control.dx)/(cos(loc_DAtime_all{iTb}(1)*(pi/180))*111);
    %max_XLAT = loc_DAtime_all{iTb}(1) + (ny/2*control.dx)/(cos(loc_DAtime_all{iTb}(1)*(pi/180))*111);
    %min_XLONG = loc_DAtime_all{iTb}(2) - (nx/2*control.dx)/111;
    %max_XLONG = loc_DAtime_all{iTb}(2) + (nx/2*control.dx)/111;
    %disp(['      min of xlong: ',num2str(min_XLONG), ', max of xlong: ',num2str(max_XLONG)]);
    %disp(['      min of xlat: ',num2str(min_XLAT), ', max of xlat: ',num2str(max_XLAT)]);
    %latitudes  = linspace(min_XLAT,max_XLAT,ny);
    %longitudes = linspace(min_XLONG,max_XLONG,nx);
    %[XLAT, XLONG] = meshgrid(latitudes,longitudes);
    
    % Note: for two ROI plans, Tb at/near different grid points are appreciated because Tb at the same point should not be picked for both large ROI plan and small ROI plan.
    % Therefore it is important to staggeredly separate gird points for small ROI plan from for large ROI plan.
    slots_x = cell(size(control.roi_oh));
    slots_y = cell(size(control.roi_oh));
    obs_index = cell(length(control.roi_oh),length(Swath_used{iTb}));
    
	% Identify the start point for each ROI plan
    filter_grid_step = control.filter_reso / control.dx;
    grid_start(2) = floor(filter_grid_step(2) / 2); % start point for one ROI plan  
    grid_start(1) = grid_start(2) + .5*(2*filter_grid_step(2) - filter_grid_step(1)); % start point for the other ROI plan 

    % For each ROI plan, find the nearest obs for each WRF grid
    % Note: for GOES_IR data, the (i,j) obs of 16 channels has the same location. Therefore, loop over channels is not needed.
	for iroi = 1:length(control.roi_oh)
		slots_x{iroi} = grid_start(iroi):filter_grid_step(iroi):nx;
        slots_y{iroi} = grid_start(iroi):filter_grid_step(iroi):ny;
        % Original grid points:         Filtered grid points:
        %       * * * * *                *   *   * 
        %       * * * * *                          
        %       * * * * *                *   *   * 
        %       * * * * *                          
        %       * * * * *                *   *   * 
        %PickRawforCRTM(lat,lon,Tb,min_XLONG,max_XLONG,min_XLAT,max_XLAT,latitudes,longitudes,slots_x{iroi},slots_y{iroi},obs_indexqcontrol);
		for it = 1:length(Swath_used{iTb})
			obs_index{iroi,it} = PickRawforCRTM(lat{it},lon{it},Tb{it},min_XLONG,max_XLONG,min_XLAT,max_XLAT,latitudes,longitudes,slots_x{iroi},slots_y{iroi},control);
            % Note: Each value in the obs_index points at a location in the vectorized array lat_col/lon_col [lat_col = reshape(lat_raw,[],1)]
		end
	end


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
			%DA_scan_angle{iroi,it} = scan_angles(scan_position)';
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
        randOrder = randperm(length(DA_Tb_perROI{ir}));
        
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


% =============================================================================================================================
% This is the main script for the microwave-observation-preprocessing system (MOPS), showing the skeleton of the system.
% =============================================================================================================================

% -------------- Configuration: Set up control variables ----------------

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  WARNING  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Please carefully evaluate each item and it is the user's responsibility to adjust the value to their needs.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

control = struct;
% ----Path
control.obs_dir = '../../raw_Obs/Microwave/'; % directory into which MW L1C raw observations are downloaded using the script Download_MWobs/get_MWobs.sh
control.obs_collect_dir = '../../raw_Obs/Collected_MW/'; % directory into which a subset of MW L1C MW raw files are collected/linked which are needed in this study 
control.bestrack_dir = '../../raw_Obs/Bestrack/'; % directory where best-track files are
control.output_dir = '../../toEnKFobs/MW/'; % directory into which the microwave-observation-preprocessing system (MOPS) outputs

% ---Storm information
control.storm_phase = {'MariaRI',}; % !!! It is recommended to process ONE storm at a time in spite of MOPS's ability to process as many storms as possible.
%control.storm_phase = ["Irma2ndRI",'JoseRI','MariaRI'};
control.period = {{'201709160000','201709180000'},}; % Date range of case study (yyyymmddHHMM)
%control.period = {{'201709030600','201709050600'},{'201709050600','201709070600'},{'201709160000','201709180000'}}; %YYYYMMDDHHmm

% --- WRF simulation setup
control.nx = 297; % number of grid points along X direction
control.ny = 297; % number of grid points along Y direction
control.dx = 3; % WRF resolution: 3 km

% ---Satellite informaiton
control.sensor = {'AMSR2','ATMS','GMI','MHS','SAPHIR','SSMI','SSMIS'}; % sensor name
control.platform = {{'GCOMW1'}, {'NPP'}, {'GPM'}, {'METOPA','METOPB','NOAA18','NOAA19'}, {'MT1'}, {'F15'}, {'F16','F17','F18'}}; % platform name (one sensor corresponds to at least one platform)
control.favFreq = {{'18.7GHzV-Pol','89GHzV-PolA-Scan','89GHzV-PolB-Scan'},{'183.31+-7GHzQH-Pol'},{'18.7GHzV-Pol','183.31+/-7GHzV-Pol'},{'190.31GHzV-Pol'},{'183.31+/-6.8GHz'},{'fcdr_tb19v','fcdr_tb85v'},{'19.35GHzV-Pol','183.31+/-6.6GHzH-Pol'}}; % frequencies of interest to this study (one sensor corresponds to at least one frequency/channel)
control.comnine_AMSR89GHz = true; % if combines two 89GHz channels on AMSR2 (89GHzV-PolA-Scan and 89GHzV-PolB-Scan)
%------------ NO NEED TO CHANGE --------------------------
control.use8xGHz = false; % if frequency 85GHz (from AMSR2) or 89GHz (from SSMI) is used. WARNING: the default value is false and this value will be adjusted automatically later in the sytem based on different conditions.
%---------------------------------------------------------

% --- Other parameters
control.domain_buffer = 1.5; % scaling factor
control.search_buffer = 0.1; % degrees: lat/lon % If the value is too large, the same observation might be picked for both ROI plans
control.filter_reso = [36;24]; % filter resolution for ROI 
% !!! Trick: Decreases the resolution if you'd like to quickly find out if MOPS is working by executing the system.
control.roi_oh = {[200,0]; [60,60]}; % ROI plans [other variables, hydrometeors]
control.obsError = [3;3];



% ===================================================== Begin the System =====================================================

date_runMOPS = date;
disp(['Starting running microwave-observation-preprocessing system (MOPS) on ', date_runMOPS]);

% ---------- Loop through each storm object -------------------
for istorm = 1:length(control.storm_phase)

    disp(['Storm (and stage): ', control.storm_phase{istorm}]);    

	% ============================================================================================================
	% Find times and center locations of the storm in the Best track file within the date range of the case study
	% (These data are only available on 00, 06, 12, 18 UTC)
	% ============================================================================================================
	
	bestrack_str = Bestrack_read(istorm, control); % (cell: {time, lat, lon})

    % ============================================================================================================
    % Collect useful MW Obs files of all sensors in all platforms into a directory every hour
	% ============================================================================================================
    
	disp('Collecting useful MW obs files for this study......');
	[Tbfile_names,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,overpass_t,singlepass_t] = Collect_MW_useful(istorm, bestrack_str, control); 
    % Tbfile_names: (string); ChIdx_all: (single); loc_DAtime_all: (cell: {double vector})
    % Swath_used, ChName_all, DAtime_all: (cell: {string,...})
    % overpass_t,singlepass_t: (string)

    % ============================================================================================================
    % Make output directory and output hourly-interpolated best-track data
    % ============================================================================================================

    % --- Make subdirectory for output
	if ~exist([control.output_dir,control.storm_phase{istorm}],'dir')
		[~, msg, ~] = mkdir(control.output_dir,control.storm_phase{istorm});
        if isempty(msg)
            disp(['Successfully created a subdirectory in ',control.output_dir,' for ',control.storm_phase{istorm}]);
        else
            error('Error: ',msg);
        end
	end

    disp('---------------------------------------------------------------');	
	% --- Output hourly best-track location and time
 	filename = strcat(control.output_dir,control.storm_phase{istorm},'/bestrack_perHour');	
	disp("Output hourly best-track location and time: " + filename);
	formatSpec = '%12s%12.3f%12.3f\n';
	fileID = fopen(filename,'w');
	for itime = 1:length(DAtime_all)
		fprintf(fileID, formatSpec, ...
                DAtime_all{itime}(1), loc_DAtime_all{1,itime}(1),loc_DAtime_all{1,itime}(2)); 
	end
    fclose(fileID);

    % ============================================================================================================
    % Output SO observation file under two situations: overpass or singlepass
	% ============================================================================================================
	% Singlepass: at one DA time, only one sensor provides MW obs of the storm
	% Overpass: ~, more than one sensor provides MW obs of the storm
	% ------------------------------------------------------------------------------------------------------------

    % --- Loop through each useful Tb file via a symbolic link ---
    Tb_dir = [control.obs_collect_dir,control.storm_phase{istorm},'/*'];
    Tb_files = strsplit(ls(Tb_dir));
    Tb_files = Tb_files(~cellfun('isempty',Tb_files)); % get rid of annoying empty cell

    disp('---------------------------------------------------------------');

    % --- Output singlepass ---
    disp('Handling single-pass Tb files......');
    for is = 1:length(singlepass_t)
        for iTb = 1:length(Tb_files)
           Tb_file = Tb_files{iTb};
            [filepath,filename,filext] = fileparts(Tb_file);            
            if contains(filename,singlepass_t(is))
                idx_collectedTb = find([filename,filext] == Tbfile_names);
                Singlepass_write(idx_collectedTb,istorm,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,Tb_file,control);
            else
                continue;
            end
        end
    end
    % Note: The order of collecting Tb files is different from the order of listing collected Tb files 
    %		- idx_collectedTb records the order of Tbs collected from different sensors & platforms in module Collect_MW_useful.m 
	%         (recall the process how useful MW obs and its parameters were collected)
    %		- iTb indicates the order of collected Tb files with ls command in a directory

    % Note: in a single-pass senario, if an AMSR2 or a SSMI Tb file exists, the low frequency (~19GHz) and the high frequency (85GHz) will be definitely used.

    disp('---------------------------------------------------------------');

	% --- Output overpass ---
    disp('Handling over-pass Tb files......');
	for io = 1:length(overpass_t) % loop through each DA time where overpass happens
		file_overpass = []; % (string)
		order_overpass = []; % (integer)
		sensor_overpass = []; % (string)
		idx_usedTb = []; % (integer)
        % for a specific time, find the overpass files     
		order_clt_io = 0;    
		for iTb = 1:length(Tb_files) % loop through all of collected Tb files
			Tb_file = Tb_files{iTb};
            [filepath,filename,filext] = fileparts(Tb_file);
            ss_info = split(filename,'.');
            sensor = ss_info{3};
			% gather names of overpass Tb files
			if contains(filename,overpass_t(io))
				order_clt_io = order_clt_io + 1; 
                order_overpass = [order_overpass, order_clt_io]; % Record the identifying order of Tb files that are overpass (e.g., 1,2,3)
                file_overpass = [file_overpass,string(Tb_file)]; % Record overpass Tb files in the order of being identified
                sensor_overpass = [sensor_overpass,string(sensor)];
                idx_collectedTb = find([filename,filext] == Tbfile_names);
				idx_usedTb(end+1) = idx_collectedTb; % Record the location of the overpass file in the collected information list
			else
				continue;
			end
		end
		
        disp('Processing over-pass Tb files: ');
		for iot = 1:length(sensor_overpass)
			disp(["  over-pass file: " + file_overpass(iot)]);
		end

		% --------------------------------- Special Treatment to 8x GHz -----------------------------------------------
        % 8x GHz will only be used if there is no 183 GHz 
		% ------------------------------------------------------------------------------------------------------------
		num_8xGHz_sensors = sum(("AMSR2" == sensor_overpass) | ("SSMI" == sensor_overpass));  % 89GHz from AMSR2 and 85GHz from SSMI

		if (num_8xGHz_sensors == 0) 
			% 8xGHz does not exist.
			disp('  None of microwave observation is from 8x (85 or 89) GHz.');
			Overpass_write(idx_usedTb,istorm,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,file_overpass,control);
		else
			% 8xGHz exists, consider if 183GHz exists.
			num_183GHz = 0;
			for iot = 1:length(sensor_overpass)
				if (sensor_overpass(iot) == "AMSR2")
					continue;
				elseif (sensor_overpass(iot) == "SSMI")
					continue;
				else
					ChName_other = ChName_all{idx_usedTb(iot)};	
					if contains(ChName_other,'183') | contains(ChName_other,'190')
						num_183GHz = num_183GHz + 1;
					end
				end
			end
			
			if num_183GHz == 0
                % 8xGHz exists, 183GHz does not exist
				control.use8xGHz = true;
				disp("  ~ 183 GHz from other sensors does not exist. Use 8x (85 or 89) GHz instead!");
                Overpass_write(idx_usedTb,istorm,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,file_overpass,control);	
			else
                % 8xGHz exists and 183GHz exists too
				control.use8xGHz = false;
				control.comnine_AMSR89GHz = false;
				disp("  ~ 183 GHz from other sensors exists. Only low frequency of AMSR2 or/and SSMI is used!");		
				Overpass_write(idx_usedTb,istorm,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,file_overpass,control);
			end

		end

    end % end loop for io = 1:length(overpass)


end % end loop for istorm = 1:length(control.storm_phase)

% ===================================================== End the System =====================================================

	

% ---------- Background Knowledge about L1C HDF5 file----------------
% Take GPM GMI as an example.
% -- Swath
% Swath S1 has nine channels (10V, 10H, 19V,19H, 23V, 37V, 37H, 89V, and 89H).
% Swath S2 has four channels (166V, 166H, 183+/-3V, and 183+/-7V).
% Relation between the swaths: Swath S2 has the same number of scans and the same number of
% pixels as swath S1. 
% -- Scan
% Each S1 scan contains nine channels sampled 221 times thus 221 pixels along the scan. 
% Each S2 scan contains four channels sampled 221 times thus 221 pixels along the scan.
% -- Granule
% A granule is a collection of scans. In this case, the size of a granule is one
% orbit, which begins and ends at the southernmost point.
% P.S. For AIRS data, a day's worth of data is divided into 240 granules,
% each of 6 min durations. Each granule consists of 135 scan lines
% containing 90 footprints/pixels per scan line; thus there are a total of
% 135*90 = 12,150 footprints per granule.
% -- Important content inside of 1C.GPM.GMI...HDF5 file
% --S1
% -latitude/longitude/Quality: per sampled footprint/pixel along per scan line.
% -SCstatus
% SCaltitude/SClatitude/SClongitude: per scan line.
% -ScanTime
% Anyvariable: per scan line.
% -Tc: per channel, per pixel along per scan line.

% There are two ways to represent text in MATLAB: Characters V.S. Strings
% Character: Single quoted. Eg. str = ['123','3456'] --> str = '1233456'
% String: Double quoted. Eg. str = ["123","456"] --> str = 1?~W2 string array "123"    "3456"



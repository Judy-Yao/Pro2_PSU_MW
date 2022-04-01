% The main script consists of the skeleton of the Process_MWobs_EnKF algorithm


% -------------- Set up control variables ----------------------
control = struct;
% ----Path
control.obs_dir = '../../Obs/Microwave/'; % directory where L1C raw observations are
control.obs_collect_dir = '../../Obs/Collected_MW/'; % directory into which L1C obs files are collected/linked 
control.bestrack_dir = '../../Obs/Bestrack/'; % directory where best-track files are
control.output_dir = '../../toEnKFobs/'; % directory where this algorithm outputs
% ---Storm information
control.storm_phase = {'Irma2ndRI',};  
%control.storm_phase = ["Irma2ndRI",'JoseRI','MariaRI'};
control.period = {{'201709030600','201709050600'},};
%control.period = {{'201709030600','201709050600'},{'201709050600','201709070600'},{'201709160000','201709180000'}}; %YYYYMMDDHHmm
% ---Satellite informaiton
%control.sensor = {'AMSR2',};
control.sensor = {'AMSR2','ATMS','GMI','MHS','SAPHIR','SSMIS'};
%control.platform = {{'GCOMW1'},};
control.platform = {{'GCOMW1'}, {'NPP'}, {'GPM'}, {'METOPA','METOPB','NOAA18','NOAA19'}, {'MT1'}, {'F16','F17','F18'}};

control.favFreq = {'18.7GHzV-Pol','19.35GHzV-Pol','89GHzV-PolA-Scan','89GHzV-PolB-Scan','183.31+/-6.6GHzH-Pol','183.31+/-6.8GHz','183.31+/-7GHzV-Pol','183.31+-7GHzH-Pol','190.31GHzV-Pol'}; % favorite frequencies
%control.favCh_sup = {'19.35GHzV-Pol','89GHzV-PolA-Scan','89GHzV-PolB-Scan','183.31+/-6.6GHzH-Pol','183.31+/-6.8GHz','190.31GHzV-Pol'};
% --- Dealing with 89GHzV-PolA-Scan and 89GHzV-PolB-Scan on AMSR2
control.comnine_AMSR89GHz = true;
control.NOTuse_AMSR89GHz = false;
% --- WRF setup
control.nx = 297; % number of grid points along X direction
control.ny = 297; % number of grid points along Y direction
control.dx = 3; % WRF resolution: 3 km
% --- Other
control.domain_buffer = 1.5; % scaling factor
control.search_buffer = 0.2; % degrees: lat/lon
control.filter_ratio = [36,24]; % speed the test
%control.filter_ratio = [6,6];
control.roi_oh = {[200,0], [60,60]}; % roi [other variables, hydrometeors]
control.obsError = [3 3];

% ---------- Loop through each storm object -------------------
for istorm = 1:length(control.storm_phase)

    % --- Find times and center locations of the storm in the Best track file (only available on 00, 06, 12, 18 UTC)
	bestrack_str = Bestrack_read(istorm, control); % (cell)

    % --- Collect useful MW Obs files of all sensors among all platforms into a directory
    disp('Collecting useful MW obs files for this study......');
	[Tbfile_names,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,overpass,singlepass] = Collect_MW_useful(istorm, bestrack_str, control); % ps: per swath
	
    % --- Loop through each useful Tb file via a symbolic link
	Tb_dir = [control.obs_collect_dir,control.storm_phase{istorm},'/*'];
	Tb_files = strsplit(ls(Tb_dir));
    Tb_files = Tb_files(cellfun(@isempty, Tb_files) == 0); % Get rid of the annyoing empty cell
	
    % --- Output file under two situations: overpass or single-pass
	if ~exist([control.obs_collect_dir,control.storm_phase{istorm}],'dir')
		[~, msg, ~] = mkdir(control.output_dir,control.storm_phase{istorm});
        if isempty(msg)
            disp(['Successfully created a subdirectory in ',control.ouput_dir,' for ',control.storm_phase{istorm}]);
        else
            error('Error: ',msg);
        end
	end

    % - Output single-pass
    disp('Handling single-pass Tb files......');
    for is = 1:length(singlepass)
        for iTb = 1:length(Tb_files)
           Tb_file = Tb_files{iTb};
            [filepath,filename,filext] = fileparts(Tb_file);            
            if contains(filename,singlepass(is))
                idx_collectedTb = find([filename,filext] == Tbfile_names)
                Singlepass_write(idx_collectedTb,istorm,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,Tb_file,control);
            else
                continue;
            end
        end

    end
    % Note: The order of collecting Tb files is different from the order of listing collected Tb files &
    % idx_collectedTb records the order of Tbs collected from different sensors & platforms in module Collect_MW_useful.m &
    % iTb indicates the order of collected Tb files with ls command in a directory

    % Note: in single-pass senario, if an AMSR2 Tb file exists, it will be used anyway

	% - Output overpass
    disp('Handling over-pass Tb files......');
	for io = 1:length(overpass)
		file_overpass = []; % (strings)
		order_overpass = [];
		sensor_overpass = [];
		idx_usedTb = []; % (integer)
        % for a specific time, find the overpass files         
		for iTb = 1:length(Tb_files)
			Tb_file = Tb_files{iTb};
            [filepath,filename,filext] = fileparts(Tb_file);
            ss_info = split(filename,'.');
            sensor = ss_info(3);
			% gather names of overpass Tb files
			if contains(filename,overpass(io))
                order_overpass = [order_overpass, io];
                file_overpass = [file_overpass,string(Tb_file)];
                sensor_overpass = [sensor_overpass,string(sensor)];
                idx_collectedTb = find([filename,filext] == Tbfile_names);
				idx_usedTb(end+1) = idx_collectedTb;
			else
				continue;
			end
		end

        % if 89GHz on AMSR2 exist, it will be only used if there is no 183 GHz
        % Below algorithm assumes that if a Tb file of AMSR2 is collected, 18.7GHzV-Pol & 89GHzV-PolA-Scan & 89GHzV-PolB-Scan all exist.
        if sum("AMSR2" == sensor_overpass) == 0 % No AMSR2 Tb
            Overpass_write(idx_usedTb,istorm,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,file_overpass,control); 
        else % AMSR2 provides one of overpasses
            idx_order_AMSR2 = find("AMSR2" == sensor_overpass);
            idx_order_other = order_overpass(order_overpass ~= idx_order_AMSR2);
        
            num_183GHz = 0;    
            ChName_other = ChName_all{idx_usedTb(idx_order_other)}    
            for it =1:length(ChName_other)
                if contains(ChName_other{it},'183') | contains(ChName_other{it},'190')  
                    num_183GHz = num_183GHz + 1;
                end
            end
        
            if num_183GHz == 0 % AMSR2 exists and only low frequency of other files are used: use all of frequencies of AMSR2. (? even low frequency?)
                control.comnine_AMSR89GHz = true;
                control.NOTuse_AMSR89GHz = false; 
                Overpass_write(idx_usedTb,istorm,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,file_overpass,control);
            else % AMSR2 exist and other files with ~ 183 GHz exist. Only low frequency of AMSR2 is used.
                control.comnine_AMSR89GHz = false;
                control.NOTuse_AMSR89GHz = true;
                Overpass_write(idx_usedTb,istorm,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,file_overpass,control);
            end
        end

    end % end loop for io = 1:length(overpass)


end % end loop for istorm = 1:length(control.storm_phase)



	

% ---------- Background Knowledge of L1C file----------------
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











































































% This script reads satvvellite MW Tbs observations and if those Tbs cover the storm at DA time then process them into records that are readable by Penn State EnKF system


% -------------- Control variables ----------------------
% By default, control variables prefer cell-array type.
control = struct;
% ----Path
control.obs_dir = '../../Obs/Microwave/';
control.obs_used_dir = '../../Obs/MW_used/'; 
control.bestrack_dir = '../../Obs/Bestrack/';
control.output_dir = '../../toEnKFobs/';
% ---Storm information
control.storm_phase = {'Irma2ndRI',};
%control.storm_phase = ["Irma2ndRI",'JoseRI','MariaRI'};
control.period = {{'201709030600','201709050600'},};
%control.period = {{'201709030600','201709050600'},{'201709050600','201709070600'},{'201709160000','201709180000'}}; %YYYYMMDDHHmm
% ---Satellite informaiton
control.sensor = {'AMSR2','ATMS','GMI','MHS','SAPHIR','SSMIS'};
control.platform = {{'GCOMW1'}, {'NPP'}, {'GPM'}, {'METOPA','METOPB','NOAA18','NOAA19'}, {'MT1'}, {'F16','F17','F18'}};

control.favCh = {'18.7GHzV-Pol','183.31+/-7GHzV-Pol','183.31+-7GHzH-Pol'}; % favorite frequencies
control.favCh_sup = {'19.35GHzV-Pol','89GHzV-PolA-Scan','89GHzV-PolB-Scan','183.31+/-6.6GHzH-Pol','183.31+/-6.8GHz','190.31GHzV-Pol'};
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


% ---------- Loop through each storm_phase -------------------
for istorm = 1:length(control.storm_phase)

    % --- Find time and location of the storm in Best track file
	bestrack_str = Bestrack_read(istorm, control);

    % --- Gather useful files of all sensors into a directory
	[Tbfile_names,Swath_used,ChIdx_ps,ChName_ps,if_swath_good,DAtime_ps,loc_DAtime_ps,overpass,singlepass] = Gather_MW_useful(istorm, bestrack_str, control); % ps: per swath

    % --- Loop through each useful Tb file via a symbolic link
	Tb_dir = [control.obs_used_dir,control.storm_phase{istorm},'/*'];
	%Tb_files = regexprep(ls(Tb_dir),'\n$', ' ');
	Tb_files = regexprep(ls(Tb_dir),'\s', '\n');
    Tb_files = regexp(Tb_files,'\n','split');
	% --- Output file under two situations: overpass or single-pass
    % The order of gathering useful Tb files is different from the order of listing gathered Tb files
    % idx_gatheredTb records the order of Tbs gathered from different sensors & platforms
    % iTb indicates the order of gathered Tb files with ls command
    % - Output single-pass
    for is = 1:length(singlepass)
        for iTb = 1:length(Tb_files)
            Tb_file = Tb_files{iTb}
            [filepath,filename,filext] = fileparts(Tb_file);            
            if contains(filename,singlepass(is))
                idx_gatheredTb = find([filename,filext] == Tbfile_names)
                Singlepass_write(idx_gatheredTb,istorm,Swath_used,ChIdx_ps,ChName_ps,if_swath_good,DAtime_ps,loc_DAtime_ps,Tb_file,control);
            else
                continue;
            end
        end

    end
	% - Output overpass
	Tb_overpass = []; % (strings)
	idx_usedTb = []; % (integer)
	for io = 1:length(overpass_mark) 
		for iTb = 1:length(Tb_files)
			Tb_file = Tb_files{iTb};
            [filepath,filename,filext] = fileparts(Tb_file);
			% gather names of overpass Tb files
			if contains(filename,overpass(io))
                idx_gatheredTb = find([filename,filext] == Tbfile_names);
				Tb_overpass = [Tb_overpass,string(Tb_file)];
				idx_usedTb(end+1) = idx_gatheredTb;
			else
				continue;
			end
		end
		
		Overpass_write(Tb_overpass,idx_usedTb,Swath_used,ChIdx_ps,ChName_ps,if_swath_good,DAtime_ps,loc_DAtime_ps,Tb_file,control);			

	end	
end




	

% ---------- Background Knowledge ----------------
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











































































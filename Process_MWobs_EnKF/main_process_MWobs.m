% This script reads satellite MW Tbs observations and if those Tbs cover the storm at DA time then gets them geographyically plotted.

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
% String: Double quoted. Eg. str = ["123","456"] --> str = 1×2 string array "123"    "3456"

% -------------- Control variables ----------------------
% By default, control variables prefer cell-array type.
control = struct;
% ----Path
control.obs_dir = '../../Obs/Microwave/';
control.obs_used_dir = '../../Obs/MW_used/';
control.bestrack_dir = '../../Obs/Bestrack/';
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
control.filter_ratio = [6,6];
control.roi_oh = {[200,0], [60,60]}; % roi [other variables, hydrometeors]


% ---------- Loop through each storm_phase -------------------
for istorm = 1:length(control.storm_phase)

    % --- Find time and location of the storm in Best track file
	bestrack_str = Bestrack_read(istorm, control);

    % --- Gather useful files of all sensors into a directory
	[Tbfile_name,Swath_used,ChIdx_ps,ChName_ps,if_swath_good,DAtime_ps,loc_DAtime_ps] = Gather_MW_useful(istorm, bestrack_str, control); % ps: per swath

    % --- loop through each useful Tb file via a symbolic link
	Tb_dir = [control.obs_used_dir,'/*'];
	Tb_files = regexprep(ls(Tb_dir),'\n$', '');
	Tb_files = regexp(Tb_files,'\n','split');
	for i = 1:length(Tb_files)
		Tb_file = Tb_files{i};
		if ~contains(Tb_file, Tbfile_name{i})
			disp('Error matching attributes with the Tb file!');
		end
		disp(['Processing level 1C file: ',Tbfile_name{i},'...............']);	
		disp(['DA Channels: ',ChName_ps{i}(if_swath_good{i})]);
		disp(['DA time: ', DAtime_ps{i}]);
		% --- Get area that we will be getting observations for 
		nx = control.nx*control.domain_buffer;
		ny = control.ny*control.domain_buffer;
		min_XLONG = loc_DAtime_ps{i}(1) - (nx/2*dx)/111;
		max_XLONG = loc_DAtime_ps{i}(1) + (nx/2*dx)/111;
		min_XLAT = loc_DAtime_ps{i}(2) - (ny/2*dx)/(cos(loc_DAtime_ps{i}(2)*(pi/180))*111);
		max_XLAT = loc_DAtime_ps{i}(2) + (ny/2*dx)/(cos(loc_DAtime_ps{i}(2)*(pi/180))*111);
		disp(['min of xlong: ',min_XLONG, ', max of xlong: ',max_XLONG]);
		disp(['min of xlat: ',min_XLAT, ', max of xlat: ',max_XLAT]);
		latitudes  = linspace(min_XLAT,max_XLAT,ny);
		longitudes = linspace(min_XLONG,max_XLONG,nx);
		[XLAT, XLONG] = meshgrid(latitudes,longitudes);	
		% --- 
	





































	end















end

























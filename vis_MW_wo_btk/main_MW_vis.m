% This script reads satellite MW Tbs observations and gets it geographyically plotted.

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

% -------------- Control variables ----------------------
control = struct;

control.obs_dir = '../../Obs/Microwave/';
control.plot_dir = '../../Visual/Obs/Microwave/';
%control.storm_phase = {'Irma1stRI'};
control.storm_phase = {'Irma1stRI','Irma2ndRI','JoseRI','MariaRI'};
control.sensor = {'AMSR2','ATMS','GMI','MHS','SAPHIR','SSMIS'};
control.platform = {{'GCOMW1'}, {'NPP'}, {'GPM'}, {'METOPA','METOPB','NOAA18','NOAA19'}, {'MT1'}, {'F16','F17','F18'}};

control.favCh = {'18.7GHzV-Pol','183.31+/-7GHz','183.31+-7GHz'}; % favorite frequencies
control.favCh_sup = {'19.35GHzV-Pol','89GHzV-PolB-Scan','183.31+/-6.6GHz','183.31+/-6.8GHz','190.31GHz'};


% ---------- Loop through each storm_phase -------------------
for istorm = 1:length(control.storm_phase)
    % Set paramters to control the area of interest
    switch(control.storm_phase{istorm})
        case 'Irma1stRI'
            min_lat = 5;
            max_lat = 20;
            min_lon = -35;
            max_lon = -20;
        case 'Irma2ndRI'
            min_lat = 10;
            max_lat = 30;
            min_lon = -60;
            max_lon = -40;
        case 'JoseRI'
            min_lat = 5;
            max_lat = 20;
            min_lon = -50;
            max_lon = -35;
        case 'MariaRI'
            min_lat = 5;
            max_lat = 20;
            min_lon = -60;
            max_lon = -45;
    end
    % --------- Loop through each sensor ---------------
    for isensor = 1:length(control.sensor)
        plfs_eachsensor = control.platform{isensor};
        % ---- Loop through each platform for each sensor ---
        for isensor_plf = 1:length(plfs_eachsensor)
            Tb_dir = [control.obs_dir, control.storm_phase{istorm}, '/', control.sensor{isensor}, '/', plfs_eachsensor{isensor_plf}, '/*'];
            Tb_files = regexprep(ls(Tb_dir),'\n$', '');
            Tb_files = regexp(Tb_files,'\n','split');
            % ---- Loop through each file for each sensor on a specific platform ---
            for i = 1:length(Tb_files)
                Tb_file = Tb_files{i};
                disp(Tb_file);
                % subroutine to obtain swaths and channel numbers 
                [Swath_used, ChIdx_perSwath, ChName_perSwath] = Swath_Channel(Tb_file, control);
                % subroutine to determine which swaths, if there is any, capture the area of interest
                [if_swath_good] =  if_area_covered(min_lat, max_lat, min_lon, max_lon, Swath_used, Tb_file, control);
                % plot Tbs
                if sum(if_swath_good) == 0
                    disp('No swath captures the area of interest! Skip this file.');
                    continue;
                else
                   MW_vis(control.storm_phase{istorm}, plfs_eachsensor{isensor_plf}, control.sensor{isensor}, min_lat, max_lat, min_lon, max_lon, Swath_used, if_swath_good, ChIdx_perSwath, ChName_perSwath, Tb_file, control);
                end 



            end
        end 
    end
end

























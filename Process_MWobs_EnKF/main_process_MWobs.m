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

% -------------- Control variables ----------------------
control = struct;

control.obs_dir = '../../Obs/Microwave/';
control.obs_used_dir = '../../Obs/MW_used/'
control.bestrack_dir = '../../Obs/Bestrack/';
control.storm_phase = {'Irma2ndRI',}
%control.storm_phase = {'Irma2ndRI','JoseRI','MariaRI'};
control.period = {{'201709030600','201709050600'},}
%control.period = {{'201709030600','201709050600'},{'201709050600','201709070600'},{'201709160000','201709180000'}}; %YYYYMMDDHHmm
control.sensor = {'AMSR2','ATMS','GMI','MHS','SAPHIR','SSMIS'};
control.platform = {{'GCOMW1'}, {'NPP'}, {'GPM'}, {'METOPA','METOPB','NOAA18','NOAA19'}, {'MT1'}, {'F16','F17','F18'}};

control.favCh = {'18.7GHzV-Pol','183.31+/-7GHzV-Pol','183.31+-7GHzH-Pol'}; % favorite frequencies
control.favCh_sup = {'19.35GHzV-Pol','89GHzV-PolA-Scan','89GHzV-PolB-Scan','183.31+/-6.6GHzH-Pol','183.31+/-6.8GHz','190.31GHzV-Pol'};


% ---------- Loop through each storm_phase -------------------
for istorm = 1:length(control.storm_phase)

    % --- Find time and location of the storm in Best track file
    bestrack_dir = [control.bestrack_dir,control.storm_phase{istorm},'/*'];
    bestrack_dir = ls(bestrack_dir);
    bestrack_file = regexprep(bestrack_dir,'\n$','');
    fid = fopen(bestrack_file);
    best_track = textscan(fid,'%s','delimiter','');
    fclose(fid);
    % For every record, separate substrings(values)
    best_track_str_vec = string(best_track{1,1}(:));
    len_record = length(best_track{1,1}(:)); % how many records there are
    % Keep necessary information
    bestrack_str = cell(len_record,3); %time,lat,lon
    for ir=1:len_record 
        record_substr = strsplit(best_track_str_vec(ir),',');
        %Time
        bestrack_str{ir,1} = append(erase(record_substr(1,3),' '),'00'); 
        %Latitude 
        bestrack_str{ir,2} = record_substr(1,7); %Latitude for the DTG: 0 - 900 tenths of degrees
        bestrack_str{ir,2} = str2double(strrep(bestrack_str{ir,2},'N',''))/10; % 0 - 90 degree
        %Longitude
        bestrack_str{ir,3} = record_substr(1,8); %Longitude for the DTG: 0 - 1800 tenths of degrees
        if contains(bestrack_str{ir,3},'W')
            bestrack_str{ir,3} =  0-str2double(strrep(bestrack_str{ir,3},'W',''))/10;
        else
            bestrack_str{ir,3} = str2double(strrep(bestrack_str{ir,3},'W',''))/10;
        end
    end 


    % --- Gather files into a directory
    Gather_MW_useful(istorm,bestrack_str,control);



    % --- loop through each useful Tb file
















end

























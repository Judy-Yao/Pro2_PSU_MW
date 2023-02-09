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

control.bestrack_dir = '../../Obs/Bestrack/';
control.obs_dir = '../../Obs/Microwave/';
control.plot_dir = '../../Visual/Obs/MW_with_btk/';
control.storm_phase = {'Irma2ndRI','JoseRI','MariaRI'};
control.period = {{'201709030600','201709050600'},{'201709050600','201709070600'},{'201709160000','201709180000'}}; %YYYYMMDDHHmm
control.sensor = {'AMSR2','ATMS','GMI','MHS','SAPHIR','SSMIS'};
control.platform = {{'GCOMW1'}, {'NPP'}, {'GPM'}, {'METOPA','METOPB','NOAA18','NOAA19'}, {'MT1'}, {'F16','F17','F18'}};

control.favCh = {'18.7GHzV-Pol','183.31+/-7GHzV-Pol','183.31+-7GHzH-Pol'}; % favorite frequencies
control.favCh_sup = {'19.35GHzV-Pol','89GHzV-PolA-Scan','89GHzV-PolB-Scan','183.31+/-6.6GHzH-Pol','183.31+/-6.8GHz','190.31GHzV-Pol'};


% ---------- Loop through each storm_phase -------------------
for istorm = 1:length(control.storm_phase)
    % Set paramters to control the display area
    switch(control.storm_phase{istorm})
        case 'Irma2ndRI'
            min_lat_dpy = 10;
            max_lat_dpy = 30;
            min_lon_dpy = -60;
            max_lon_dpy = -40;
        case 'JoseRI'
            min_lat_dpy = 5;
            max_lat_dpy = 20;
            min_lon_dpy = -50;
            max_lon_dpy = -35;
        case 'MariaRI'
            min_lat_dpy = 5;
            max_lat_dpy = 20;
            min_lon_dpy = -60;
            max_lon_dpy = -45;
    end

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

    % --------- Loop through each sensor ---------------
    for isensor = 1:length(control.sensor)
        plfs_eachsensor = control.platform{isensor};
        % ---- Loop through each platform for each sensor ---
        for isensor_plf = 1:length(plfs_eachsensor)
            disp(['Processing sensor: ', control.sensor{isensor}, ' on platform ', plfs_eachsensor{isensor_plf}]);
            
            Tb_dir = [control.obs_dir, control.storm_phase{istorm}, '/', control.sensor{isensor}, '/', plfs_eachsensor{isensor_plf}, '/*'];
            Tb_files = regexprep(ls(Tb_dir),'\n$', '');
            Tb_files = regexp(Tb_files,'\n','split');
            % ---- Loop through each file for each sensor on a specific platform ---
            for i = 1:length(Tb_files)
                Tb_file = Tb_files{i};  
                disp(Tb_file);
                % subroutine to judge if this Tb is within the period of interest
                [use_Tb_file] = Filter_file_out_period(istorm, Tb_file, control);        
                if use_Tb_file == 0
                    disp('Microwave observations are not within the period of interest! Skip this file.');
                    continue;  
                else
                    % subroutine to obtain swaths and channel indices under each swath
                    [Swath_used, ChIdx_perSwath, ChName_perSwath] = Swath_Channel(Tb_file, control);  
                    % subroutine to identify the best DA time
                    [if_swath_good, DAtime_perSwath, loc_storm_DAtime] = Find_DAtime_loc(bestrack_str,Swath_used,Tb_file, control); 
                    if sum(if_swath_good) == 0
                        disp('Microwave observations do not exist in the area of interest at DA time! Skip this file.');
                        continue;
                    else
                        % display necessary information
                        for item=1:length(Swath_used)
                            if if_swath_good(item)
                                disp([ChName_perSwath{item}, ' can be DAed at ', DAtime_perSwath{item}]);
                            end
                        end 
                        
                        % For each Tb file, read the time from header
                        fileheader = split(h5readatt(Tb_file,'/','FileHeader'),';');
                        for i = 1:length(fileheader)
                           if contains(fileheader{i},'StartGranuleDateTime')
                              time_start = erase(fileheader{i},'StartGranuleDateTime=');
                           elseif contains(fileheader{i},'StopGranuleDateTime')
                              time_end = erase(fileheader{i},'StopGranuleDateTime=');
                           end
                        end                    
                        time_str1 = erase(extractBefore(strtrim(time_start),17),'-');
                        time_str2 = erase(extractBetween(strtrim(time_end),6,16),'-');
                        time = append(time_str1, '-', time_str2);  

                       % ----- Special treatment to AMSR2 ----
                        if contains(control.sensor{isensor}, 'AMSR2')
                           [if_swath_good, ChName_perSwath] = Combine_AMSR2(control.storm_phase{istorm}, plfs_eachsensor{isensor_plf}, control.sensor{isensor}, min_lat_dpy, max_lat_dpy, min_lon_dpy, max_lon_dpy, time, loc_storm_DAtime, Swath_used, if_swath_good, DAtime_perSwath, ChIdx_perSwath, ChName_perSwath, Tb_file, control);
                        end  

                        %  plot Tbs with best-track location at DA time
                        % -------- Loop through each swath -----
                        for i_sw = 1:length(if_swath_good)                        
                            if if_swath_good(i_sw) == false
                                continue;
                            else
                                % ----- Read Variables --------
                                % -Tbs of a certain channel
                                Tb_char = [Swath_used{i_sw},'/Tc'];
                                Tb_Chs = h5read(Tb_file,Tb_char);
                                idx_Tb_oneCh = ChIdx_perSwath(i_sw);
                                Tb_oneCh = squeeze(Tb_Chs(idx_Tb_oneCh,:,:));

                                % -lat/lon for a certian swath
                                lat_char = [Swath_used{i_sw},'/Latitude'];
                                xlat = h5read(Tb_file,lat_char);
                                lon_char = [Swath_used{i_sw},'/Longitude'];
                                xlon = h5read(Tb_file,lon_char);

                                % Masked Tbs, lat, lon based on the area of interest
                                idx_mask = (xlon < max_lon_dpy) & (xlon > min_lon_dpy) & (xlat < max_lat_dpy) & (xlat > min_lat_dpy);
                                lat_need = xlat(idx_mask);
                                lon_need = xlon(idx_mask);
                                Tb_need = Tb_oneCh(idx_mask);
                                % Plot Tbs
                                Plot(control.storm_phase{istorm}, plfs_eachsensor{isensor_plf}, control.sensor{isensor}, min_lat_dpy, max_lat_dpy, min_lon_dpy, max_lon_dpy, time, loc_storm_DAtime{i_sw}, lat_need, lon_need, Tb_need, DAtime_perSwath{i_sw}, ChName_perSwath{i_sw}, control);
                           end
                       end 
                    end
                end
            end 
        end
    end
end

























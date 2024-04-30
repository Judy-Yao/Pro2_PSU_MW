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

control.bestrack_dir = '../../Preprocess_Obs/raw_Obs/Bestrack/';
control.obs_dir = '../../Preprocess_Obs/raw_Obs/Microwave/';
control.plot_dir = '../../Preprocess_Obs/Visual/raw_Obs/MW_with_btk/';
control.storm_phase = {'IRMA'};
control.period = {{'201709030000','201709050000'}}; %YYYYMMDDHHmm
%control.sensor = {'SSMI',};
control.sensor = {'AMSR2','ATMS','GMI','MHS','SAPHIR','SSMI','SSMIS'}; %sensor name
%control.platform = {{'F15'},};
control.platform = {{'GCOMW1'}, {'NPP'}, {'GPM'}, {'METOPA','METOPB','NOAA18','NOAA19'}, {'MT1'}, {'F15'}, {'F16','F17','F18'}}; % platform name (one sensor corresponds to at least one platform)
%control.favFreq = {{'fcdr_tb19v','fcdr_tb85v'},};
control.favFreq = {{'18.7GHzV-Pol','89GHzV-PolA-Scan','89GHzV-PolB-Scan'},{'183.31+-7GHzQH-Pol'},{'18.7GHzV-Pol','183.31+/-7GHzV-Pol'},{'190.31GHzV-Pol'},{'183.31+/-6.8GHz'},{'fcdr_tb19v','fcdr_tb85v'},{'19.35GHzV-Pol','183.31+/-6.6GHzH-Pol'}};


% ---------- Loop through each storm_phase -------------------
for istorm = 1:length(control.storm_phase)
    % Set paramters to control the display area
    switch(control.storm_phase{istorm})
        case 'IRMA'
            min_lat_dpy = 10;
            max_lat_dpy = 30;
            min_lon_dpy = -60;
            max_lon_dpy = -40;
        %case 'JoseRI'
        %    min_lat_dpy = 5;
        %    max_lat_dpy = 20;
        %    min_lon_dpy = -50;
        %    max_lon_dpy = -35;
        case 'MARIA'
            min_lat_dpy = 5;
            max_lat_dpy = 20;
            min_lon_dpy = -60;
            max_lon_dpy = -40;
    end


    % ============================================================================================================
    % Find times and center locations of the storm in the Best track file within the date range of the case study
    % (These data are only available on 00, 06, 12, 18 UTC)
    % ============================================================================================================

    bestrack_str = Bestrack_read(istorm, control); % (cell: {time, lat, lon})


    % ============================================================================================================
    % Loop1: each sensor; Loop2: each platform; Loop3: each Tb file; Loop4: each channel
    % Plot the geographical distribution of Tb
    % ============================================================================================================

    % --------- Loop through each sensor ---------------
    for isensor = 1:length(control.sensor)
        plfs_eachsensor = control.platform{isensor};
        % ---- Loop through each platform for each sensor ---
        for isensor_plf = 1:length(plfs_eachsensor)
            disp(['Processing sensor: ', control.sensor{isensor}, ' on platform ', plfs_eachsensor{isensor_plf}]);

            if strcmp(control.sensor{isensor}, "SSMI")
                Tb_dir = [control.obs_dir, control.storm_phase{istorm}, '/', control.sensor{isensor}, '/', plfs_eachsensor{isensor_plf}, '/*.nc'];
            else
                Tb_dir = [control.obs_dir, control.storm_phase{istorm}, '/', control.sensor{isensor}, '/', plfs_eachsensor{isensor_plf}, '/*.HDF5'];
            end

            Tb_files = strsplit(ls(Tb_dir));
            Tb_files = Tb_files(cellfun(@isempty, Tb_files) == 0); % Get rid of the annyoing empty cell
            
            % ---- Loop through each file for each sensor on a specific platform ---
            for iTb = 1:length(Tb_files)
                Tb_file = Tb_files{iTb};  
                disp(Tb_file);
                disp(['  Examining Tb file: ',Tb_file]);
                % subroutine to judge if this Tb is within the period of interest
                [use_Tb_file] = Filter_file_out_period(istorm, Tb_file, control);        
                if use_Tb_file == 0
                    disp('Microwave observations are not within the period of interest! Skip this file.');
                    continue;  
                else
                    % obtain swaths, channel/frequency index under each swath, and frequency(ies) name(s) of interest
                    [Swath_used, ChIdx_perCh, ChName_perCh] = Match_Freq(isensor, Tb_file, control); % (strings) (integer) (strings)
                    % subroutine to identify the best DA time for each item
                    [if_swath_good, DAtime_perCh, loc_storm_DAtime] = Find_DAtime_loc(bestrack_str,Swath_used,Tb_file, control); % (logical) (strings) (cell: {double vector})
                    if sum(if_swath_good) == 0
                        disp('    Microwave observations do not exist in the area of interest at DA time! Skip this file.');
                        continue;
                    else
                        DAtime = unique(DAtime_perCh(DAtime_perCh ~= ""));
                        if length(DAtime) > 1
                            error('      DA time is not unique!');
                        end
                        disp(['      (More) microwave observations exist in the area of interest at DA_time ',DAtime{1},'!']);
 
                        % For each Tb file, read the time from header
                        [~,~,filext] = fileparts(Tb_file);
                        if contains(filext,"HDF5")
                            fileheader = split(h5readatt(Tb_file,'/','FileHeader'),';');
                            for i = 1:length(fileheader)
                                if contains(fileheader{i},'StartGranuleDateTime')
                                    time_start = erase(fileheader{i},'StartGranuleDateTime=');
                                elseif contains(fileheader{i},'StopGranuleDateTime')
                                    time_end = erase(fileheader{i},'StopGranuleDateTime=');
                                end
                            end                    
                        elseif contains(filext,"nc")
                            time_start = ncreadatt(Tb_file,'/','time_coverage_start');
                            time_end = ncreadatt(Tb_file,'/','time_coverage_end');
                        end
                        time_str1 = erase(extractBefore(strtrim(time_start),17),'-');
                        time_str2 = erase(extractBetween(strtrim(time_end),6,16),'-');
                        time = append(time_str1, '-', time_str2);  

                       % ----- Special treatment to AMSR2 ----
                        if contains(control.sensor{isensor}, 'AMSR2')
                           [if_swath_good, ChName_perCh] = Combine_AMSR2(control.storm_phase{istorm}, plfs_eachsensor{isensor_plf}, control.sensor{isensor}, min_lat_dpy, max_lat_dpy, min_lon_dpy, max_lon_dpy, time, loc_storm_DAtime, Swath_used, if_swath_good, DAtime_perCh, ChIdx_perCh, ChName_perCh, Tb_file, control);
                        end  

                        %  plot Tbs with best-track location at DA time
                        % -------- Loop through each swath -----
                        for i_sw = 1:length(if_swath_good)                        
                            if if_swath_good(i_sw) == false
                                continue;
                            else
                                % ----- Read Variables --------
                                if contains(filext,"HDF5")
                                    % -Tbs of a certain channel
                                    Tb_char = [Swath_used{i_sw},'/Tc'];
                                    Tb_Chs = h5read(Tb_file,Tb_char);
                                    idx_Tb_oneCh = ChIdx_perCh(i_sw);
                                    Tb_oneCh = squeeze(Tb_Chs(idx_Tb_oneCh,:,:));
                                    % -lat/lon for a certian swath
                                    lat_char = [Swath_used{i_sw},'/Latitude'];
                                    xlat = h5read(Tb_file,lat_char);
                                    lon_char = [Swath_used{i_sw},'/Longitude'];
                                    xlon = h5read(Tb_file,lon_char);

                                elseif contains(filext,"nc")
                                    % -Tbs of a certain channel
                                    Tb_oneCh = ncread(Tb_file, ChName_perCh(i_sw)); % npixel, nscan
                                    % -lat/lon for a certian swath
                                    xlat = ncread(Tb_file, ['lat_' + Swath_used(i_sw)]); % npixel, nscan
                                    xlon = ncread(Tb_file, ['lon_' + Swath_used(i_sw)]); % npixel, nscan
                               end
 
                                % Masked Tbs, lat, lon based on the area of interest
                                idx_mask = (xlon < max_lon_dpy) & (xlon > min_lon_dpy) & (xlat < max_lat_dpy) & (xlat > min_lat_dpy);
                                lat_need = xlat(idx_mask);
                                lon_need = xlon(idx_mask);
                                Tb_need = Tb_oneCh(idx_mask);
                                % Plot Tbs
                                Plot(control.storm_phase{istorm}, plfs_eachsensor{isensor_plf}, control.sensor{isensor}, min_lat_dpy, max_lat_dpy, min_lon_dpy, max_lon_dpy, time, loc_storm_DAtime{i_sw}, lat_need, lon_need, Tb_need, DAtime_perCh{i_sw}, ChName_perCh{i_sw}, control);
                           end
                       end 
                    end
                end
            end 
        end
    end
end

























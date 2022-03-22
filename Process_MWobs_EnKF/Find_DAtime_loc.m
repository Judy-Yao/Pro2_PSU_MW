function [if_swath_good,DAtime_perSwath,loc_storm_DAtime] = Find_DAtime_loc(bestrack_str,Swath_used,Tb_file, control)
% Identify the DA_time (YYYYYMMDDHH00) for a Tb observation file 

% ------------------Algorithm--------------------
% Designed by Zhu (Judy) Yao in February 2022---------

% ------------------- Step 1 ----------------------------
%           start_datetime       end_datetime
% |__________|____!______|__________|_!_________|
% h-1     h-0.5          h        h+0.5        h+1
% For example, for a Tb obvervation file, if StartGranuleDateTime is 2017089170434 and EndGranuleDateTime is
% 2017089170616, the candidates for DA_time are 2017089170500 and 2017089170600.
% Genericlly, for DA_time = 05UTC, effective Tb observations should be limited to those with scan time 
% ranging from 4:30 to 5:30 with one-hour interval. However in this case,for DA_time = 05UTC, effective 
% observations are limited from 4:34 to 5:30; for DA_time = 06UTC, effective observations are limited from 5:30 to 6:16.

    % -read StartGranuleDateTime and StopGranuleDateTime 
    fileheader = split(h5readatt(Tb_file,'/','FileHeader'),';');
    for i = 1:length(fileheader)
        if contains(fileheader{i},'StartGranuleDateTime') 
            time_start = extractBefore(erase(fileheader{i},'StartGranuleDateTime='),18); % only save year-month-date-hour-minutes
        elseif contains(fileheader{i},'StopGranuleDateTime') 
            time_end = extractBefore(erase(fileheader{i},'StopGranuleDateTime='),18); % only save year-month-date-hour-minutes
        end
    end

    % -convert saved time strings to scalar datetime arrays
    start_datetime = datetime(strtrim(time_start),'InputFormat','yyyy-MM-dd''T''HH:mm', 'TimeZone','UTC');
    end_datetime = datetime(strtrim(time_end),'InputFormat','yyyy-MM-dd''T''HH:mm', 'TimeZone','UTC');

    % -Get every possible DA_time (especially hour) and its effective Tb observations
    if day(start_datetime) == day(end_datetime) % at the same day
        DA_year = year(start_datetime);
        DA_month = month(start_datetime);
        DA_day = day(start_datetime);
        % list the candidates for DA_hour which should lie between
        % StartGranuleDateTime and StopGranuleDateTime
        start_DAhh_cand = hour(start_datetime)+1; 
        end_DAhh_cand = hour(end_datetime);
        DAhh_cand = start_DAhh_cand:end_DAhh_cand; 

        DA_per_hhcand = cell(length(DAhh_cand),3); %*****

        for i_cand = 1:length(DAhh_cand)
            DA_per_hhcand{i_cand,1} = datetime(DA_year,DA_month,DA_day,DAhh_cand(i_cand),0,0,'TimeZone','UTC'); % DA_time candidate

            % start_obs_time
            if i_cand == 1 % DA_time candidate nearest to StartGranuleDateTime
                if minute(start_datetime) > 30 
                    DA_per_hhcand{i_cand,2} = start_datetime;
                else
                    DA_per_hhcand{i_cand,2} =DA_per_hhcand{i_cand,1}-minutes(30);
                end
            else
                DA_per_hhcand{i_cand,2} = DA_per_hhcand{i_cand,1}-minutes(30);
            end   

            % end_obs_time
            if i_cand == length(DAhh_cand)
                if minute(end_datetime) < 30
                    DA_per_hhcand{i_cand,3} = end_datetime;
                else
                    DA_per_hhcand{i_cand,3} = DA_per_hhcand{i_cand,1}+minutes(30);
                end        
            else
                DA_per_hhcand{i_cand,3} = DA_per_hhcand{i_cand,1}+minutes(30);
            end
        end
    else %not at the same day 
    % It happens around midnight; It defenitely span two consecutive days; It might span two months
        DA_year = year(start_datetime);
        % list the candidates for DA_hour which should lie between
        % StartGranuleDateTime and StopGranuleDateTime
        start_DAhh_cand = hour(start_datetime)+1; 
        end_DAhh_cand = hour(end_datetime)+24; % 00Z --> 24;01Z --> 25 
        DAhh_cand = start_DAhh_cand:end_DAhh_cand;

        DA_per_hhcand = cell(length(DAhh_cand),3);

        for i_cand = 1:length(DAhh_cand)
            if DAhh_cand >= 24
              DA_month = month(end_datetime); 
              DA_day = day(end_datetime);
              DA_hour = DAhh_cand(i_cand) - 24;
            else
              DA_month = month(start_datetime); 
              DA_day = day(start_datetime);
              DA_hour = DAhh_cand(i_cand);
            end  

            DA_per_hhcand{i_cand,1} = datetime(DA_year,DA_month,DA_day,DA_hour,0,0,'TimeZone','UTC'); % DA_time candidate

            % start_obs_time
            if i_cand == 1 % DA_time candidate nearest to StartGranuleDateTime
                if minute(start_datetime) > 30 
                    DA_per_hhcand{i_cand,2} = start_datetime;
                else
                    DA_per_hhcand{i_cand,2} = DA_per_hhcand{i_cand,1}-minutes(30);
                end
            else
                DA_per_hhcand{i_cand,2} = DA_per_hhcand{i_cand,1}-minutes(30);
            end   

            % end_obs_time
            if i_cand == length(DAhh_cand)
                if minute(end_datetime) < 30
                    DA_per_hhcand{i_cand,3} = end_datetime;
                else
                    DA_per_hhcand{i_cand,3} = DA_per_hhcand{i_cand,1}+minutes(30);
                end        
            else
                DA_per_hhcand{i_cand,3} = DA_per_hhcand{i_cand,1}+minutes(30);
            end
        end
    end

    % -display necessary information
    for i_cand = 1:length(DAhh_cand)
        disp(['DA_time candidate: ', datestr(DA_per_hhcand{i_cand,1},'yyyymmddhhMM')]);
        disp(['Effective observations are from ', datestr(DA_per_hhcand{i_cand,2},'yyyymmddHH:MM'),' to ',datestr(DA_per_hhcand{i_cand,3},'yyyymmddHH:MM')]);
    end

    %% ------------------- Step 2 ----------------------------
    % Best-track locations are only available at 00 UTC, 06 UTC, 12 UTC and 18
    % UTC, therefore we need to linearly interpolate the available locations to the DA_time candidates.
    % For each DA_time candidate, we calculate how many observations lie within the sqaure area that is 
    % centered at the best-track location of a DA_time candidate.
    % Determine the best DA_time as the DA_time candidate with most observations in
    % that square area.

    num_useful_scan = [];
    if_swath_good = []; % (logical)
    DAtime_perSwath = {}; % (cell > strings)
    loc_storm_DAtime = {}; % lat,lon

    for i_sw = 1:length(Swath_used)
        % -For a used Swath, get all scans 
        scan_Year_char = [Swath_used(i_sw) + '/ScanTime/Year'];
        scan_Month_char = [Swath_used(i_sw) + '/ScanTime/Month'];
        scan_Day_char = [Swath_used(i_sw) + '/ScanTime/DayOfMonth']; 
        scan_Hour_char = [Swath_used(i_sw) + '/ScanTime/Hour']; 
        scan_Minute_char = [Swath_used(i_sw) + '/ScanTime/Minute']; 

        scan_Year = h5read(Tb_file,scan_Year_char);
        scan_Month = h5read(Tb_file,scan_Month_char);
        scan_Day = h5read(Tb_file,scan_Day_char);
        scan_Hour = h5read(Tb_file,scan_Hour_char);
        scan_Minute = h5read(Tb_file,scan_Minute_char);  

        all_scan_time = datetime(scan_Year,scan_Month,scan_Day,scan_Hour,scan_Minute,0,'TimeZone','UTC');

        num_useful_scan = zeros(length(DAhh_cand),1);
        loc_storm_DAtime_percand = zeros(length(DAhh_cand),2);  
        % ---- loop through DA_time candidate
        for i_cand = 1:length(DAhh_cand)
            DA_datetime = DA_per_hhcand{i_cand,1};
            % Find index_time between start_obs_time and end_obs_time
            idx_all_scan = isbetween(all_scan_time, DA_per_hhcand{i_cand,2},DA_per_hhcand{i_cand,3});
            % ---- Get useful observations ---------------- 
            %find the start time and end time of the 6h interval where this DA_time is within
            DA_year = year(DA_datetime);
            DA_month = month(DA_datetime);
            DA_day = day(DA_datetime);
            DA_hour = hour(DA_datetime);

            ith_6h = fix(DA_hour/6);

            start_6h_hh = ith_6h*6; 
            start_6h_str = append(datestr(DA_datetime,'yyyymmdd'),num2str(start_6h_hh,'%02.f')); 
            if ith_6h == 3
                end_6h = datetime(DA_year,DA_month,DA_day,start_6h_hh,0,0)+ hours(6); % move the date forward for one day and/or one month
                end_6h_str = datestr(end_6h,'yyyymmddhh');

                end_6h_hh = hour(end_6h)+24; % for the convenience of linear interpolation
            else
                end_6h_hh = (ith_6h+1)*6; 
                end_6h_str = append(datestr(DA_datetime,'yyyymmdd'),num2str(end_6h_hh,'%02.f'));
            end
            % find the corresponding best_track locations at start_6h and end_6h
            best_track_times = [start_6h_hh, end_6h_hh];
            for ir=1:size(bestrack_str,1)
                if contains(bestrack_str{ir,1},start_6h_str)
                    ir_start_6h = ir;
                end

                if contains(bestrack_str{ir,1},end_6h_str)
                    ir_end_6h = ir;
                end
            end
            best_track_locations = [bestrack_str{ir_start_6h,2}, bestrack_str{ir_end_6h,2}; bestrack_str{ir_start_6h,3}, bestrack_str{ir_end_6h,3}];
            % linearly interpolate best_track locations at start_6h and end_6h to DA_time
            best_track_weights = [abs((DA_hour)-best_track_times(2)); abs((DA_hour)-best_track_times(1))] / abs(diff(best_track_times));
            loc_storm_DAtime_percand(i_cand,:) = best_track_locations*best_track_weights;
            % find a square area centered at the best-track location at DA_time
            min_lat = loc_storm_DAtime_percand(i_cand,1)-0.1;
            max_lat = loc_storm_DAtime_percand(i_cand,1)+0.1;
            min_lon = loc_storm_DAtime_percand(i_cand,2)-0.1;
            max_lon = loc_storm_DAtime_percand(i_cand,2)+0.1;
            % mask observation pixels that are between the start_obs_time of and
            % the end_obs_time of a DA_time candidate (idx_all_scan)
            pixel_lat_char = [Swath_used(i_sw) + '/Latitude'];
            pixel_lon_char = [Swath_used(i_sw) + '/Longitude']; 
            pixel_lat = h5read(Tb_file,pixel_lat_char);
            pixel_lat = pixel_lat(:,idx_all_scan);
            pixel_lon = h5read(Tb_file,pixel_lon_char);
            pixel_lon = pixel_lon(:,idx_all_scan);
            % determine if any masked observations are within the square area
            % centered at the best-track location at this DA_time candidate
            num_useful_scan(i_cand) = sum((pixel_lon(:) < max_lon) & (pixel_lon(:) > min_lon) & (pixel_lat(:) < max_lat) & (pixel_lat(:) > min_lat)); 
            %num_useful_scan = [num_useful_scan, useful_scan_i_cand];
        end

        % store values
        % if_swath_good(length(Swath_used))
        % DAtime_perSwath(length(Swath_used))
        if sum(num_useful_scan) == 0
            if_swath_good = [if_swath_good, false];
            DAtime_perSwath{end+1} = "";
            loc_storm_DAtime{end+1} = NaN; 
        else
            if_swath_good = [if_swath_good, true];
            idx_DAtime = find(num_useful_scan == max(num_useful_scan));
            if length(idx_DAtime) > 1
                idx_DAtime = idx_DAtime(1);
            end
            DAtime_perSwath{end+1} = datestr(DA_per_hhcand{idx_DAtime,1},'yyyymmddhhMM');
            loc_storm_DAtime{end+1} = loc_storm_DAtime_percand(idx_DAtime,:);
        end   

    end
    
    DAtime_perSwath = string(DAtime_perSwath);
    % sanity check
    if length(if_swath_good) ~= length(Swath_used) | length(DAtime_perSwath) ~= length(Swath_used) | length(loc_storm_DAtime) ~= length(Swath_used)
        disp('Error identifying DA time for each swath/channel!');
    end

end

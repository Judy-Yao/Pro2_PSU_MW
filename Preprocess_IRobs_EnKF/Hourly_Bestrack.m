% Find times and center locations of the storm in Best track file (only available on 00, 06, 12, 18 UTC) &
% linearly interpolate locations at these times to other o'clock

function [ hr_btk,start_time,end_time ] = Hourly_Bestrack(istorm, control)

    % ----- Read relevant information into memory ------------------
    % scan the best-track file 
    bestrack_dir = [control.bestrack_dir,control.storm_phase{istorm},'/*'];
    bestrack_dir = ls(bestrack_dir);
    bestrack_file = regexprep(bestrack_dir,'\n$','');
    fid = fopen(bestrack_file);
    btk_file = textscan(fid,'%s','delimiter','');
    fclose(fid);
    % For every record, separate substrings(values)
    btk_file_str_vec = string(btk_file{1,1}(:));
    len_record = length(btk_file{1,1}(:)); % how many records there are
    % Read necessary information into an array: time_string, time_datetime, lat, lon
    bestrack = cell(len_record,4); 
    for ir=1:len_record
        record_substr = strsplit(btk_file_str_vec(ir),',');
        %Time_string
        bestrack{ir,1} = append(erase(record_substr(1,3),' '),'00'); % %YYYYMMDDHHmm
        %Time_datetime
        bestrack{ir,2} = datetime(bestrack{ir,1},'InputFormat','yyyyMMddHHmm','TimeZone','UTC');
        %Latitude 
        bestrack_lat = record_substr(1,7); %Latitude for the DTG: 0 - 900 tenths of degrees (e.g.,161N)
        if contains(bestrack_str{ir,3},'N')
            bestrack_str{ir,3} = str2double(strrep(bestrack_lat,'N',''))/10; % 0 - 90 degree
        else
            bestrack_str{ir,3} = 0-str2double(strrep(bestrack_lat,'S',''))/10;
        end
        %Longitude
        bestrack_lon = record_substr(1,8); %Longitude for the DTG: 0 - 1800 tenths of degrees (e.g.,269W)
        if contains(bestrack_lon,'W')
            bestrack{ir,4} =  0-str2double(strrep(bestrack_lon,'W',''))/10; 
        else
            bestrack{ir,4} = str2double(strrep(bestrack_lon,'E',''))/10; 
        end
    end

    % ----- Linearly interpolate storm locations to every hour within the period of interest  ---------
    % generate the sequence of time of interest
    time_st_oi = datetime(control.period{istorm}{1},'InputFormat','yyyyMMddHHmm','TimeZone','UTC');
    time_end_oi = datetime(control.period{istorm}{2},'InputFormat','yyyyMMddHHmm','TimeZone','UTC');
    duration_hh =  hours(time_end_oi - time_st_oi);
    DA_datetime_ss = time_st_oi + hours(0:duration_hh); % time with every hour within the period of interest (dateime) 
   
    hr_btk = cell(length(DA_datetime_ss),3);
    for ih = 1:length(DA_datetime_ss)
        DA_datetime = DA_datetime_ss(ih); 
        %find the start time and end time of the 6h interval where this DA_time is within
        DA_year = year(DA_datetime);
        DA_month = month(DA_datetime);
        DA_day = day(DA_datetime);
        DA_hour = hour(DA_datetime);

        ith_6h = fix(DA_hour/6); % 0 <= ith_6h <= 3 

        start_6h_hh = ith_6h*6; 
        start_6h_str = append(datestr(DA_datetime,'yyyymmdd'),num2str(start_6h_hh,'%02.f')); 
        if ith_6h == 3
            end_6h = datetime(DA_year,DA_month,DA_day,start_6h_hh,0,0)+ hours(6); % move the date forward for one day and/or one month
            end_6h_str = datestr(end_6h,'yyyymmddhh');

            end_6h_hh = hour(end_6h)+24; % for the convenience of linear interpolation, we change uadrovigesimal to decimal measuring hours
        else
            end_6h_hh = (ith_6h+1)*6; 
            end_6h_str = append(datestr(DA_datetime,'yyyymmdd'),num2str(end_6h_hh,'%02.f'));
        end
        % match start_6h and end_6h with records in the best-track file
        best_track_times = [start_6h_hh, end_6h_hh];
        for ir=1:size(bestrack,1)
            if bestrack{ir,2} == DA_datetime
               ir_DA_time = ir; 
            end

            if contains(bestrack{ir,1},start_6h_str)
                ir_start_6h = ir;
            end

            if contains(bestrack{ir,1},end_6h_str)
                ir_end_6h = ir;
            end
        end
        % find the center location at start_6h and end_6h 
        best_track_locations = [bestrack{ir_start_6h,3}, bestrack{ir_end_6h,3}; bestrack{ir_start_6h,4}, bestrack{ir_end_6h,4}];
        % linearly interpolate best_track locations at start_6h and end_6h to DA_time
        best_track_weights = [abs((DA_hour)-best_track_times(2)); abs((DA_hour)-best_track_times(1))] / abs(diff(best_track_times));        
        hr_btk{ih,1} = datestr(DA_datetime,'yyyymmddHHMM');% time_string
        hr_btk{ih,2} = DA_datetime; % datetime
        hr_btk{ih,3} = best_track_locations*best_track_weights; 
    end

    start_time = hr_btk{1,1};
    end_time = hr_btk{length(DA_datetime_ss),1};


end

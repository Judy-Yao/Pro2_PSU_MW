% This function simply filters out MW L1C files that don't lie within the period of interest 

function use_Tb_file = Filter_file_out_period(idx_storm,Tb_file,control)

    % -read StartGranuleDateTime and StopGranuleDateTime for a Tb file
    fileheader = split(h5readatt(Tb_file,'/','FileHeader'),';');
    for i = 1:length(fileheader)
        if contains(fileheader{i},'StartGranuleDateTime') 
            time_start = extractBefore(erase(fileheader{i},'StartGranuleDateTime='),18); % only save year-month-date-hour-minutes
			disp(['    Start Granule Time: ',strtrim(time_start)]);
		elseif contains(fileheader{i},'StopGranuleDateTime') 
            time_end = extractBefore(erase(fileheader{i},'StopGranuleDateTime='),18); % only save year-month-date-hour-minutes
			disp(['    End Granule Time: ',strtrim(time_end)]);
		end
    end
    % -convert saved time strings to scalar datetime arrays
    start_datetime = datetime(strtrim(time_start),'InputFormat','yyyy-MM-dd''T''HH:mm', 'TimeZone','UTC');
    end_datetime = datetime(strtrim(time_end),'InputFormat','yyyy-MM-dd''T''HH:mm', 'TimeZone','UTC');
    
    % -convert the start and the end of the period for the storm to scalar datetime arrays
    pd_start_dt = datetime(control.period{idx_storm}{1},'InputFormat','yyyyMMddHHmm','TimeZone','UTC');
    pd_end_dt = datetime(control.period{idx_storm}{2},'InputFormat','yyyyMMddHHmm','TimeZone','UTC');

    % Skip the Tb file that is not within the period of interest
    if (end_datetime < pd_start_dt) || (start_datetime > pd_end_dt) % || expr2 is not evaluated if expr1 is logical 1 (true).
        use_Tb_file = false;
    else
        use_Tb_file = true;
    end

end

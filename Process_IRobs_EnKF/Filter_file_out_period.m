% This function simply filters out MW L1C files that don't lie within the period of interest 

function use_Tb_file = Filter_file_out_period(idx_storm,Tb_file,control)

    % -read time_coverage_start and time_coverage_end from global attributes from a Tb file
    time_start = ncreadatt(Tb_file,'/','time_coverage_start');
    time_end =  ncreadatt(Tb_file,'/','time_coverage_end'); 
   
     % -convert saved time strings to scalar datetime arrays
    start_datetime = replace(extractBefore(time_start,17),'T',' ');
    end_datetime =  replace(extractBefore(time_end,17),'T',' ');
    
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

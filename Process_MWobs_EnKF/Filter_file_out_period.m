% ----------------------------------------------------------------------------------------- %
% This function simply filters out MW L1C files that are not within the period of interest 
% ----------------------------------------------------------------------------------------- %

function use_Tb_file = Filter_file_out_period(idx_storm,Tb_file,control)

    % Convert the start and the end of the period of interest for the storm to scalar datetime arrays
    pd_start_dt = datetime(control.period{idx_storm}{1},'InputFormat','yyyyMMddHHmm','TimeZone','UTC');
    pd_end_dt = datetime(control.period{idx_storm}{2},'InputFormat','yyyyMMddHHmm','TimeZone','UTC');

	% Determine if the file should be read by ncread or h5read
	[~,~,filext] = fileparts(Tb_file);
	
	% Read the beginning and ending of the time coverage for this file
	% HDF5 file
	if contains(filext,"HDF5")
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
		% convert saved time strings to scalar datetime arrays
		start_time = datetime(strtrim(time_start),'InputFormat','yyyy-MM-dd''T''HH:mm', 'TimeZone','UTC');
		end_time = datetime(strtrim(time_end),'InputFormat','yyyy-MM-dd''T''HH:mm', 'TimeZone','UTC'); 
	% NC file
	elseif contains(filext,"nc")
        start_time = datetime(extractBefore(ncreadatt(Tb_file,'/','time_coverage_start'),17),'InputFormat','yyyy-MM-dd''T''HH:mm', 'TimeZone','UTC');
		end_time = datetime(extractBefore(ncreadatt(Tb_file,'/','time_coverage_end'),17),'InputFormat','yyyy-MM-dd''T''HH:mm', 'TimeZone','UTC');
	end

    % Skip the Tb file that is not within the period of interest
    if (end_time < pd_start_dt) || (start_time > pd_end_dt) % || expr2 is not evaluated if expr1 is logical 1 (true).
        use_Tb_file = false;
    else
        use_Tb_file = true;
    end

end

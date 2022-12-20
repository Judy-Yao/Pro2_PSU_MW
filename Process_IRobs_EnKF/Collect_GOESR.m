function [] = Collect_GOESR(istorm, bestrack_str, control)

    % Determine if a subdirectory for a storm object exists
    if ~exist([control.obs_collect_dir,control.storm_phase{istorm}],'dir')
        [~, msg, ~] = mkdir(control.obs_collect_dir,control.storm_phase{istorm});% create a subdirectory for a storm object under Collected_MW/
        if isempty(msg)
            disp(['  Successfully created a subdirectory in ',control.obs_collect_dir,' for ',control.storm_phase{istorm}]); 
        else
            error('  Error: ',msg);
        end
        destination = [control.obs_collect_dir,control.storm_phase{istorm},'/'];
    else
        % Clean existed symbolic files
        destination = [control.obs_collect_dir,control.storm_phase{istorm},'/'];
        [status,~] = system(['rm ',destination,'*']); % Tricky part: it seems that it can't delete unsuccessful links!!
        if status ~= 0
            disp('  Warning: cleaning existed symbolic files failed!');
        end
    end 
    
    all_Tbfile_name = []; % (strings)
    time_btk_str = strings(size(bestrack_str,1),1);
    for ibr = 1:size(bestrack_str,1)
        time_btk_str(ibr) = convertCharsToStrings(bestrack_str{ibr,1});
    end
    % Loop through files for each channel
    for ich = 1:length(control.favCH) 
        Tb_dir = [control.obs_dir, control.storm_phase{istorm},'/*.nc'];
        %Tb_dir = [control.obs_dir, control.storm_phase{istorm},'/CH',num2str(control.favCH(ich)),'/*.nc'];
        Tb_files = strsplit(ls(Tb_dir));
        Tb_files = Tb_files(cellfun(@isempty, Tb_files) == 0); % Get rid of the annyoing empty cell

        for i = 1:length(Tb_files)
            Tb_file = Tb_files{i};
            disp(['  Examining Tb file: ',Tb_file]);
            % determine if this Tb is within the period of interest
            [use_Tb_file] = Filter_file_out_period(istorm, Tb_file, control); %(logical)
            if use_Tb_file == 0
                continue;
            else
                % read time_coverage_start from global attributes from the Tb file
                time_start = ncreadatt(Tb_file,'/','time_coverage_start');
                time_start = convertCharsToStrings(erase(extractBefore(time_start,17),["T","-",":"]));
                % find the location of this date/hour in bestrack_str
                ir = find(time_start == time_btk_str);
                if isempty(ir) 
                    continue;
                else
                    [filepath,filename,filext] = fileparts(Tb_file);
                    source_file = erase(Tb_file,'raw_Obs/');
                    newfile_name = ['DAt' + time_btk_str(ir) + '_GOESR_ABI_L2_CMIPF_CH' + num2str(control.favCH(ich)) + '.nc'];
                    command = ["ln -s " + source_file + " " + destination + newfile_name];% Symbolic link's path -> source file. The relatively path of source file is relative to symbolic link's path.
                    [status,~] = system(command);
                    if status == 0
                        disp(['  ',filename,filext, ' is collected.']);
                    else
                        disp('  Warning: collecting the Tb file might have failed!');
                    end
                    
                end 
            end
        end 
    end






end

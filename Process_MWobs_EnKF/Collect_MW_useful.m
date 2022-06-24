% This function loops through all level 1c files directly downloaded from NASA GES DISC website for a storm object 
% It collects files that meet the requiremens and symbolically links these files to a directory (/Obs/Collected_MW/)
% It also returns useful attributes for each useful Tb files

function [all_Tbfile_name,all_Swath_used,all_ChIdx_perCh,all_ChName_perCh,all_DAtime_perCh,all_loc_storm_DAtime,overpass,singlepass] = Collect_MW_useful(istorm, bestrack_str, control)

    % Determine if a subdirectory for a storm object exists
	if ~exist([control.obs_collect_dir,control.storm_phase{istorm}],'dir')
	    [~, msg, ~] = mkdir(control.obs_collect_dir,control.storm_phase{istorm});% create a subdirectory for a storm object under Collected_MW/
		if isempty(msg)
			disp(['Successfully created a subdirectory in ',control.obs_collect_dir,' for ',control.storm_phase{istorm}]); 
		else
			error('Error: ',msg);
		end
		destination = [control.obs_collect_dir,control.storm_phase{istorm},'/'];
	else
		% Clean existed symbolic files
		destination = [control.obs_collect_dir,control.storm_phase{istorm},'/'];
		[status,~] = system(['rm ',destination,'*']); % Tricky part: it seems that it can't delete unsuccessful links!!
		if status ~= 0
			disp('Warning: cleaning existed symbolic files failed!');
		end
	end
 
	% Define empty cells (final length of each cell is the number of useful Tbs)
	all_Tbfile_name = []; % (strings)
	all_Swath_used = {};
	all_ChIdx_perCh = {};
	all_ChName_perCh = {};
	all_if_swath_good = {};
	all_DAtime_perCh = {};
	all_loc_storm_DAtime = {};
    
    % --------- Loop through each sensor ---------------
    for isensor = 1:length(control.sensor)
        plfs_eachsensor = control.platform{isensor};
        % ---- loop through each platform for a sensor ---
        for isensor_plf = 1:length(plfs_eachsensor)
            disp(['Processing sensor: ', control.sensor{isensor}, ' on platform ', plfs_eachsensor{isensor_plf}]);

            Tb_dir = [control.obs_dir, control.storm_phase{istorm}, '/', control.sensor{isensor}, '/', plfs_eachsensor{isensor_plf}, '/*'];
            Tb_files = strsplit(ls(Tb_dir));
            Tb_files = Tb_files(cellfun(@isempty, Tb_files) == 0); % Get rid of the annyoing empty cell
            %Tb_files = regexprep(ls(Tb_dir),'\n$', '');
            %Tb_files = regexp(Tb_files,'\n','split');
            % ---- Loop through each file for a sensor on a specific platform ---
            for i = 1:length(Tb_files)
                Tb_file = Tb_files{i};
                disp(['  Examining Tb file: ',Tb_file]);
                % determine if this Tb is within the period of interest
                [use_Tb_file] = Filter_file_out_period(istorm, Tb_file, control); %(logical)
                if use_Tb_file == 0
                    disp('  Microwave observations are not within the period of interest! Skip this file.');
                    continue;
                else
                    % obtain swaths, channel/frequency index under each swath, and frequency(ies) name(s) of interest
                    [Swath_used, ChIdx_perCh, ChName_perCh] = Match_Freq(isensor, Tb_file, control); % (strings) (double) (strings)
					% subroutine to identify the best DA time for each item
                    [if_swath_good, DAtime_perCh, loc_storm_DAtime] = Find_DAtime_loc(bestrack_str,Swath_used,Tb_file, control); % (logical) (strings) (cell{double})
                    if sum(if_swath_good) == 0
                        disp('	Microwave observations do not exist in the area of interest at DA time! Skip this file.');
                        continue;
                    else
						% ----- Gather the useful Tb file of all sensors !!
                        source_file = erase(Tb_file,'raw_Obs/');
                        % A potential BUG exists: the current algorithm assumes that for all channels of a L1C MW file the best-track locations and DA times are the same
						if (length(DAtime_perCh) > 1) & (DAtime_perCh(1) ~= DAtime_perCh(2))
							disp('	Error collecting the Tb file! Potential risk exists!');
						end
						[filepath,filename,filext] = fileparts(Tb_file);
						if contains(filext,"HDF5")
							newfile_name = ['DAt' + DAtime_perCh(1) + '_1C.' + control.platform{isensor}{isensor_plf} + '.' + control.sensor{isensor} + '.HDF5'];
						elseif contains(filext,"nc")
							newfile_name = ['DAt' + DAtime_perCh(1) + '_1C.' + control.platform{isensor}{isensor_plf} + '.' + control.sensor{isensor} + '.nc'];	
						end                    
						command = ["ln -s " + source_file + " " + destination + newfile_name];% Symbolic link's path -> source file. The relatively path of source file is relative to symbolic link's path.
                        [status,~] = system(command);
                        if status == 0
                            disp(['  ',filename,filext, ' is collected.']);
                        else
                            disp('  Error collecting the Tb file!');
                        end
						% store useful information with correction from if_swath_good
                        all_Tbfile_name = [all_Tbfile_name,newfile_name]; % (strings)
						if_swath_good = logical(if_swath_good);
						all_Swath_used{end+1} = Swath_used(if_swath_good); % (cell{strings})
						all_ChIdx_perCh{end+1} = ChIdx_perCh(if_swath_good); % (cell{single})
						all_ChName_perCh{end+1} =  ChName_perCh(if_swath_good); % (cell{strings})
						all_DAtime_perCh{end+1} = DAtime_perCh(if_swath_good); % (cell{strings})
						all_loc_storm_DAtime{end+1} = loc_storm_DAtime{if_swath_good}; % (cell{cell{double}})
                    end
                end
            end
        end
    end
		

	% Sanity check
	if (length(all_Tbfile_name) ~= length(all_Swath_used)) | (length(all_Tbfile_name) ~= length(all_ChName_perCh))
		disp('  Error collecting useful attributes for Tb files!');
	end

	% --- Mark satellite overpass and single-pass
	% Only works only if Each Tb file has only one DAtime 
	overpass = [];
	singlepass = [];
	strs_date = [];
	for f = 1:length(all_Tbfile_name)
		strs_date =  [strs_date,all_DAtime_perCh{f}(1)];
	end
	unique_date = unique(strs_date);
	for id = 1:length(unique_date)
		repeate_times = sum(strs_date == unique_date(id));
		if repeate_times == 1
			singlepass = [singlepass,unique_date(id)];
		else
			overpass = [overpass,unique_date(id)];
		end
	end


end

% Warning! 
% If a string is only consisted of number 
% {test = [];test(end+1) = string} will automatically convert string to the number
% In this case, it is better to use {test = [];test = [test,string]}













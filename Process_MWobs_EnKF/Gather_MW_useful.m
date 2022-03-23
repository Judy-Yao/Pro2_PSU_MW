function [all_Tbfile_name,all_Swath_used,all_ChIdx_perSwath,all_ChName_perSwath,all_if_swath_good,all_DAtime_perSwath,all_loc_storm_DAtime,overpass,singlepass] = Gather_MW_useful(istorm, bestrack_str, control)

% This function loops through all level 1c files directly downloaded from NASA GES DISC website 
% It picks files that meet the requirement and symbolically links these files to a directory
% It also returns 7 useful attributes for each useful Tb files
 
    % Clean existing symbolic files 
    destination = [control.obs_used_dir,control.storm_phase{istorm},'/']; 
    [status,~] = system(['rm ',destination,'*']); % Tricky part: it seems that it can't delete unsuccessful links
	% Define empty cells (final length of each cell is the number of useful Tbs)
	all_Tbfile_name = []; % (strings)
	all_Swath_used = {};
	all_ChIdx_perSwath = {};
	all_ChName_perSwath = {};
	all_if_swath_good = {};
	all_DAtime_perSwath = {};
	all_loc_storm_DAtime = {};
    
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
                % subroutine to determine if this Tb is within the period of interest
                [use_Tb_file] = Filter_file_out_period(istorm, Tb_file, control); %(logical)
                if use_Tb_file == 0
                    disp('Microwave observations are not within the period of interest! Skip this file.');
                    continue;
                else
                    % subroutine to obtain swaths and channel indices under each swath
                    [Swath_used, ChIdx_perSwath, ChName_perSwath] = Swath_Channel(Tb_file, control); % (strings) (double) (strings)
                    % subroutine to identify the best DA time for each item
                    [if_swath_good, DAtime_perSwath, loc_storm_DAtime] = Find_DAtime_loc(bestrack_str,Swath_used,Tb_file, control); % (logical) (strings) (cell{double})
                    if sum(if_swath_good) == 0
                        disp('Microwave observations do not exist in the area of interest at DA time! Skip this file.');
                        continue;
                    else
						% ----- Gather the useful Tb file of all sensors !!
                        [filepath,filename,filext] = fileparts(Tb_file);
                        source_file = erase(Tb_file,'Obs/');
                        % A potential BUG exists which doesn't need to be taken care of at this point
						% The current code assumes that for all channels the best-track location and DA time are the same
						if (length(DAtime_perSwath) > 1) & (DAtime_perSwath(1) ~= DAtime_perSwath(2))
							disp('Error renaming the Tb file! Potential risk exists!');
						end
						newfile_name = ['DAt' + DAtime_perSwath(1) + '_1C.' + control.platform{isensor}{isensor_plf} + '.' + control.sensor{isensor} + '.HDF5'];
                        command = ["ln -s " + source_file + " " + destination + newfile_name];% Symbolic link's path -> source file. The relatively path of source file is relative to symbolic link's path.
                        [status,~] = system(command);
                        if status == 0
                            disp([filename,filext, ' is gathered.']);
                        else
                            disp('Error gathering Tb file!');
                        end
						% store useful information
                        all_Tbfile_name = [all_Tbfile_name,newfile_name]; % (strings)
						all_Swath_used{end+1} = Swath_used; % (cell{strings})
						all_ChIdx_perSwath{end+1} = ChIdx_perSwath; % (cell{single})
						all_ChName_perSwath{end+1} = ChName_perSwath; % (cell{strings})
						all_if_swath_good{end+1} = if_swath_good; % (cell{logical})
						all_DAtime_perSwath{end+1} = DAtime_perSwath; % (cell{strings})
						all_loc_storm_DAtime{end+1} = loc_storm_DAtime; % (cell{cell{double}})
                    end
                end
            end
        end
    end

	% Sanity check
	if (length(all_Tbfile_name) ~= length(all_Swath_used)) | (length(all_Tbfile_name) ~= length(all_if_swath_good))		
		disp('Error gathering useful attributes for Tb files!');
	end

	% --- Mark satellite overpass and single-pass
	% Only works only if Each Tb file has only one DAtime 
	overpass = [];
	singlepass = [];
	strs_date = [];
	for f = 1:length(all_Tbfile_name)
		strs_date =  [strs_date,all_DAtime_perSwath{f}(1)];
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













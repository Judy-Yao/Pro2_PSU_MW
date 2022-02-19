function [] = Gather_MW_useful(istorm, bestrack_str, control)
    
    % -------- Clean existing symbolic files -----------
    destination = [control.obs_used_dir,control.storm_phase{istorm},'/']; 
    [status,~] = system(['rm ',destination,'*']); % Tricky part: it seems that it can't delete unsuccessful links
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
                        [filepath,filename,filext] = fileparts(Tb_file);
                        source_file = erase(Tb_file,'Obs/');
                        command = ['ln -s ',source_file,' ',destination,filename,filext]; % Symbolic link's path -> source file. The relatively path of source file is relative to symbolic link's path.
                        [status,~] = system(command);
                        if status == 0
                            disp([filename,filext, ' is gathered.']);
                        else
                            disp('Error gathering Tb file!');
                        
                        end

                        % display necessary information
                       % for item=1:length(Swath_used)
                       %     if if_swath_good(item)
                       %         disp([ChName_perSwath{item}, ' can be DAed at ', DAtime_perSwath{item}]);
                       %     end
                       % end

                       % ----- Special treatment to AMSR2 ----
                        %if contains(control.sensor{isensor}, 'AMSR2')
                        %[if_swath_good, ChName_perSwath] = Combine_AMSR2(control.storm_phase{istorm}, plfs_eachsensor{isensor_plf}, control.sensor{isensor}, min_lat_dpy, max_lat_dpy, min_lon_dpy, max_lon_dpy, time, loc_storm_DAtime, Swath_used, if_swath_good, DAtime_perSwath, ChIdx_perSwath, ChName_perSwath, Tb_file, control);
                        %end
                    end
                end
            end
        end
    end





end















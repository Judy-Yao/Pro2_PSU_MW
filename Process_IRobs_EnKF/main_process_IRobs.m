% -------------- Set up control variables ----------------------
control = struct;
% ----Path

control.bestrack_dir = '../../raw_Obs/Bestrack/'; % directory where best-track files are
control.output_dir = '../../toEnKFobs/GOESR_IR/'; % directory where this algorithm outputs
% ---Storm information
control.storm_phase = {'MariaRI',};  
%control.storm_phase = ["Irma2ndRI",'JoseRI','MariaRI'};
control.period = {{'201709160000','201709180000'},};
%control.period = {{'201709030600','201709050600'},{'201709050600','201709070600'},{'201709160000','201709180000'}}; %YYYYMMDDHHmm
% --- WRF setup
control.nx = 297; % number of grid points along X direction
control.ny = 297; % number of grid points along Y direction
control.dx = 3; % WRF resolution: 3 km
% --- Other
control.domain_buffer = 1.5; % scaling factor
%control.filter_ratio = [36,24]; % only for testing!!
%control.filter_ratio = [6,6];
control.roi_oh = {[200,0], [30,30]}; % roi [other variables, hydrometeors]
control.obsError = [3 3];

% ---------- Loop through each storm object -------------------
for istorm = 1:length(control.storm_phase)

    % --- Make subdirectory for output
    if ~exist([control.output_dir,control.storm_phase{istorm}],'dir')
        [~, msg, ~] = mkdir(control.output_dir,control.storm_phase{istorm});
        if isempty(msg)
            disp(['Successfully created a subdirectory in ',control.output_dir,' for ',control.storm_phase{istorm}]);
        else
            error('Error: ',msg);
        end
    end

    % --- Generate hours of interest and obtain the best-track locations at these hours
    bestrack_str = Hourly_Bestrack(istorm, control); % (cell)
    
    filename = strcat(control.output_dir,control.storm_phase{istorm},'/bestrack_perHour');
    disp("Output hourly best-track location and time: "+filename);
    formatSpec = '%12s%12.3f%12.3f\n';
    fileID = fopen(filename,'w');
    for itime = 1:size(bestrack_str,1)
        fprintf(fileID, formatSpec, ...
            bestrack_str{itime,1}, bestrack_str{itime,3}(1), bestrack_str{itime,3}(2));
    end
    fclose(fileID);

   % --- Collect GOESR Obs files within the period of interest into a directory
    disp('Collecting useful GOES-R obs files for this study......');
    [] = Collect_GOESR(istorm, bestrack_str, control);
    %[Tbfile_names,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,overpass,singlepass] = Collect_MW_useful(istorm, bestrack_str, control); % ps: per swath
 

   % % --- Loop through each useful Tb file via a symbolic link
   % Tb_dir = [control.obs_collect_dir,control.storm_phase{istorm},'/*'];
   % Tb_files = strsplit(ls(Tb_dir));
   % Tb_files = Tb_files(~cellfun('isempty',Tb_files)); % get rid of annoying empty cell    

   % for iTb = 1:length(Tb_files)

    



    %end

end


































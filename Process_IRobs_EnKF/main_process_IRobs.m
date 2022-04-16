% Skeleton of the algorithm producing GOES-R obs that is readable by EnKF


% -------------- Set up control variables ----------------------
control = struct;
% ----Path
control.obs_dir = '../../raw_Obs/GOESR_IR/'; % directory where original GOESR files are
control.bestrack_dir = '../../raw_Obs/Bestrack/'; % directory where best-track files are
control.output_dir = '../../toEnKFobs/GOESR_IR/'; % directory where this algorithm outputs
control.obs_collect_dir = '../../raw_Obs/Collected_IR/'; % 
% ---Storm information
control.storm_phase = {'JoseRI',};  
%control.storm_phase = ["Irma2ndRI",'JoseRI','MariaRI'};
control.period = {{'201709050600','201709070600'},};
%control.period = {{'201709030600','201709050600'},{'201709050600','201709070600'},{'201709160000','201709180000'}}; %YYYYMMDDHHmm
% ---Satellite information
control.favCH = [8,];
control.facWL = {'6.2um', };
% --- WRF setup
control.nx = 297; % number of grid points along X direction
control.ny = 297; % number of grid points along Y direction
control.dx = 3; % WRF resolution: 3 km
% --- Other
control.domain_buffer = 1.05; % scaling factor (Enlarge the domain 3 a bit so that enough obs will be obtained)
control.filter_reso = [18;12]; 
control.search_buffer = 0.1; % degrees: lat/lon
% 18-by-18 km box with a 200-km radius of inﬂuence for non-hydro variables; 12-by-12 km box with a 30-km radius of inﬂuence for all variables
control.roi_oh = {[200,0]; [30,30]}; % roi [other variables, hydrometeors]
control.obsError = 3;
control.Sat_alt = 35000; % km


% ---------- Loop through each storm object -------------------
for istorm = 1:length(control.storm_phase)

    % --- Make subdirectory for output
    if ~exist([control.output_dir,control.storm_phase{istorm},],'dir')
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
   Collect_GOESR(istorm, bestrack_str, control);

   % --- Loop through each DAtime and output obs
   Tb_dir = [control.obs_collect_dir,control.storm_phase{istorm},'/*'];
   Tb_files = strsplit(ls(Tb_dir));
   Tb_files = Tb_files(~cellfun('isempty',Tb_files)); % get rid of annoying empty cell    

   for it = 1%:size(bestrack_str,1)
        files = [];
        for iTb = 1%:length(Tb_files)
            Tb_file = Tb_files{iTb};
            if contains(Tb_file, bestrack_str{it,1})
                files = [files, string(Tb_file)];
            else
                continue;
            end
        end
        
        DAtime = bestrack_str{it,1};
        DA_btk = bestrack_str{it,3}; % best-track location at the DA time
        disp(["Dealing with file(s) at DAtime: " + DAtime]);

        if length(files) == 1
            SingleCH_write(istorm, DAtime, DA_btk, files, control);                
        else
            MultiCh_write(istorm, DAtime, DA_btk, files, control);
        end
   end

end


































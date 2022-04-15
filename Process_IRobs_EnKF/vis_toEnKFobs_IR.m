% -------------- Set up control variables ----------------------
control = struct;
% ----Path
control.obs_dir = '../../toEnKFobs/GOESR_IR/';
control.output_dir = '../../Visual/toEnKFobs/GOESR_IR/';
% ---Storm information
control.storm_phase = 'Irma2ndRI';  
%control.storm_phase = {'Irma2ndRI','JoseRI','MariaRI'};
% --- WRF setup
control.nx = 297; % number of grid points along X direction
control.ny = 297; % number of grid points along Y direction
control.dx = 3; % WRF resolution: 3 km
% --- Other
%control.domain_buffer = 1.05; % scaling factor (Enlarge the domain 3 a bit so that enough obs will be obtained)
control.filter_reso = [18;12];
% 18-by-18 km box with a 200-km radius of inﬂuence for non-hydro variables; 12-by-12 km box with a 30-km radius of inﬂuence for all variables
control.roi_oh = {[200,0], [30,30]}; % roi [other variables, hydrometeors]

obs_dir = [control.obs_dir,control.storm_phase,'/*_so'];
bestrack_dir = [control.obs_dir,control.storm_phase,'/bestrack_perHour'];
obs_files = strsplit(ls(obs_dir));
obs_files = obs_files(~cellfun('isempty',obs_files));

disp(['Plotting Tb of storm: ',control.storm_phase]);

if ~exist([control.output_dir,control.storm_phase],'dir')
    [~, msg, ~] = mkdir(control.output_dir,control.storm_phase);
     if isempty(msg)
        disp(['Successfully created a subdirectory in ',control.output_dir,' for ',control.storm_phase]);
     else
        error('Error: ',msg);
     end
end

% -------------- Read raw best-track file ------------
fid = fopen(bestrack_dir);
bt_record = textscan(fid,'%s','delimiter','');
fclose(fid);
bt_str_all = string(bt_record{1}(:));

loc_storm = strings(length(bt_record), 3);
for ir = 1:length(bt_record)
    bt_str_per = strsplit(bt_str_all(ir));
    loc_storm(ir,:) = bt_str_per(1,:); % time,lat,lon
end


% -------------- Handle MW so file ------------
for iso = 1:length(obs_files)
    DA_time = loc_storm(iso,1);
    so_file = obs_files{iso};
    disp(so_file);
    % - Read the so file
    fid = fopen(so_file);
    obs_record = textscan(fid,'%s','delimiter','');
    fclose(fid);
    obs_str_all = string(obs_record{1}(:));
    len_record = length(obs_record{1}(:)); % how many records there are per so file

    obs_3 = zeros(len_record,6);
    for ir = 1:len_record
        obs_str_per = strsplit(obs_str_all(ir));
        obs_3(ir,:) =  str2double(obs_str_per(1,4:7)); % latitude, longitude, Tb value, ROI for hydro
    end

    % - Separate obs with large ROI from small ROI
    idx_largeROI =  obs_str_3(:,4) == 0;
    idx_smallROI = obs_str_3(:,4) == 30;
    obs = cell(length(control.roi_oh),1);
    obs{1,1} = obs_str_3(idx_largeROI,1:3);
    obs{2,1} = obs_str_3(idx_smallROI,1:3);

    lat_bt = loc_storm(iso,2); lon_bt =  loc_storm(iso,3);
    % map boundaries
    min_lat = lat_bt - (control.ny/2*control.dx)/(cos(lat_bt*(pi/180))*111);
    max_lat = lat_bt + (control.ny/2*control.dx)/(cos(lat_bt*(pi/180))*111);
    min_lon = lon_bt - (control.nx/2*control.dx)/111;
    max_lon = lon_bt + (control.nx/2*control.dx)/111;  
    % limits of the map
    min_xlat = double(min_lat-0.2);
    min_xlon = double(min_lon-0.2);
    max_xlat = double(max_lat+0.2);
    max_xlon = double(max_lon+0.2);
   
    % ---------- Plot the figure -----------
    for iroi = 1:length(control.roi_oh)
        figure;
        hFig=gcf;
        set(hFig, 'Position', [0 0 750 800]);
        % scatter Tbs on a projected map
        m_proj('mercator','lon',[min_xlon max_xlon],'lat',[min_xlat max_xlat]);    
        H = m_scatter(obs{iroi,1}(:,2), obs{iroi,1}(:,1),pointsize,obs{iroi,1}(:,3),'o','filled');
        hold on;
        m_scatter(lon_bt,lat_bt,50,0, '*'); % 0 is represented by black color in this colormap
        % use the customized colormap
        cmap = IR_colormap(0.5); caxis([185 325]); 
        cb = colorbar;
        set(cb,'Fontsize', 23);
        cb.Label.String = 'Brightness Temperature (K)';
        % add coastline 
        m_coast('color','k');
        % grid lines 
        lon_range = round(min_lon:2:max_lon);
        lat_range = round(min_lat:2:max_lat);
        m_grid('xtick',lon_range,'ytick',lat_range,'tickdir','out','fontsize',22);
        xlh = xlabel(['DA time: ', DA_time], 'Fontsize',22,'fontweight','bold');
        xlh.Position(2) = xlh.Position(2) - 0.01;  % move the label 0.02 data-units further down 
        % title
        title = [];


    end















































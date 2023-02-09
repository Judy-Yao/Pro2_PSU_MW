% This function only displays the toEnKFobs over a smaller area 

% -------------- Set up control variables ----------------------
control = struct;
% ----Path
control.obs_dir = '../../toEnKFobs/GOESR_IR/'; 
control.output_dir = '../../Visual/toEnKFobs/GOESR_IR/';
% ---Storm information
control.storm_phase = {'MARIA',};  
%control.storm_phase = {'Irma2ndRI','JoseRI','MariaRI'};
% --- WRF setup
control.nx = 20; % number of grid points along X direction
control.ny = 20; % number of grid points along Y direction
control.dx = 3; % WRF resolution: 3 km
% --- Other
control.domain_buffer = 1.05; % scaling factor (Enlarge the domain 3 a bit so that enough obs will be obtained)
control.filter_reso = [18;12]; 
% 18-by-18 km box with a 200-km radius of inﬂuence for non-hydro variables; 12-by-12 km box with a 30-km radius of inﬂuence for all variables
control.roi_oh = {[200,0], [30,30]}; % roi [other variables, hydrometeors]
% --- file of interest
control.ir = 1;

for istorm = 1:length(control.storm_phase)
    obs_dir = [control.obs_dir,control.storm_phase{istorm},'/*_so'];
    obs_files = strsplit(ls(obs_dir));
    obs_files = obs_files(~cellfun('isempty',obs_files));
    bestrack_dir = [control.obs_dir,control.storm_phase{istorm},'/bestrack_perHour'];

    disp(['Validating the pre-processing IR algorithm on strom ',control.storm_phase{istorm}]);

    % make the output directory
    if ~exist([control.output_dir,control.storm_phase{istorm}],'dir')
        [~, msg, ~] = mkdir(control.output_dir,control.storm_phase{istorm});
        if isempty(msg)
            disp(['Successfully created a subdirectory in ',control.output_dir,' for ',control.storm_phase{istorm}]);
        else
            error('Error: ',msg);
        end
    end 
    
     %-------------- Read raw best-track file ------------
     fid = fopen(bestrack_dir);
     bt_record = textscan(fid,'%s','delimiter','');
     fclose(fid);
     bt_str_all = string(bt_record{1}(:));
   
     loc_storm = strings(length(control.ir), 3); 
     bt_str_per = strsplit(bt_str_all(control.ir));
     loc_storm = bt_str_per; % time,lat,lon
     
     DA_time = loc_storm(1)
     DA_btk = str2double(loc_storm(2:3)); 

     % -------------- Handle so file ------------
     so_file = obs_files{control.ir};
     disp(so_file);

     fid = fopen(so_file);
     obs_record = textscan(fid,'%s','delimiter','');
     fclose(fid);
     obs_str_all = string(obs_record{1}(:));
     len_record = length(obs_record{1}(:)); % how many records there are per so file

     obs_str_3 = zeros(len_record,3); 
     for ir = 1:len_record
        obs_str_per = strsplit(obs_str_all(ir));
        obs_str_3(ir,:) =  str2double(obs_str_per(1,[4,5,7])); % latitude, longitude, Tb value
     end
    
     idx_largeROI =  obs_str_3(:,3) == 0;
     idx_smallROI = obs_str_3(:,3) == 30;
     obs = cell(length(control.roi_oh),1);
     obs{1,1} = obs_str_3(idx_largeROI,:);
     obs{2,1} = obs_str_3(idx_smallROI,:);

     % Geolocation for WRF domain
    nx = control.nx*control.domain_buffer; % zoom out
    ny = control.ny*control.domain_buffer; % zoom out
    % below algorithm works if only for all frequencies of interest, the DA_time are the same
    min_WRF_lat = DA_btk(1) - (ny/2*control.dx)/(cos(DA_btk(1)*(pi/180))*111);
    max_WRF_lat = DA_btk(1) + (ny/2*control.dx)/(cos(DA_btk(1)*(pi/180))*111);
    min_WRF_lon = DA_btk(2) - (nx/2*control.dx)/111;
    max_WRF_lon = DA_btk(2) + (nx/2*control.dx)/111;
    disp(['      min of xlong: ',num2str(min_WRF_lon), ', max of xlong: ',num2str(max_WRF_lon)]);
    disp(['      min of xlat: ',num2str(min_WRF_lat), ', max of xlat: ',num2str(max_WRF_lat)]);
    latitudes  = linspace(min_WRF_lat,max_WRF_lat,ny);
    longitudes = linspace(min_WRF_lon,max_WRF_lon,nx);
    [WRF_lat, WRF_lon] = meshgrid(latitudes,longitudes);
    
    % -------------- Find the masked indices for each ROI --------------------------------
    obs_idx = cell(length(control.roi_oh),1); 
    for iroi = 1:length(control.roi_oh)
        lat = obs{iroi,1}(:,1);
        lon = obs{iroi,1}(:,2);
        obs_indices{iroi,1} = (lat >= min_WRF_lat) & (lat <= max_WRF_lat) & (lon >= min_WRF_lon) & (lon <= max_WRF_lon);
    end

    % ------------------- Figure --------------------------------
    % Create a figure window    
    figure;
    set(gcf,'PaperPositionMode', 'auto'); 
    hFig=gcf;
    set(hFig, 'Position', [0 0 1200 1000]);
    m_proj('mercator','lon',[(min_WRF_lon-0.1) (max_WRF_lon+0.1)],'lat',[(min_WRF_lat-0.1) (max_WRF_lat+0.1)]);
    m_scatter(DA_btk(2),DA_btk(1),100,'black','d','filled');
    hold on;
    m_scatter(WRF_lon(:), WRF_lat(:), 50, 'black','+','linewidth',2);
    hold on;
    obs_idx = obs_indices{1,1};
    m_scatter(obs{1,1}(obs_idx,2), obs{1,1}(obs_idx,1), 50, 'red', 'filled','linewidth',2);
    hold on;
    obs_idx = obs_indices{2,1};
    m_scatter(obs{2,1}(obs_idx,2), obs{2,1}(obs_idx,1), 50, 'blue','linewidth',2);
    % grid lines 
    m_grid('off','fancy','tickdir','out','fontsize',30);

    t1 = m_text(min_WRF_lon-0.08, max_WRF_lat+0.08, '\diamondsuit Best-track Location');
    t1.Color = 'black'; t1.FontSize = 20;

    t2 = m_text(min_WRF_lon-0.08, max_WRF_lat+0.06, '+ WRF Grids ');
    t2.Color = 'black'; t2.FontSize = 20;

    t3 = m_text(max_WRF_lon - 0.2, max_WRF_lat+0.08, '\bullet Large ROI(18KM)');
    t3.Color = 'red'; t3.FontSize = 20;

    t = m_text(max_WRF_lon- 0.2, max_WRF_lat+0.05, '\circ Small ROI (12KM)');
    t.Color = 'blue'; t.FontSize = 20;

    t4 = m_text(min_WRF_lon-0.08, min_WRF_lat-0.03, 'Model Resolution: 3 KM');
    t4.Color = 'black'; t4.FontSize = 16;

    t5 = m_text(min_WRF_lon-0.08, min_WRF_lat-0.06, 'Note: pointsize does not represent any FOV-size information.');
    t5.Color = 'black'; t5.FontSize = 16;

    title([control.storm_phase{istorm},':nearest-neighbor method'],'Fontsize',25);

    saveas(gcf,[control.obs_dir,control.storm_phase{istorm},'/',control.storm_phase{istorm},'_GOES_thin_validation.png']);


end

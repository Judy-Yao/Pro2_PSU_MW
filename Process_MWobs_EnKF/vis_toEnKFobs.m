% -------------- Set up control variables ----------------------
control = struct;
% ----Path
control.obs_dir = '../../toEnKFobs/MW/';
control.output_dir = '../../Visual/toEnKFobs/MW/';
% ---Storm information
control.storm_phase = 'MariaRI';
%control.storm_phase = ["Irma2ndRI",'JoseRI','MariaRI'};
control.sensor_list = ["amsr2","atms","gmi","mhs","saphir","ssmis"];
control.ch_list= {[7,13],18,[3,13],5,5,[9,13]};
control.freq_list = {{'18.7GHzV-Pol','89GHzV-Pol'},{'183.31+-7GHzQH-Pol'},... 
         {'18.7GHzV-Pol','183.31+/-7GHzV-Pol'},{'190.31GHzV-Pol'},{'183.31+/-6.8GHz'},{'19.35GHzV-Pol','183.31+/-6.6GHzH-Pol'}};

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
len_BTrecord = length(bt_record{1}(:));

loc_storm = strings(len_BTrecord, 3); 
for ir = 1:len_BTrecord
    bt_str_per = strsplit(bt_str_all(ir));
    loc_storm(ir,:) = bt_str_per(1,:); % time,lat,lon
end

% -------------- Handle MW so file ------------
for iso = 1:length(obs_files)
    % - Read so file
    disp(obs_files{iso});
    fid = fopen(obs_files{iso});
    obs_record = textscan(fid,'%s','delimiter','');
    fclose(fid);
    obs_str_all = string(obs_record{1}(:));
    len_record = length(obs_record{1}(:)); % how many records there are per so file

    obs_str = strings(len_record/2, 5); % The first half and second half shares the same Tb record except the ROI combinaiton is different
    for ir = 1:len_record/2
        obs_str_per = strsplit(obs_str_all(ir));
        obs_str(ir,:) = obs_str_per(1,2:6); % sensor name, channel number, latitude, longitude, Tb value
    end
    % -----------------------
    chNum = unique(obs_str(:,2));
    idx_perCh = cell(length(chNum));

    [~,filename,~] = fileparts(obs_files{iso});
    file_info = split(filename,'_'); DA_time = file_info{3};

    idx_inBT = DA_time == loc_storm(:,1);

    for ich = 1:length(chNum)
        idx_perCh{ich} = obs_str(:,2) == chNum(ich);

        sensor_plf = unique(obs_str(idx_perCh{ich},1));
        chNum_char = unique(obs_str(idx_perCh{ich},2));
        chNum_int = str2double(chNum_char);
        sensor = split(sensor_plf,"_"); sensor = sensor{1};
        idx_sensor = sensor == control.sensor_list;
        idx_ch = chNum_int == control.ch_list{idx_sensor};
        freq = control.freq_list{idx_sensor}{idx_ch};

        lat = str2double(obs_str(idx_perCh{ich},3));
        lon = str2double(obs_str(idx_perCh{ich},4));
        Tb = str2double(obs_str(idx_perCh{ich},5));
        % - Plot the figure
        figure;
        hFig=gcf;
        set(hFig, 'Position', [0 0 750 800]);
        
        min_lat = min(lat);
        max_lat = max(lat);
        min_lon = min(lon);
        max_lon = max(lon);
        % limits of the map
        min_xlat = double(min_lat-0.1);
        min_xlon = double(min_lon-0.1);
        max_xlat = double(max_lat+0.1);
        max_xlon = double(max_lon+0.1);
        % scatter Tbs on a projected map
        m_proj('mercator','lon',[min_xlon max_xlon],'lat',[min_xlat max_xlat]);
        % Pointsize varies with sensor
        if contains(sensor_plf,'gmi')
            pointsize = 10;
        elseif contains(sensor_plf,'amsr2')
            pointsize = 10;
        elseif contains(sensor_plf,'ssmis')
            pointsize = 50;
        elseif contains(sensor_plf,'mhs')
            pointsize = 50;
        elseif contains(sensor_plf,'atms')
            pointsize = 50; 
        elseif contains(sensor_plf,'saphir')
            pointsize = 10;
        end
        H = m_scatter(lon,lat,pointsize,Tb,'o','filled');
        hold on;
        if contains(control.freq_list{idx_sensor}{idx_ch}, '18.7GHzV-Pol') | contains(control.freq_list{idx_sensor}{idx_ch},'19.35GHzV-Pol')
            m_scatter(str2double(loc_storm(idx_inBT,3)),str2double(loc_storm(idx_inBT,2)),50,0, '*'); % 0 is represented by black color in this colormap
        else
            m_scatter(str2double(loc_storm(idx_inBT,3)),str2double(loc_storm(idx_inBT,2)),50,299, '*'); % 299 is represented by red color in this colormap
        end
        % customize colormap
        max_T = 300;
        min_T = 80;
        min_Jet = 180;
        jetLength = (max_T) - (min_Jet) + 1;
        notJetLength = (min_Jet - 1) - (min_T + 1) + 1;
        theJet = jet(jetLength);
        jet_red = theJet(:,1);
        jet_green = theJet(:,2);
        jet_blue = theJet(:,3);
        jetAdd_red = (1.0:((jet_red(1)-1)/notJetLength)  :jet_red(1))';
        jetAdd_green = (1.0:((jet_green(1)-1)/notJetLength):jet_green(1))';
        jetAdd_blue  = (1.0:((jet_blue(1)-1)/notJetLength) :jet_blue(1))';
        cm_red   = [0.0; jetAdd_red;   jet_red;   0.0];
        cm_green = [0.0; jetAdd_green; jet_green; 0.0];
        cm_blue  = [0.0; jetAdd_blue;  jet_blue;  0.0];
        myColormap = ([cm_red cm_green cm_blue]);
        % use the customized colormap
        colormap(myColormap);
        caxis([min_T,max_T])
        cb = colorbar;
        set(cb,'Fontsize', 23);
        cb.Label.String = 'Brightness Temperature (K)';
        % add coastline 
        m_coast('color','k');
        % grid lines 
        lon_range = round(min_lon:2:max_lon);
        lat_range = round(min_lat:2:max_lat);
        m_grid('xtick',lon_range,'ytick',lat_range,'tickdir','out','fontsize',22);
        % add DA time
        xlh = xlabel(['DA time: ', DA_time], 'Fontsize',22,'fontweight','bold');
        xlh.Position(2) = xlh.Position(2) - 0.01;  % move the label 0.02 data-units further down 
        % add title
        if contains(sensor,"amsr2")
            sensor = "AMSR2\_GCOM-W1";
        elseif contains(sensor_plf,"atms")
            sensor_plf = "ATMS\_NPP";
        elseif contains(sensor_plf,"gmi_gpm")
            sensor_plf = "GMI\_GPM";
        elseif contains(sensor_plf,"mhs_metop-a")
            sensor_plf = "MHS\_METOP-A";
        elseif contains(sensor_plf,"mhs_metop-b")
            sensor_plf = "MHS\_METOP-B";
        elseif contains(sensor_plf,"mhs_n18")
            sensor_plf = "MHS\_N18";
        elseif contains(sensor_plf,"mhs_n19")
            sensor_plf = "MHS\_N19";
        elseif contains(sensor_plf,"saphir")
            sensor_plf = "SAPHIR\_MEGHAT";
        elseif contains(sensor_plf,"ssmis_f16")
            sensor_plf = "SSMIS\_F16";
        elseif contains(sensor_plf,"ssmis_f17")
            sensor_plf = "SSMIS\_F17";
        elseif contains(sensor_plf,"ssmis_f18")
            sensor_plf = "SSMIS\_F18";
        end
        title_char1 = control.storm_phase + sensor_plf ;
        title_char2 = 'Channel ' + chNum_char + ": "+freq ;
        title_char = [title_char1, title_char2];
        title(title_char,'Fontsize',20);

        save_dir = [control.output_dir+string(control.storm_phase)+'/'+string(DA_time)+'_ch'+chNum_char+'.png'];
        saveas(gcf,save_dir);
    end
    
    

end

% =============================================================================================================================
% Script to plot MW observations that are processed by MOPS 
% =============================================================================================================================

% -------------- Set up control variables ----------------------
control = struct;
% ----Path
control.obs_dir = '../../toEnKFobs/MW/';
control.output_dir = '../../Visual/toEnKFobs/MW/';
% ---Storm information
control.storm_phase = 'MariaRI';
control.sensor_list = ["amsr2","atms","gmi","mhs","saphir","ssmi","ssmis"];
control.freq_list = {{'18.7GHzV-Pol','89GHzV-PolA-Scan','89GHzV-PolB-Scan'},...
        {'183.31+-7GHzQH-Pol'},{'18.7GHzV-Pol','183.31+/-7GHzV-Pol'},{'190.31GHzV-Pol'},{'183.31+/-6.8GHz'},{'fcdr\_tb19v','fcdr\_tb85v'},{'19.35GHzV-Pol','183.31+/-6.6GHzH-Pol'}};
control.ch_list= {[7,13],18,[3,13],5,5,[1,6],[13,9]}; % Reference to SourceCode/output_sensor_facts.py
control.numPerscan = {[243,486],96,[221,221],90,182,[64,128],[90,180]}; % Reference to SourceCode/output_sensor_facts.py (be careful with AMSR2 89GHz-Pol!)
control.nx = 297;
control.ny = 297;
control.dx = 3;

control.roi_oh = {[200,0]; [60,60]};
% --- Customize colormap
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

% -------------- Obtain files ----------------------
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


% -------------- Read and process raw best-track file ------------
fid = fopen(bestrack_dir);
bt_record = textscan(fid,'%s','delimiter','');
fclose(fid);
bt_str_all = string(bt_record{1}(:));
[bt_unique_all,idx] = unique(bt_str_all(:,1)); % discards repeated times
len_unique_BTrecord = length(bt_unique_all);
bt_str_all = bt_str_all(idx);

loc_storm = strings(len_unique_BTrecord, 3); 
for ir = 1:len_unique_BTrecord
    bt_str_per = strsplit(bt_str_all(ir));
    loc_storm(ir,:) = bt_str_per(1,:); % time,lat,lon
end

% -------------- Handle MW so file ------------
for iso = 1:length(obs_files)
    so_file = obs_files{iso};
    disp(so_file);
    % --- Read the so file
    fid = fopen(so_file);
    obs_record = textscan(fid,'%s','delimiter','');
    fclose(fid);
    % --- Separate records for different ROI plans
    obs_str_all = string(obs_record{1}(:));
    len_record = length(obs_record{1}(:)); % how many records there are per so file
    
    obs_split_all = strings(len_record, 17);
    for ir = 1:len_record
        obs_split_all(ir,:) = strsplit(obs_str_all(ir));
    end

    index_ROI = cell(length(control.roi_oh),1);
    index_ROI{1,1} = obs_split_all(:,7) == repmat("0",len_record,1); % mask out the records using first ROI plan
    index_ROI{2,1} = obs_split_all(:,7) == repmat("60",len_record,1); % mask out the records using second ROI plan
    % --- For different ROI plans, read values of variables used for plotting
    obs_str = cell(length(control.roi_oh),1);
    for iroi = 1:length(index_ROI)
        obs_str{iroi} = strings(sum(index_ROI{iroi} == 1),7);
    end

    for iroi = 1:length(index_ROI)
         obs_str{iroi}(:,1:5) = obs_split_all(index_ROI{iroi},2:6); 
         obs_str{iroi}(:,6:7) = obs_split_all(index_ROI{iroi},10:11);
    end
    % --- Identify location of the storm at a DA time
    [~,filename,~] = fileparts(so_file);
    file_info = split(filename,'_'); DA_time = file_info{3};

    idx_inBT = DA_time == loc_storm(:,1);
    lat_bt = str2double(loc_storm(idx_inBT,2));
    lon_bt = str2double(loc_storm(idx_inBT,3));

    % --- Allocate indices for records for each sensor and each channel
    for iroi = 1:length(index_ROI)
        sensors_plf = unique(obs_str{iroi,1}(:,1)); % unique sensor names
        
        chNum = cell(length(sensors_plf),1); % channel number
        idx_perChperSS = cell(length(sensors_plf),1);
        for iss = 1:length(sensors_plf)
            idx_perSS = sensors_plf(iss) == obs_str{iroi}(:,1);
            chNum{iss,1} = unique(obs_str{iroi,1}(idx_perSS,2));
            idx_perChperSS{iss,1} = cell(length(chNum{iss,1}),1);
        end
    
        % --- Plot each kind of record
        for iss = 1:length(sensors_plf)
            for ich = 1:length(idx_perChperSS{iss,1})
                sensor_plf = sensors_plf(iss);
                sensor = split(sensor_plf,"_"); sensor = sensor{1};

                chNum_char = chNum{iss,1}(ich);
                chNum_int = str2double(chNum_char);
                iss_list = sensor == control.sensor_list;
                ich_list = chNum_int == control.ch_list{iss_list};
                freq = control.freq_list{iss_list}{ich_list};

                idx_ss_obs_str = obs_str{iroi,1}(:,1) == sensor_plf;
                ichNum_obs_str = obs_str{iroi,1}(:,2) == chNum_char;
                idx_perChperSS{iss,1}{ich,1} = idx_ss_obs_str & ichNum_obs_str;
             
                lat = str2double(obs_str{iroi,1}(idx_perChperSS{iss,1}{ich,1},3));
                lon = str2double(obs_str{iroi,1}(idx_perChperSS{iss,1}{ich,1},4));
                Tb = str2double(obs_str{iroi,1}(idx_perChperSS{iss,1}{ich,1},5));
                % - Plot the figure
                figure;
                hFig=gcf;
                set(hFig, 'Position', [0 0 750 800]);
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
                % scatter Tbs on a projected map
                m_proj('mercator','lon',[min_xlon max_xlon],'lat',[min_xlat max_xlat]);
                % Pointsize is inverse proportional to number per scan for each channel
                pointsize = control.numPerscan{iss_list}(ich_list)*(-0.114)+60.26; % linear interpolation
                % number per scan: 486  <--> pointsize = 5
                % number per scan: 90   <--> pointsize = 50
                H = m_scatter(lon,lat,pointsize,Tb,'o','filled');
                hold on;
                if contains(control.freq_list{iss_list}{ich_list}, '18.7GHzV-Pol') | contains(control.freq_list{iss_list}{ich_list},'19.35GHzV-Pol') | contains(control.freq_list{iss_list}{ich_list},'fcdr_tb19v')
                    m_scatter(lon_bt,lat_bt,50,0, '*'); % 0 is represented by black color in this colormap
                else
                    m_scatter(lon_bt,lat_bt,50,299, '*'); % 299 is represented by red color in this colormap
                end
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
                    sensor_plf = "AMSR2\_GCOM-W1";
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
                elseif contains(sensor_plf,"ssmi_f15")
                    sensor_plf = "SSMI\_F15"; 
                elseif contains(sensor_plf,"ssmis_f16")
                    sensor_plf = "SSMIS\_F16";
                elseif contains(sensor_plf,"ssmis_f17")
                    sensor_plf = "SSMIS\_F17";
                elseif contains(sensor_plf,"ssmis_f18")
                    sensor_plf = "SSMIS\_F18";
                end
                title_char1 = control.storm_phase + " " + sensor_plf; 
                title_char2 = "Ch" + chNum_char + ": "+freq ;
                title_char = [title_char1, title_char2];
                title(title_char,'Fontsize',20)
         
                if iroi == 1
                    if contains(control.freq_list{iss_list}{ich_list}, '18.7GHzV-Pol') | contains(control.freq_list{iss_list}{ich_list},'19.35GHzV-Pol') | contains(control.freq_list{iss_list}{ich_list},'fcdr\_tb19v')
                        save_dir = [control.output_dir+string(control.storm_phase)+'/BigROI/lf/'+string(DA_time)+'_ch'+chNum_char+'.png'];  
                    else
                        save_dir = [control.output_dir+string(control.storm_phase)+'/BigROI/hf/'+string(DA_time)+'_ch'+chNum_char+'.png'];
                    end
                else
                    if contains(control.freq_list{iss_list}{ich_list}, '18.7GHzV-Pol') | contains(control.freq_list{iss_list}{ich_list},'19.35GHzV-Pol') | contains(control.freq_list{iss_list}{ich_list},'fcdr\_tb19v')
                        save_dir = [control.output_dir+string(control.storm_phase)+'/SmallROI/lf/'+string(DA_time)+'_ch'+chNum_char+'.png'];
                    else
                        save_dir = [control.output_dir+string(control.storm_phase)+'/SmallROI/hf/'+string(DA_time)+'_ch'+chNum_char+'.png'];
                    end
                end
               saveas(gcf, save_dir); 
            end
        end
    end    
end

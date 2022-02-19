function [] = MW_vis(storm_phase, platform, sensor, min_lat, max_lat, min_lon, max_lon, Swath_used, if_swath_good, ChIdx_perSwath, ChName_perSwath, Tb_file, control)
% plot MW Tbs

% -------- Loop through each swath -----
for i_sw = 1:length(if_swath_good)
    if if_swath_good(i_sw) == false
        continue;
    else
        % ----- Read Variables --------
        % -time
        fileheader = split(h5readatt(Tb_file,'/','FileHeader'),';');
        for i = 1:length(fileheader)
            if contains(fileheader{i},'StartGranuleDateTime') 
                time_start = erase(fileheader{i},'StartGranuleDateTime=');
            elseif contains(fileheader{i},'StopGranuleDateTime') 
                time_end = erase(fileheader{i},'StopGranuleDateTime=');
            end
        end
        % -Tbs of a certain channel
        Tb_char = [Swath_used{i_sw},'/Tc'];
        Tb_Chs = h5read(Tb_file,Tb_char);
        idx_Tb_oneCh = ChIdx_perSwath(i_sw);
        Tb_oneCh = squeeze(Tb_Chs(idx_Tb_oneCh,:,:));
        
        % -lat/lon for a certian swath
        lat_char = [Swath_used{i_sw},'/Latitude'];
        xlat = h5read(Tb_file,lat_char);
        lon_char = [Swath_used{i_sw},'/Longitude'];
        xlon = h5read(Tb_file,lon_char);

        % ---------- Operations on Variables -----------------
        % time
        time_str1 = erase(extractBefore(strtrim(time_start),17),'-');
        time_str2 = erase(extractBetween(strtrim(time_end),6,16),'-');
        time = append(time_str1, '-', time_str2); 

        time_save = replace(erase(time,':'),'-','_');
        % Masked Tbs, lat, lon based on the area of interest
        idx_mask = (xlon < max_lon) & (xlon > min_lon) & (xlat < max_lat) & (xlat > min_lat);
        lat_need = xlat(idx_mask);
        lon_need = xlon(idx_mask);
        Tb_need = Tb_oneCh(idx_mask);
        
        % ------------ Plotting --------------------
        figure;
        set(gcf,'PaperPositionMode', 'auto'); 
        hFig=gcf;
        set(hFig, 'Position', [0 0 900 950]);
        set(gcf,'PaperPositionMode', 'auto');

        % limits of the map
        min_xlat = double(min_lat-0.2);
        min_xlon = double(min_lon-0.2);
        max_xlat = double(max_lat+0.2);
        max_xlon = double(max_lon+0.2);
        % scatter Tbs on a projected map
        m_proj('mercator','lon',[min_xlon max_xlon],'lat',[min_xlat max_xlat]);
        % Pointsize varies with sensor
        if contains(sensor,'GMI')
            pointsize = 10;
        elseif contains(sensor,'AMSR2')
            pointsize = 10;
        elseif contains(sensor,'SSMIS')
            pointsize = 50;
        elseif contains(sensor,'MHS')
            pointsize = 50;
        elseif contains(sensor,'ATMS')
            pointsize = 50; 
        elseif contains(sensor,'SAPHIR')
            pointsize = 10;
        end
        H = m_scatter(lon_need,lat_need,pointsize,Tb_need,'o','filled');
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
        m_grid('off','fancy','tickdir','out','fontsize',22);
        % title
        disp(ChName_perSwath(i_sw));
        title_char1 = cellstr(append(storm_phase,': ',platform,'-',sensor));
        title_char = [title_char1,ChName_perSwath(i_sw),time];
        title(title_char,'Fontsize',20);
        if contains(ChName_perSwath(i_sw),'18.7GHzV-Pol') | contains(ChName_perSwath(i_sw),'19.35GHzV-Pol') 
            disp([control.plot_dir,'/',storm_phase,'/',sensor,'/low_f/','MW_lf_', char(time_save),'.png']);
            saveas(gcf, [control.plot_dir,'/',storm_phase,'/',sensor,'/low_f/','MW_lf_', char(time_save),'.png']);
        else
            disp([control.plot_dir,'/',storm_phase,'/',sensor,'/high_f/','MW_hf_', char(time_save),'.png']);
            saveas(gcf, [control.plot_dir,'/',storm_phase,'/',sensor,'/high_f/','MW_hf_', char(time_save),'.png']);
        end
     end
end

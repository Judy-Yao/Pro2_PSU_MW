function [] = Plot_AMSR2scans(storm_phase, platform, sensor, min_lat_dpy, max_lat_dpy, min_lon_dpy, max_lon_dpy, time, loc_storm, lat_need, lon_need, Tb_need, DAtime, control)

        % ------------ Plotting --------------------
        figure;
        set(gcf,'PaperPositionMode', 'auto');
        hFig=gcf;
        set(hFig, 'Position', [0 0 800 850]);
        set(gcf,'PaperPositionMode', 'auto');

        % limits of the map
        min_xlat_dpy = double(min_lat_dpy-0.2);
        min_xlon_dpy = double(min_lon_dpy-0.2);
        max_xlat_dpy = double(max_lat_dpy+0.2);
        max_xlon_dpy = double(max_lon_dpy+0.2);
        % scatter Tbs on a projected map
        m_proj('mercator','lon',[min_xlon_dpy max_xlon_dpy],'lat',[min_xlat_dpy max_xlat_dpy]);
        % Pointsize varies with sensor
        pointsize = 10;        
        % scatter points
        m_scatter(lon_need,lat_need,pointsize,Tb_need,'o','filled');
        hold on;
        % indicate the best_track location of storm
        m_scatter(loc_storm(2),loc_storm(1),50,299, '*'); % 299 is represented by red color in this colormap
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
        m_coast('color','k','linewidth',2);
        m_grid('off','fancy','tickdir','out','fontsize',22);
        % title
        xlh = xlabel(['DA time: ', DAtime], 'Fontsize',22,'fontweight','bold');        
        xlh.Position(2) = xlh.Position(2) - 0.02;  % move the label 0.02 data-units further down

        %disp('89GHzV-pol');
        title_char1 = cellstr(append(storm_phase,': ',platform,'-',sensor));
        title_char = [title_char1, '89GHzV-pol', time];
        title(title_char,'Fontsize',20);
        %subtitle(['DA time: ', DAtime]);
        
        time_save = replace(erase(time,':'),'-','_');
        disp([control.plot_dir,'/',storm_phase,'/hf/','MW_hf_', char(time_save),'.png']);
        saveas(gcf, [control.plot_dir,'/',storm_phase,'/hf/','MW_hf_', char(time_save),'.png']);


end

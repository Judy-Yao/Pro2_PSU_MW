function [] = IR_vis(istorm, min_lat, max_lat, min_lon, max_lon, Tb_file, control)

    % ----------------- Read Variables -----------------------------
    % time
    time_start = ncreadatt(Tb_file,'/','time_coverage_start');
    time_end =  ncreadatt(Tb_file,'/','time_coverage_end');
    % Tbs of Ch 8
    Tbs = ncread(Tb_file, 'CMI');
    % Vairables about lat/lon
    x=ncread(Tb_file,'x');
    y=ncread(Tb_file,'y');

    % ----------------- Operations on Variables ------------------------------
    % time
    time_str1 = extractBefore(time_start,17);
    time_str2 = extractBetween(time_end,12,16);
    time = append(time_str1, '-', time_str2);
     
    time_save_str1 = extractBefore(time_start,5);
    time_save_str2 = extractBetween(time_start,6,7);
    time_save_str3 = extractBetween(time_start,9,10);
    time_save_str4 = extractBetween(time_start,12,13);
    time_save_str5 = extractBetween(time_start,15,16);
    time_save = append(time_save_str1,time_save_str2,time_save_str3,time_save_str4,time_save_str5);
    % Lat/lon;
    req=6378137;
    rpol=6356752.31414;
    H=42164160;
    l0=ncread(Tb_file,'nominal_satellite_subpoint_lon')*pi/180;
    [y2d,x2d]=meshgrid(y,x);
    a=(sin(x2d)).^2+(cos(x2d)).^2.*((cos(y2d)).^2+(req/rpol.*sin(y2d)).^2);
    b=-2.*H.*cos(x2d).*cos(y2d);
    c=H.^2-req.^2;
    r=(-b-sqrt(b.^2-4.*a.*c))/2./a;
    sx=r.*cos(x2d).*cos(y2d);
    sy=-r.*sin(x2d);
    sz=r.*cos(x2d).*sin(y2d);
    lat=atan(req.^2./rpol.^2.*sz./sqrt((H-sx).^2+sy.^2));
    lon=l0-atan(sy./(H-sx));
    lat=lat*180/pi;
    lon=lon*180/pi;
    lat(find(imag(lat)))=NaN;
    lon(find(imag(lon)))=NaN;
    % -------------- Mask ------------------------------
    % Get rid of points with Nan values
    idx_noNan = ~isnan(lat);
    lat_noNan = lat(idx_noNan);
    lon_noNan = lon(idx_noNan);
    Tbs_needed = Tbs(idx_noNan);
    % Get rid of points outside the wanted area
    idx_mask_geo = (lon_noNan < max_lon) & (lon_noNan > min_lon) & (lat_noNan < max_lat) & (lat_noNan > min_lat);
    lat_needed = lat_noNan(idx_mask_geo);
    lon_needed = lon_noNan(idx_mask_geo);
    Tbs_needed = Tbs_needed(idx_mask_geo);

    %-------------- Plot Figures ----------------------
    % Create a figure window    
    figure;
    set(gcf,'PaperPositionMode', 'auto'); 
    hFig=gcf;
    set(hFig, 'Position', [0 0 1000 1000]);
    set(gcf,'PaperPositionMode', 'auto');

    % Limits of the map
    min_lat_plt = min(lat_needed)-0.3;
    min_lon_plt = min(lon_needed)-0.3;
    max_lat_plt = max(lat_needed)+0.3;
    max_lon_plt = max(lon_needed)+0.3;

    % Initialize the projection 
    m_proj('mercator','lon',[min_lon_plt max_lon_plt],'lat',[min_lat_plt max_lat_plt]);
    % Plot Tbs distribution by m_scatter
    H = m_scatter(lon_needed,lat_needed,Tbs_needed, Tbs_needed, 'o','filled');
    
    % Set up my own colormap
    max_T = 250;
    min_T = 190;
    min_Jet = 210;
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

    colormap(myColormap);
    caxis([min_T,max_T])
    cb = colorbar;
    set(cb,'Fontsize', 23);
    cb.Label.String = 'Brightness Temperature (K)';

    % Coast line
    m_coast('color','k');
    % Grid lines  
    m_grid('off','fancy','tickdir','out','fontsize',22);

    % title
    if string(control.storm_phase{istorm}) == 'Irma_1stRI'
        title(['CH8 1st RI of Irma: ', time],'Fontsize',20);
    elseif string(control.storm_phase{istorm}) == 'Irma_2ndRI'
        title(['CH8 2nd RI of Irma: ', time],'Fontsize',20);
    elseif string(control.storm_phase{istorm}) == 'Jose_RI'
        title(['CH8 RI of Jose: ', time],'Fontsize',20);
    elseif string(control.storm_phase{istorm}) == 'Maria_RI'
        title(['CH8 RI of Maria: ', time],'Fontsize',20);
    end
    savename = append(control.plot_dir,control.storm_phase{istorm},'/','IR_CH8_',time_save,'.png');
    saveas(gcf, string(savename));

end




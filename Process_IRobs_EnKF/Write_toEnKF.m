function [] = Write_toEnKF(Tb_file,control)

    % ----------------Reading GOES DATA----------------------------------
    % Tbs of Ch 8
    Tbs = ncread(Tb_file, 'CMI');
    % Vairables about lat/lon
    x=ncread(Tb_file,'x');
    y=ncread(Tb_file,'y');
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
    % Get rid of points with Nan values
    idx_noNan = ~isnan(lat);
    lat_noNan = lat(idx_noNan);
    lon_noNan = lon(idx_noNan);
    %Tbs_needed = Tbs(idx_noNan);
    
    % ------------------- GRID INFO --------------------------------
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
    
    % -------------------  --------------------------------
    idx_d3 = (lat_noNan >= min_WRF_lat) & (lat_noNan <= max_WRF_lat) & (lon_noNan >= min_WRF_lon) & (lon_noNan <= max_WRF_lon);
    lat_d3 = lat_noNan(idx_d3); lon_d3 = lon_noNan(idx_d3); 
    
    filter_grid_step = control.filter_reso / control.dx; 
    grid_start(2) = floor(filter_grid_step(2) / 2);
    grid_start(1) = grid_start(2) + .5*(2*filter_grid_step(2) - filter_grid_step(1));
    
    slots_x = cell(size(control.roi_oh));
    slots_y = cell(size(control.roi_oh));
    obs_index = cell(size(control.roi_oh));
    
    for iroi = 1:length(control.roi_oh) % ROI first references [200,0] second [30,30]
        obs_index{iroi} = grid_start(iroi):filter_grid_step(iroi):length(lat_d3);   
    end






end

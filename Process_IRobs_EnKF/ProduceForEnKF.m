% nearest-neighbor method

function [myTimes,mySat_name,myChNum,myLat,myLon,myTb,myROI_hydro,myROI_other,myObsErr,mySat_alt] = ProduceForEnKF(DAtime, DA_btk, Ch, Tb_file, control)

    % ---------------------------------------------------------------------
    % ------------- Read and sligtly re-process lat, lon, Tb from one GOESR file ----------
    % ---------------------------------------------------------------------
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
    lat_col = lat(idx_noNan);
    lon_col = lon(idx_noNan);
    Tb_col = Tbs(idx_noNan);  
    
    % ---------------------------------------------------------------------
    % ---- Define area of interest for WRF simulation AND
    % ---- Select the raw obs for every grid point for EnKF assimilation
    % ---------------------------------------------------------------------
    % Prepare area: best-track location followed
    nx = control.nx*control.domain_buffer; % zoom out
    ny = control.ny*control.domain_buffer; % zoom out
    % Note: below algorithm works if only for all channels of interest, the DA_time are the same
    min_XLAT = DA_btk(1) - (ny/2*control.dx)/(cos(DA_btk(1)*(pi/180))*111);
    max_XLAT = DA_btk(1) + (ny/2*control.dx)/(cos(DA_btk(1)*(pi/180))*111);
    min_XLONG = DA_btk(2) - (nx/2*control.dx)/111;
    max_XLONG = DA_btk(2) + (nx/2*control.dx)/111;
    disp(['    min of xlong: ',num2str(min_XLONG), ', max of xlong: ',num2str(max_XLONG)]);
    disp(['    min of xlat: ',num2str(min_XLAT), ', max of xlat: ',num2str(max_XLAT)]);
    latitudes  = linspace(min_XLAT,max_XLAT,ny);
    longitudes = linspace(min_XLONG,max_XLONG,nx);
    [XLAT, XLONG] = meshgrid(latitudes,longitudes);
    
    % Separate gird points for ROI plan 1 from ROI plan 2
    % Original grid points:         Filtered grid points for one ROI:
    %       * * * * *                *   *   * 
    %       * * * * *                          
    %       * * * * *                *   *   * 
    %       * * * * *                          
    %       * * * * *                *   *   * 
    slots_x = cell(size(control.roi_oh));
    slots_y = cell(size(control.roi_oh));
    obs_index = cell(size(control.roi_oh));      

    % Locate the start point for each ROI plan
    filter_grid_step = control.filter_reso / control.dx;
    grid_start(2) = floor(filter_grid_step(2) / 2); % start point for ROI plan 1 
    grid_start(1) = grid_start(2) + .5*(2*filter_grid_step(2) - filter_grid_step(1)); % start point for ROI plan 2 
   
    % loop through each ROI plan to filter grid points and obtain GOESR obs for filtered points
    for iroi = 1:length(control.roi_oh)
        slots_x{iroi} = grid_start(iroi):filter_grid_step(iroi):nx;
        slots_y{iroi} = grid_start(iroi):filter_grid_step(iroi):ny;
        obs_index{iroi} = PickRawforCRTM(lat_col,lon_col,Tb_col,min_XLONG,max_XLONG,min_XLAT,max_XLAT,latitudes,longitudes,slots_x{iroi},slots_y{iroi},control);
    end
 
    % ---------- Process obs for each ROI plan --------------
    % Preallocating memory
    myTimes = cell(size(control.roi_oh));
    mySat_name = cell(size(control.roi_oh));
    myChNum = cell(size(control.roi_oh));
    myLat = cell(size(control.roi_oh));
    myLon = cell(size(control.roi_oh));
    myTb  = cell(size(control.roi_oh));
    myROI_hydro = cell(size(control.roi_oh));
    myROI_other = cell(size(control.roi_oh));
    myObsErr = cell(size(control.roi_oh));
    mySat_alt = cell(size(control.roi_oh));

    for iroi = 1:length(control.roi_oh) % ROI first references [200,0] second [30,30]
        obs_index_array = obs_index{iroi};
        obs_index_1d = obs_index_array(obs_index_array(:) == obs_index_array(:)); % get rid of obs_index with value NaN
        lat = lat_col(obs_index_1d);
        lon = lon_col(obs_index_1d);
        Tb = Tb_col(obs_index_1d);
        % randomize IR records for this ROI
        randOrder = randperm(length(cat(1,Tb(:))));
        
        myTimes{iroi} = repmat([DAtime+"00"], length(Tb), 1);
        mySat_name{iroi} = repmat("abi_gr", length(Tb), 1);
        myChNum{iroi} = repmat(Ch, length(Tb), 1);  
        myLat{iroi} = lat(randOrder); myLon{iroi} = lon(randOrder); myTb{iroi} = Tb(randOrder);
        myROI_hydro{iroi} = repmat(control.roi_oh{iroi}(2), length(Tb), 1); 
        myROI_other{iroi} = repmat(control.roi_oh{iroi}(1), length(Tb), 1);
        myObsErr{iroi} =  repmat(control.obsError, length(Tb), 1);   
        mySat_alt{iroi} = repmat(control.Sat_alt, length(Tb), 1);    
    end
    











end

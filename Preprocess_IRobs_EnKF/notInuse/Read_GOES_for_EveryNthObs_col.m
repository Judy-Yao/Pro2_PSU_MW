% Read and sligtly re-process lat, lon, Tb from one GOESR file

function [lat_noNan,lon_noNan,Tb_noNan] = Read_GOES(Tb_file)

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
    Tb_noNan = Tbs(idx_noNan);

end

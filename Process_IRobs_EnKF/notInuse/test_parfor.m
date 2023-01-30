
Tb_file = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/raw_Obs/GOESR_IR/MARIA/OR_ABI-L2-CMIPF-M3C08_G16_s20172570000354_e20172570011121_c20172570011191.nc';
geo_file = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/Preprocess_Domain/MARIA/201709140000/geo_em.d03.nc';
roi_oh = {[200,0]; [30,30]};
filter_reso = [18;12];
dx = 3;

% read values from observation files
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

xlat_m = ncread(geo_file,'XLAT_M');
xlon_m = ncread(geo_file,'XLONG_M');
min_xlat = min(xlat_m,[],'all')-0.05;
max_xlat = max(xlat_m,[],'all')+0.05;
min_xlon = min(xlon_m,[],'all')-0.05;
max_xlon = max(xlon_m,[],'all')+0.05;


idx_inArea = (lat_col >= min_xlat) & (lat_col <= max_xlat) & (lon_col >= min_xlon) & (lon_col <= max_xlon);
lat_obs = lat_col(idx_inArea);
lon_obs = lon_col(idx_inArea);
Tb_obs = Tb_col(idx_inArea);
obs_index = cell(size(roi_oh));

% Locate the start point for each ROI plan
filter_grid_step = filter_reso / dx;
grid_start(2) = floor(filter_grid_step(2) / 2); % start point for ROI plan 1 
grid_start(1) = grid_start(2) + .5*(2*filter_grid_step(2) - filter_grid_step(1));

tic
parpool(48);
for iroi = 1:length(roi_oh)
    obs_index{iroi} = Pick(lat_obs,lon_obs,xlat_m,xlon_m,grid_start(iroi),filter_grid_step(iroi));
    [v, w] = unique(obs_index{iroi},'stable' );
    duplicate_idx = setdiff(1:numel(obs_index{iroi}),w );
    disp(['Number of repeated obs for the same ROI selection: ', num2str(length(duplicate_idx))]);
end

% find the repeated data
[Cdata] = intersect(obs_index{1},obs_index{2});
disp(['Number of repeated obs across two ROI selections: ', num2str(length(Cdata))]);
toc

function [idx_getObs] = Pick(lat_obs,lon_obs,xlat_m,xlon_m,start_grid,step_grid)

    idx_step = start_grid:step_grid:size(xlat_m,1);  % xlat_m/xlon_m has same dimension value along x and y axis 
    xlat = xlat_m(idx_step, idx_step);
    xlon = xlon_m(idx_step, idx_step);
    xlat_col = xlat(:);
    xlon_col = xlon(:);

    % For each filtered WRF grid, select the nearest obs 
    idx_getObs = nan(length(xlat_col),1);
   
    for id =1:length(xlon_col)
        dis_ig = distance(xlat_col(id),xlon_col(id),lat_obs,lon_obs);
        idx = find( dis_ig == min(dis_ig) );
        idx_getObs(id,1) = idx;
    end
    
end





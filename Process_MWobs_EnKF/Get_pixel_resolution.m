function [crossTrack, alongTrack] = Get_pixel_resolution(scantype,ch_num,fov_alongTrack,fov_crossTrack,lat,lon,zenith,sat_lat,sat_lon,sat_alt)

Re = 6378.1;

if (strcmp(scantype,'Conical'))
    
    crossTrack = ones(length(lat),1,'double')*double(fov_crossTrack);
    alongTrack = ones(length(lat),1,'double')*double(fov_alongTrack);
    
elseif (strcmp(scantype,'Cross-track'))
    
    crossTrack_degrees = ones(length(lat),1,'double')*fov_crossTrack;
    alongTrack_degrees = ones(length(lat),1,'double')*fov_alongTrack;
        
    [xEast,yNorth,zUp] = geodetic2enu(sat_lat,sat_lon,1000*sat_alt,lat,lon,0,referenceEllipsoid('WGS 84'));
    slant_distance = sqrt(xEast.^2 + yNorth.^2 + zUp.^2) / 1e3;

    alongTrack = 2.*slant_distance.*tand(alongTrack_degrees / 2);
    
    crossTrack_gifov_notProj = 2.*slant_distance.*tand(crossTrack_degrees / 2);
    crossTrack = crossTrack_gifov_notProj .* sind(90 + crossTrack_degrees/2) ./ sind(90 - (zenith + crossTrack_degrees/2));
end

% [xEast,yNorth,zUp] = geodetic2enu(lat,lon,h,lat0,lon0,h0,spheroid) transforms the geodetic coordinates specified by lat, lon, and h to the local east-north-up (ENU) Cartesian coordinates specified by xEast, yNorth, and zDown. Specify the origin of the local ENU system with the geodetic coordinates lat0, lon0, and h0. 

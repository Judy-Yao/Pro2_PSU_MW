function [if_swath_good] = if_area_covered(min_lat, max_lat, min_lon, max_lon, Swath_used, Tb_file, control)
% Determine if this observation file should be used or not based on if it captures the area of interest
% Number of scan lat/lon to judge if the area is covered

if_swath_good = [];

for i_sw = 1:length(Swath_used)
    scan_lat_char = [Swath_used{i_sw},'/SCstatus/SClatitude'];
    scan_lon_char = [Swath_used{i_sw},'/SCstatus/SClongitude']; 
    
    scan_lat = h5read(Tb_file,scan_lat_char);
    scan_lon = h5read(Tb_file,scan_lon_char);
    % determine if this area is covered enough
    num_useful_scan = sum((scan_lon < max_lon) & (scan_lon > min_lon) & (scan_lat < max_lat) & (scan_lat > min_lat));

    if num_useful_scan >= 20
        if_swath_good = [if_swath_good, true];
    else
        if_swath_good = [if_swath_good, false];
    end


end
    
% sanity check
if length(Swath_used) ~= length(if_swath_good)
    disp('Error determining if swath covers the area of interest');
end

end

% This function picks a L1C MW obs for each grid point in a separated grid mesh for a frequency
function [obs_index] = PickRawforCRTM(lat_obs,lon_obs,Tb_obs,min_XLONG,max_XLONG,min_XLAT,max_XLAT,xlon_1d,xlat_1d,index_x,index_y,control)

    lat_obs_col = reshape(lat_obs,[],1);
    lon_obs_col = reshape(lon_obs,[],1);
    Tb_obs_col = reshape(Tb_obs,[],1);

    obs_distances = nan(length(index_x),length(index_y));
    obs_index = nan(length(index_x),length(index_y));
    obs_n_candidates = nan(length(index_x),length(index_y));
    obs_was_pushed_from_top_candidate = zeros(size(obs_n_candidates));

    n_points_noObs = 0;

    % ------------- Select several possible raw-obs candidates for each grid point (x,y) -------------
    % loop through grids of interest along xlon_1d
    for x = index_x 
        % find the index of the grid in the subset of grids of interest along xlon_1d
        idx_subset_x = find(index_x==x);
        % find the lon of the focal point
        focus_lon = xlon_1d(x);
        % search obs in the range by lontitude bounds & flag obs that are in the range
        lg_x = (lon_obs_col > max(min_XLONG, focus_lon - control.search_buffer)) & (lon_obs_col < min(max_XLONG, focus_lon + control.search_buffer));
        
        % loop through grids of interest along xlat_1d 
        for y = index_y
            clear Distance_obs_inArea
            % find the index of the grid in the subset of grids of interest along xlat_1d
            idx_subset_y = find(index_y==y);
             % find the lat of the focal point
            focus_lat = xlat_1d(y);
            % search obs in the range by latitude bounds & flag obs that are in the range
            lg_y = (lat_obs_col > max(min_XLAT, focus_lat - control.search_buffer*1.5)) & (lat_obs_col <  min(max_XLAT, focus_lat + control.search_buffer*1.5));

            % flag obs that meet both requirements
            lg_good = ((lg_x + lg_y) == 2);
            if max(lg_good) == 0
                continue;
            end
            % find the index of obs that meet both requirements
            idx_obs_good = find(lg_good);
            % find the lon/lat/Tb of obs that meet both requirements
            lat_obs_good = lat_obs_col(idx_obs_good);
            lon_obs_good = lon_obs_col(idx_obs_good);
            % For each possible good obs, calculate its distance from the focal point
            for i = 1:length(idx_obs_good)
                Distance_obs_inArea(i) = distance(lat_obs_good(i),lon_obs_good(i),focus_lat,focus_lon);
            end
            
            % ------------- Identify the best raw-obs candidate -------------
            if sum(idx_obs_good) > 0
                my_obs_index = NaN;
                % --- Keep looping the raw-obs candidates for the ith grid in (idx_subset_x,idx_subset_y) until there is none
                while any(~isnan(Distance_obs_inArea)) 
                    % best candidate: nearest raw obs to the ith grid
                    [obs_distances(idx_subset_x,idx_subset_y), nearest_obs_index] = min(Distance_obs_inArea);
                    my_obs_index_candidate = idx_obs_good(nearest_obs_index); % ith raw obs
                    % use this obs candidate if it hasn't been selected for use yet
                    if (~any(obs_index==my_obs_index_candidate))
                        my_obs_index = my_obs_index_candidate;
                        break
                    % otherwise keep searching for other raw-obs candidates
                    else
                        Distance_obs_inArea(nearest_obs_index) = NaN;
                        obs_was_pushed_from_top_candidate(idx_subset_x,idx_subset_y) = obs_was_pushed_from_top_candidate(idx_subset_x,idx_subset_y) + 1;
%                       disp(strcat('ob did not take top candidate at x = ',num2str(idx_subset_x),' , y = ',num2str(idx_subset_y)));
                    end
                end
                % Other conditions
                if control.min_Tb_threshold && (Tb_obs_col(my_obs_index) < 0 )
                    obs_index(idx_subset_x,idx_subset_y) = NaN;
                elseif control.max_Tb_threshold && (Tb_obs_col(my_obs_index) > 999)
                    obs_index(idx_subset_x,idx_subset_y) = NaN;
                else
                    obs_index(idx_subset_x,idx_subset_y) = my_obs_index; %% !
                end
                % Record the number of raw-obs candidates for ith grid
                obs_n_candidates(idx_subset_x,idx_subset_y) = length(Distance_obs_inArea);
            end
        % ------------------------------------------------------------

        % Flag if no raw-obs candidates are used 
        %disp(strcat('    obs did not take ANY candidate at x = ',num2str(x),' , y = ',num2str(y)));

        end % end loop: for y = index_y
    end % end loop: for x = index_x




end

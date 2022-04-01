function [obs_index] = PickRawforCRTM(lat_raw,lon_raw,Tb_raw,min_XLONG,max_XLONG,min_XLAT,max_XLAT,latitudes_crtm,longitudes_crtm,slots_x,slots_y,control)

    lat_col = reshape(lat_raw,[],1);
    lon_col = reshape(lon_raw,[],1);
    Tb_col = reshape(Tb_raw,[],1);

    obs_distances = nan(length(slots_x),length(slots_y));
    obs_index = nan(length(slots_x),length(slots_y));
    obs_n_candidates = nan(length(slots_x),length(slots_y));
    obs_was_pushed_from_top_candidate = zeros(size(obs_n_candidates));

    n_points_noObs = 0;

    for x = slots_x % operation on grids on x direction
        %x
        idx_obs_slots_x = find(slots_x==x);
        % reduce by-distance search of possible closest location with longitude bounds
        myLon = longitudes_crtm(x);
        lg_x = (lon_col > max(min_XLONG, myLon - control.search_buffer)) & (lon_col < min(max_XLONG, myLon + control.search_buffer));
        
        for y = slots_y
            clear myTb_inArea_distance
            idx_obs_slots_y = find(slots_y==y);
            % reduce by-distance search of possible closest location with latitude bounds
            myLat = latitudes_crtm(y);
            lg_y = (lat_col > max(min_XLAT, myLat - control.search_buffer*1.5)) & (lat_col <  min(max_XLAT, myLat + control.search_buffer*1.5));

            lg_good = ((lg_x + lg_y) == 2);
            if max(lg_good) == 0
                continue;
            end

            idx_search_area = find(lg_good);
            lat_good = lat_col(idx_search_area);
            lon_good = lon_col(idx_search_area);
            Tb_good = Tb_col(idx_search_area);
            
            % For the grid point of interest in the loop, calculate the distance between it and grid points in the search area   
            for i = 1:length(idx_search_area)
                myTb_inArea_distance(i) = distance(lat_good(i),lon_good(i),myLat,myLon);
            end
            
            % ------------- Select the best raw-obs candidate -------------
            %if (~isempty(idx_search_area))
            if sum(idx_search_area) > 0
                my_obs_index = NaN;
                % --- Keep searching for the raw-obs candidates for ith grid
                % ---(idx_obs_slots_x,idx_obs_slots_y) until there is none
                while any(~isnan(myTb_inArea_distance)) 
                    % best candidate: nearest raw obs to the ith grid
                    [obs_distances(idx_obs_slots_x,idx_obs_slots_y), nearest_obs_index] = min(myTb_inArea_distance);
                    my_obs_index_candidate = idx_search_area(nearest_obs_index); % ith raw obs
                    % use this obs candidate if it hasn't been selected yet
                    %if (~any(any(obs_index==my_obs_index_candidate)))
                    if (~any(obs_index==my_obs_index_candidate))
                        my_obs_index = my_obs_index_candidate;
                        break
                    % otherwise keep searching for other raw-obs candidates
                    else
                        myTb_inArea_distance(nearest_obs_index) = NaN;
                        obs_was_pushed_from_top_candidate(idx_obs_slots_x,idx_obs_slots_y) = obs_was_pushed_from_top_candidate(idx_obs_slots_x,idx_obs_slots_y) + 1;
%                                 disp(strcat('ob did not take top candidate at x = ',num2str(idx_obs_slots_x),' , y = ',num2str(idx_obs_slots_y)));
                    end
                end
                % Flag if no raw-obs candidates are used 
                if (isnan(my_obs_index))
                  if ( mod(n_points_noObs, 100) == 0)
                    disp(strcat('ob did not take ANY candidate at x = ',num2str(idx_obs_slots_x),' , y = ',num2str(idx_obs_slots_y)));
                    n_points_noObs = n_points_noObs + 1;
                  end 
                end
                % Other conditions
                if (exist('min_threshold') && Tb_good(my_obs_index) < min_threshold)
                    obs_index(idx_obs_slots_x,idx_obs_slots_y) = NaN;
                elseif (exist('max_threshold') && Tb_good(my_obs_index) > max_threshold)
                    obs_index(idx_obs_slots_x,idx_obs_slots_y) = NaN;
                else
                    obs_index(idx_obs_slots_x,idx_obs_slots_y) = my_obs_index; %% !
                end
                % Record the number of raw-obs candidates for ith grid
                obs_n_candidates(idx_obs_slots_x,idx_obs_slots_y) = length(myTb_inArea_distance);
            end
            % ------------------------------------------------------------
        end % end loop: for y = slots_y
    end % end loop: for x = slots_x

end

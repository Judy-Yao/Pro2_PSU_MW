function [if_swath_good, ChName_perSwath] = Combine_AMSR2(storm_phase, platform, sensor, min_lat_dpy, max_lat_dpy, min_lon_dpy, max_lon_dpy, time, loc_storm, Swath_used, if_swath_good, DAtime_perSwath, ChIdx_perSwath, ChName_perSwath, Tb_file, control)

    num_89GHz_scan = 0;
    sw_89GHz_scan =[];
    % identify if both 89 GHz scans are used
    for i_sw = 1:length(Swath_used)
        if contains(ChName_perSwath{i_sw},'89GHzV-PolA-Scan') | contains(ChName_perSwath{i_sw},'89GHzV-PolB-Scan')
            num_89GHz_scan = num_89GHz_scan + 1;
            sw_89GHz_scan = [sw_89GHz_scan, i_sw];
        end
    end
    
    if num_89GHz_scan == 2
        if DAtime_perSwath{sw_89GHz_scan(1)} == DAtime_perSwath{sw_89GHz_scan(2)}
            % if DA_time is the same for both scans, combine them
            disp(['Combine A scan and B scan of 89GHzV-pol together...']);
            for i_sw = 1:length(Swath_used)
                if contains(ChName_perSwath{i_sw},'89GHzV-PolA-Scan')
                    xlat_Ascan = h5read(Tb_file,[Swath_used{i_sw},'/Latitude']);
                    xlon_Ascan = h5read(Tb_file,[Swath_used{i_sw},'/Longitude']);
                    idx_Tb_oneCh = ChIdx_perSwath(i_sw);
                    Tb_Ascan_Chs = h5read(Tb_file,[Swath_used{i_sw},'/Tc']);
                    Tb_Ascan = squeeze(Tb_Ascan_Chs(idx_Tb_oneCh,:,:));
                elseif contains(ChName_perSwath{i_sw},'89GHzV-PolB-Scan')
                    xlat_Bscan = h5read(Tb_file,[Swath_used{i_sw},'/Latitude']);
                    xlon_Bscan = h5read(Tb_file,[Swath_used{i_sw},'/Longitude']);
                    idx_Tb_oneCh = ChIdx_perSwath(i_sw);
                    Tb_Bscan_Chs = h5read(Tb_file,[Swath_used{i_sw},'/Tc']); 
                    Tb_Bscan = squeeze(Tb_Bscan_Chs(idx_Tb_oneCh,:,:));
                else
                    continue;
                end
            end
            xlat = [xlat_Ascan,xlat_Bscan];
            xlon = [xlon_Ascan, xlon_Bscan];
            Tb_oneCh = [Tb_Ascan, Tb_Bscan];
            % masked Tbs, lat, lon based on the area of interest
            idx_mask = (xlon < max_lon_dpy) & (xlon > min_lon_dpy) & (xlat < max_lat_dpy) & (xlat > min_lat_dpy);
            sum(idx_mask, 'all')
            lat_need = xlat(idx_mask);
            lon_need = xlon(idx_mask);
            Tb_need = Tb_oneCh(idx_mask);
            % plot Tbs 
            %disp(['best-track location:', num2str(loc_storm{i_sw})]);
            Plot_AMSR2scans(storm_phase, platform, sensor, min_lat_dpy, max_lat_dpy, min_lon_dpy, max_lon_dpy, time, loc_storm{i_sw}, lat_need, lon_need, Tb_need, DAtime_perSwath{i_sw},control);
            % remove these two scans
            if_swath_good(sw_89GHz_scan) = false;
        else
            disp('DA time of scans are different. Skip the combintaion.');
            return; % exit a function
        end        
    else
        disp(['Only a scan of 89GHz exist. Skip the combination.']);
        return;
    end



end

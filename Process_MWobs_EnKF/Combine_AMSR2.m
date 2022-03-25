function [Swath_used,ChIdx_ps,ChName_ps,lat,lon,Tb,zenith,sat_lat,sat_lon,sat_alt,sat_azimuth,outime] = Combine_AMSR2(iTb,Swath_used,ChIdx_ps,ChName_ps,DAtime_ps,lat,lon,Tb,zenith,sat_lat,sat_lon,sat_alt,sat_azimuth,outime)

% The combination of 89GHz scans in AMSR2 sensor is based on the assumption that all parameters (especially zenith angle) of B scan are the same as of A scan
% In other words, under such an assumption, 89GHz B scan just provides more measurements of the same kind data to 89GHz A scan

    % double check 
    num_89GHz_scan = 0;
    idx_89GHz =[];
    % identify if both 89 GHz scans are used
    for i_sw = 1:length(Swath_used{iTb})
        if contains(ChName_ps{iTb}(i_sw),'89GHzV-PolA-Scan') | contains(ChName_ps{iTb}(i_sw),'89GHzV-PolB-Scan')
            num_89GHz_scan = num_89GHz_scan + 1;
            idx_89GHz = [idx_89GHz, i_sw];
        end
    end  

    if num_89GHz_scan ~= 2
        disp(['Only a scan of 89GHz exist. Skip the combination.']);
        return 
    else
        if DAtime_ps{iTb}(idx_89GHz(1)) == DAtime_ps{iTb}(idx_89GHz(2))
            % if DA_time is the same for both scans, combine them
            disp(['Combine A scan and B scan of 89GHzV-pol together...']);
            Clat = [lat{idx_89GHz(1)},lat{idx_89GHz(2)}];
            lat{idx_89GHz(1)} = Clat; lat{idx_89GHz(2)} = [];
            Clon = [lon{idx_89GHz(1)},lon{idx_89GHz(2)}];
            lon{idx_89GHz(1)} = Clon; lon{idx_89GHz(2)} = [];
            CTb = [Tb{idx_89GHz(1)},Tb{idx_89GHz(2)}];
            Tb{idx_89GHz(1)} = CTb; Tb{idx_89GHz(2)} = [];
            Czenith = [zenith{idx_89GHz(1)},zenith{idx_89GHz(2)}];
            zenith{idx_89GHz(1)} = Czenith; zenith{idx_89GHz(2)} = [];
            Csat_lat = [sat_lat{idx_89GHz(1)},sat_lat{idx_89GHz(2)}];
            sat_lat{idx_89GHz(1)} = Csat_lat; sat_lat{idx_89GHz(2)} = [];
            Csat_lon = [sat_lon{idx_89GHz(1)},sat_lon{idx_89GHz(2)}];
            sat_lon{idx_89GHz(1)} = Csat_lon; sat_lon{idx_89GHz(2)} = [];
            Csat_alt = [sat_alt{idx_89GHz(1)},sat_alt{idx_89GHz(2)}];
            sat_alt{idx_89GHz(1)} = Csat_alt; sat_alt{idx_89GHz(2)} = [];
            Csat_azimuth = [sat_azimuth{idx_89GHz(1)},sat_azimuth{idx_89GHz(2)}];
            sat_azimuth{idx_89GHz(1)} = Csat_azimuth; sat_azimuth{idx_89GHz(2)} = [];
            Coutime = [outime{idx_89GHz(1)},outime{idx_89GHz(2)}];
            outime{idx_89GHz(1)} = Coutime; outime{idx_89GHz(2)} = [];
            % find indices of 89 A scan and B scan among all possible channels in the AMSR2
            for ith = 1:length(idx_89GHz)
                if contains(ChName_ps{iTb}{idx_89GHz(ith)},'89GHzV-PolB-Scan')
                    idx_Bscan = idx_89GHz(ith);
                end
            end 
            %idx_Ascan = idx_89GHz(idx_89GHz ~= idx_Bscan);
            % Get rid of information about B scan
            swath_AB = Swath_used{iTb};
            swath_AB(idx_Bscan) = [];
            Swath_used{iTb} = swath_AB;

            ChIdx_AB = ChIdx_ps{iTb};
            ChIdx_AB(idx_Bscan) = [];
            ChIdx_ps{iTb} = ChIdx_AB;
            
            ChName_AB = ChName_ps{iTb};
            ChName_AB(idx_Bscan) = [];
            ChName_ps{iTb} = ChName_AB;
    end

    size(Clat)

end

% This function combines the two 89GHz scans in AMSR2 sensor
% Note: The combination of 89GHz scans in AMSR2 sensor is based on the assumption that all parameters (especially zenith angle) of B scan are the same as of A scan. &
% In other words, under such an assumption, 89GHz B scan just provides more measurements of the same kind data to 89GHz A scan.

function [Swath_used,ChIdx_all,ChName_all,lat,lon,Tb,zenith,sat_lat,sat_lon,sat_alt,sat_azimuth,outime] = Handle_AMSR2(iTb,Swath_used,ChIdx_all,ChName_all,DAtime_all,lat,lon,Tb,zenith,sat_lat,sat_lon,sat_alt,sat_azimuth,outime,control)

    % Identify how many 89 GHz scans are used  
    num_89GHz_scan = 0;
    idx_89GHz =[];
    for i_sw = 1:length(Swath_used{iTb})
        if contains(ChName_all{iTb}(i_sw),'89GHzV-PolA-Scan') | contains(ChName_all{iTb}(i_sw),'89GHzV-PolB-Scan')
            num_89GHz_scan = num_89GHz_scan + 1;
            idx_89GHz = [idx_89GHz, i_sw];
        end
    end 
	
    if control.NOTuse_AMSR89GHz & (~control.comnine_AMSR89GHz) % 89-GHz data of AMSR2 is not in use.
        lat(idx_89GHz) = [];
		lon(idx_89GHz) = [];
        Tb(idx_89GHz) = [];
        zenith(idx_89GHz) = [];
        sat_lat(idx_89GHz) = [];
        sat_lon(idx_89GHz) = [];
        sat_alt(idx_89GHz) = [];
        sat_azimuth(idx_89GHz) = [];
        outime(idx_89GHz) = [];
			
        Swath_used{iTb}(idx_89GHz) = [];
        ChIdx_all{iTb}(idx_89GHz) = [];
        ChName_all{iTb}(idx_89GHz) = [];
        disp('        89 GHz of AMSR2 is not in use!');
    elseif (~control.NOTuse_AMSR89GHz) & control.comnine_AMSR89GHz % 89-GHz data of AMSR2 is in use.
        if num_89GHz_scan ~= 2
            disp('        Only a scan of 89GHz exists. Skip the combination.');
            return 
        else % If both 89 GHz scans are used
            if DAtime_all{iTb}(idx_89GHz(1)) == DAtime_all{iTb}(idx_89GHz(2)) %if DA_time is the same for both scans, combine them
                disp(['        (sanity check) size of one 89 GHZ scan: ',num2str(size(lat{idx_89GHz(1)}))]);
                disp(['        (sanity check) size of the other 89 GHZ scan: ',num2str(size(lat{idx_89GHz(2)}))]);
                disp(['        Combine A scan and B scan of 89GHzV-pol together...']);
                Clat = [lat{idx_89GHz(1)},lat{idx_89GHz(2)}];
                lat{idx_89GHz(1)} = Clat; lat(idx_89GHz(2)) = [];
                Clon = [lon{idx_89GHz(1)},lon{idx_89GHz(2)}];
                lon{idx_89GHz(1)} = Clon; lon(idx_89GHz(2)) = [];
                CTb = [Tb{idx_89GHz(1)},Tb{idx_89GHz(2)}];
                Tb{idx_89GHz(1)} = CTb; Tb(idx_89GHz(2)) = [];
                Czenith = [zenith{idx_89GHz(1)},zenith{idx_89GHz(2)}];
                zenith{idx_89GHz(1)} = Czenith; zenith(idx_89GHz(2)) = [];
                Csat_lat = [sat_lat{idx_89GHz(1)},sat_lat{idx_89GHz(2)}];
                sat_lat{idx_89GHz(1)} = Csat_lat; sat_lat(idx_89GHz(2)) = [];
                Csat_lon = [sat_lon{idx_89GHz(1)},sat_lon{idx_89GHz(2)}];
                sat_lon{idx_89GHz(1)} = Csat_lon; sat_lon(idx_89GHz(2)) = [];
                Csat_alt = [sat_alt{idx_89GHz(1)},sat_alt{idx_89GHz(2)}];
                sat_alt{idx_89GHz(1)} = Csat_alt; sat_alt(idx_89GHz(2)) = [];
                Csat_azimuth = [sat_azimuth{idx_89GHz(1)},sat_azimuth{idx_89GHz(2)}];
                sat_azimuth{idx_89GHz(1)} = Csat_azimuth; sat_azimuth(idx_89GHz(2)) = [];
                Coutime = [outime{idx_89GHz(1)},outime{idx_89GHz(2)}];
                outime{idx_89GHz(1)} = Coutime; outime(idx_89GHz(2)) = [];
                disp(['        (sanity check) size of the final 89 GHZ scan: ',num2str(size(lat{idx_89GHz(1)}))]);
                % Find indices of 89 A scan and B scan among all possible channels in the AMSR2
                for ith = 1:length(idx_89GHz)
                    if contains(ChName_all{iTb}{idx_89GHz(ith)},'89GHzV-PolB-Scan')
                        idx_Bscan = idx_89GHz(ith);
                    end
                end 
                %idx_Ascan = idx_89GHz(idx_89GHz ~= idx_Bscan);
                % get rid of B scan
                swath_AB = Swath_used{iTb};
                swath_AB(idx_Bscan) = [];
                Swath_used{iTb} = swath_AB;

                ChIdx_AB = ChIdx_all{iTb};
                ChIdx_AB(idx_Bscan) = [];
                ChIdx_all{iTb} = ChIdx_AB;
            
                ChName_AB = ChName_all{iTb};
                ChName_AB(idx_Bscan) = [];
                ChName_all{iTb} = ChName_AB;
        end
    end

end

function [] = Singlepass_write(idx_usedTb,istorm,Swath_used,loc_DAtime_ps,Tb_file,control)

	[filepath,filename,filext] = fileparts(Tb_file); sensor_info = split(filename,'.'); platform = sensor_info(2);

	[sat_name,myLat,myLon,myTb,mySat_lat,mySat_lon,mySat_alt,mySat_azimuth,myScan_angle,myZenith_angle,myFov_crossTrack,myFov_alongTrack,myTimes,myChNum,myRoi_hydro,myRoi_otherVars,myObsErr] = ProduceforEnKF(idx_usedTb,Swath_used,filename,loc_DAtime_ps,control); 

   % modify satellite-and-sensor name so that it is consistent with what is used in the CRTM package 
	if (matches(platform,'GCOMW1'))
        sat_name = 'amsr2_gcom-w1';
    elseif (matches(platform,'NPP'))
        sat_name = 'atms_npp';
    elseif (matches(platform,'GPM'))
        sat_name = 'gmi_gpm';
    elseif (matches(platform,'METOPA'))
        sat_name = 'mhs_metop-a';
    elseif (matches(platform,'METOPB'))
        sat_name = 'mhs_metop-b';
    elseif (matches(platform,'NOAA18'))
        sat_name = 'mhs_n18';
    elseif (matches(platform,'NOAA19'))
        sat_name = 'mhs_n19';
    elseif (matches(platform,'MT1'))
        sat_name = 'saphir_meghat';
    elseif (matches(platform,'F16'))
        sat_name = 'ssmis_f16';
    elseif (matches(platform,'F17'))
        sat_name = 'ssmis_f17';
    elseif (matches(platform,'F18'))
        sat_name = 'ssmis_f18';
    end
    % create text file name
    filename = strcat(control.output_dir,control.storm_phase(istorm),'/microwave_d03_',DAtime_ps{idx_usedTb}(1),'_so');
    formatSpec = '%12s%16s%12i%12.3f%12.3f%12.3f%12i%12i%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n';
    fileID = fopen(filename,'w');
    % reshape values into columns
    out_lat = cat(1, myLat{:}); out_lon = cat(1, myLon{:}); out_Tb = cat(1, myTb{:});
    out_sat_lat = cat(1, mySat_lat{:}); out_sat_lon = cat(1, mySat_lon{:}); out_sat_alt = cat(1, mySat_alt{:});
    out_Sat_azimuth = cat(1, mySat_azimuth{:}); out_Scan_angle = cat(1, myScan_angle{:});
    out_Zenith_angle = cat(1, myZenith_angle{:}); out_Fov_crossTrack = cat(1, myFov_crossTrack{:}); out_Fov_alongTrack = cat(1, myFov_alongTrack{:});
    out_times = cat(1, myTimes{:}); out_chNum = cat(1, myChNum{:}); out_obsErr = cat(1, myObsErr{:});
    out_ROI_other = cat(1, myRoi_otherVars{:}); out_ROI_hydro = cat(1, myRoi_hydro{:});
    % write values to a file
    for rd = 1:length(out_Tb)
        if (strcmp(sat_name,'gmi_gpm'))
            if (out_chNum(rd) < 9)
                out_sat_name = 'gmi_gpm_lf';
            else
                out_sat_name = 'gmi_gpm_hf';
            end
        else
            out_sat_name = sat_name;
        end
        fprintf(fileID, formatSpec, ...
                out_times(rd), out_sat_name, out_chNum(rd),...
                out_lat(rd), out_lon(rd), out_Tb(rd), ...
                out_ROI_hydro(rd), out_ROI_other(rd), out_obsErr(rd), ...
                out_Fov_crossTrack(rd), out_Fov_alongTrack(rd),...
                out_Scan_angle(rd), out_Zenith_angle(rd), out_Sat_azimuth(rd),...
                out_sat_lat(rd), out_sat_lon(rd), out_sat_alt(rd));
    end
    
    fclose(fileID);





end

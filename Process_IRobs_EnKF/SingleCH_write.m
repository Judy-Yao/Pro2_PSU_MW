function [] = SingleCH_write(istorm, DAtime, Tb_file, control)

    disp(['  Processing obs for file: '+ Tb_file]);
    
    [~, filename, ~] = fileparts(Tb_file);
    file_info = split(filename,'_'); Ch = str2num(erase(file_info(6),'CH'));
    [myTimes,mySat_name,myChNum,myLat,myLon,myTb,myROI_hydro,myROI_other,myObsErr,mySat_alt] = ProduceForEnKF(istorm, DAtime, Ch, Tb_file, control);
     
    % create text file name
    filename = strcat(control.output_dir,control.storm_phase{istorm},'/radiance_d03_',DAtime,'_so');
    disp(['  Output processed GOESR obs file: ', filename]);
    formatSpec = '%12s%12s%12i%12.3f%12.3f%12.3f%12i%12i%12.3f%12.3f\n';
    fileID = fopen(filename,'w');

    % reshape values into columns
    out_times = cat(1, myTimes{:});
    out_sat_name = cat(1, mySat_name{:}); out_chNum = cat(1, myChNum{:});
    out_lat = cat(1, myLat{:}); out_lon = cat(1, myLon{:}); out_Tb = cat(1, myTb{:});
    out_ROI_hydro = cat(1, myROI_hydro{:}); out_ROI_other = cat(1, myROI_other{:});
    out_obsErr = cat(1, myObsErr{:}); out_Sat_alt = cat(1, mySat_alt{:});

    % write values to a file
    for rd = 1:length(out_Tb)
        fprintf(fileID, formatSpec, ...
                out_times(rd), out_sat_name(rd), out_chNum(rd),...
                out_lat(rd), out_lon(rd), out_Tb(rd), ...
                out_ROI_hydro(rd), out_ROI_other(rd), out_obsErr(rd), out_Sat_alt(rd));
    end
    fclose(fileID);

end

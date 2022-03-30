% Find time and center location of the storm in Best track file 

function [bestrack_str] = Bestrack_read(istorm, control)

    bestrack_dir = [control.bestrack_dir,control.storm_phase{istorm},'/*'];
    bestrack_dir = ls(bestrack_dir);
    bestrack_file = regexprep(bestrack_dir,'\n$','');
    fid = fopen(bestrack_file);
    best_track = textscan(fid,'%s','delimiter','');
    fclose(fid);
    % For every record, separate substrings(values)
    best_track_str_vec = string(best_track{1,1}(:));
    len_record = length(best_track{1,1}(:)); % how many records there are
    % Keep necessary information
    bestrack_str = cell(len_record,3); %time,lat,lon
    for ir=1:len_record
        record_substr = strsplit(best_track_str_vec(ir),',');
        %Time
        bestrack_str{ir,1} = append(erase(record_substr(1,3),' '),'00'); % %YYYYMMDDHHmm
        %Latitude 
        bestrack_str{ir,2} = record_substr(1,7); %Latitude for the DTG: 0 - 900 tenths of degrees (e.g.,161N)
        bestrack_str{ir,2} = str2double(strrep(bestrack_str{ir,2},'N',''))/10; % 0 - 90 degree
        %Longitude
        bestrack_str{ir,3} = record_substr(1,8); %Longitude for the DTG: 0 - 1800 tenths of degrees (e.g.,269W)
        if contains(bestrack_str{ir,3},'W')
            bestrack_str{ir,3} =  0-str2double(strrep(bestrack_str{ir,3},'W',''))/10; 
        else
            bestrack_str{ir,3} = str2double(strrep(bestrack_str{ir,3},'E',''))/10; 
        end
    end

end

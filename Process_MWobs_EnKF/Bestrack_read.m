% =============================================================================================================================
% Read times and locations of the storm from the best-track file
% =============================================================================================================================

function [bestrack_str] = Bestrack_read(istorm, control)
    
    % List best-track file(s) in the directory and separate the names
    bestrack_dir = [control.bestrack_dir,control.storm_phase{istorm},'/*'];
    bestrack_dir = ls(bestrack_dir);
    bestrack_file = regexprep(bestrack_dir,'\n$','');
    % Read the formatted data as character into a cell array 
    fid = fopen(bestrack_file);
    best_track = textscan(fid,'%s','delimiter','');
    fclose(fid);
    % Count the number of records
    len_record = length(best_track{1}(:)); 
    % Convert data type from character to string 
    best_track_str_vec = string(best_track{1}(:));

    % Keep necessary information
    bestrack_str = cell(len_record,3); %time,lat,lon
    for ir=1:len_record
        % Split one record into different parts/variables
        record_substr = strsplit(best_track_str_vec(ir),',');
        % Read time
        bestrack_str{ir,1} = append(erase(record_substr(1,3),' '),'00'); % %YYYYMMDDHHmm
        % Read latitude 
        bestrack_str{ir,2} = record_substr(1,7); %Latitude for the DTG: 0 - 900 tenths of degrees (e.g.,161N)
        if contains(bestrack_str{ir,2},'N')
            bestrack_str{ir,2} = str2double(strrep(bestrack_str{ir,2},'N',''))/10; % 0 - 90 degree
        else
            bestrack_str{ir,2} = 0-str2double(strrep(bestrack_str{ir,2},'S',''))/10; 
        end
        %Longitude
        bestrack_str{ir,3} = record_substr(1,8); %Longitude for the DTG: 0 - 1800 tenths of degrees (e.g.,269W)
        if contains(bestrack_str{ir,3},'W')
            bestrack_str{ir,3} =  0-str2double(strrep(bestrack_str{ir,3},'W',''))/10; 
        else
            bestrack_str{ir,3} = str2double(strrep(bestrack_str{ir,3},'E',''))/10; 
        end
    end

end

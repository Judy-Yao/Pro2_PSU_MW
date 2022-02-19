function [Swath_used, ChIdx_perSwath, ChName_perSwath] = Swath_Channel(Tb_file, control)
% For the frequencies of interest, obtain the swaths which include the channels of interest and the channel numbers under these swathes.
    file_info = h5info(Tb_file);
    num_swath = length(file_info.Groups); % total number of swaths per file

    Swath_used= {}; % names of swaths which include the channels of interest
    ChIdx_perSwath = []; % index of channels of interest under these swaths
    ChName_perSwath = {}; % names of channels of interest under these swaths


    % First go over favoriate Channels
    for i_sw = 1:num_swath
        AllCh_name_bad = file_info.Groups(i_sw).Datasets(4).Attributes(6).Value; % several lines with white spaces and extra charaters
        AllCh_name_good = char(join(erase(splitlines(erase(strtrim(AllCh_name_bad),'and')), ' '))); % one single line

        for i_ch = 1:length(control.favCh)
            if contains(AllCh_name_good, control.favCh{i_ch})
               oneCh_name_idx_start = strfind(AllCh_name_good,control.favCh{i_ch});
               idx_one_Ch = cellfun(@str2num, extractBetween(AllCh_name_good,oneCh_name_idx_start-2,oneCh_name_idx_start-2));
               % Save
               Swath_used{end+1} = file_info.Groups(i_sw).Name;
               ChName_perSwath{end+1} = control.favCh{i_ch};
               ChIdx_perSwath = [ChIdx_perSwath,idx_one_Ch];
            else
               continue;
            end       
        end
    end

    % Then go over similar Channels
    if length(ChName_perSwath) < 2

        for i_sw = 1:num_swath
            AllCh_name_bad = file_info.Groups(i_sw).Datasets(4).Attributes(6).Value;
            AllCh_name_good = char(join(erase(splitlines(erase(strtrim(AllCh_name_bad),'and')), ' ')));

            for i_ch = 1:length(control.favCh_sup)
                if contains(AllCh_name_good, control.favCh_sup{i_ch})
                   oneCh_name_idx_start = strfind(AllCh_name_good, control.favCh_sup{i_ch});
                   idx_one_Ch = cellfun(@str2num, extractBetween(AllCh_name_good,oneCh_name_idx_start-2,oneCh_name_idx_start-2));
                   % Save
                   Swath_used{end+1} = file_info.Groups(i_sw).Name;
                   ChName_perSwath{end+1} =  control.favCh_sup{i_ch};
                   ChIdx_perSwath = [ChIdx_perSwath,idx_one_Ch];
                else
                   continue;
                end       
            end
        end

    end
    
    disp('Frequencies of interest are:');
    for item = 1:length(ChName_perSwath)
        disp(['        ',ChName_perSwath{item}]);
    end

end

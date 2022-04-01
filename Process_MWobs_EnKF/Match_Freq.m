% This function reads frequencies in a L1C MW obs file & 
% loops the frequenies of interest listed in control.favF &
% finds the frequency(ies) in the L1C file that matches the frequency(ies) of interest to the study &
% returns the location of the frequency(ies) of interest (which swath and which index under that swath)

% Note: in the current study, we are only interested in a low frequency around 18 GHz and a high frequency around 183 GHz (if 183 GHz doesn' t exist, we will use 89 GHz instead)
% Note: this algorithm only works if every frequency name is unique. 
% Note: each frequency is treated as an object/item as well as its swath and index under that swath 

function [Swath_used, ChIdx_perSwath, ChName_perSwath] = Match_Freq(Tb_file, control)
    
    file_info = h5info(Tb_file);
    num_swath = length(file_info.Groups); % total number of swaths per L1C file

    Swath_used= []; % (strings) name(s) of swaths which include the frequency(ies) of interest
    ChName_perSwath = []; % (strings) name(s) of frequency(ies) of interest under the swath(s)
    ChIdx_perSwath = []; % (double) index(ices) of frequency(ies) of interest under the swath(s)

    % loop through each swath in L1C file
    for i_sw = 1:num_swath
        % Note: name(s) of frequency(ies) is(are) buried in the attributes of the Tc variable
        AllCh_name_bad = file_info.Groups(i_sw).Datasets(4).Attributes(6).Value; % several lines with white spaces and extra charaters
        AllCh_name_good = char(join(erase(splitlines(erase(strtrim(AllCh_name_bad),'and')), ' '))); % one single line
        % loop through each frequency in control.favCh
        for i_ch = 1:length(control.favFreq)
            if contains(AllCh_name_good, control.favFreq{i_ch})
               % Locate the characters [name(s) of frequency(ies)]
               oneCh_name_idx_start = strfind(AllCh_name_good,control.favFreq{i_ch});
               idx_one_Ch = cellfun(@str2num, extractBetween(AllCh_name_good,oneCh_name_idx_start-2,oneCh_name_idx_start-2));   
               % Save
               Swath_used = [Swath_used,string(file_info.Groups(i_sw).Name)]; %file_info.Groups(i_sw).Name is Character. It is converted to strings in this case. E.g., '/S1' is converted to "/S1"
               ChName_perSwath = [ChName_perSwath, string(control.favCh{i_ch})];
               ChIdx_perSwath = [ChIdx_perSwath,idx_one_Ch];
            else
               continue;
            end       
        end
    end


    % Sanity Check
    if length(Swath_used) ~= length(ChName_perSwath) || length(Swath_used) ~= length(ChIdx_perSwath)
        disp('Error getting used swathes and channels!');
    end
    
end

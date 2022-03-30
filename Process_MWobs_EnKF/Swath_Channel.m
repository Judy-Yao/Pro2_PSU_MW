% This function reads frequencies in a L1C MW obs file & 
% loops the frequenies of interest listed in control.favCh/control.favCh_sup &
% finds the frequency(ies) in the L1C file that matches the frequency(ies) of interest to the study &
% returns the location of the frequency(ies) of interest (which swath and which index under that swath)
% Note: in the current study, we are only interested in a low frequency around 18 GHz and a high frequency 183 GHz (if 183 GHz doesn' t exist, we will use 89 GHz instead)

function [Swath_used, ChIdx_perSwath, ChName_perSwath] = Swath_Channel(Tb_file, control)
    
    file_info = h5info(Tb_file);
    num_swath = length(file_info.Groups); % total number of swaths per file

    Swath_used= []; % (strings) name(s) of swaths which include the frequency(ies) of interest
    ChName_perSwath = []; % (strings) name(s) of frequency(ies) of interest under the swath(s)
    ChIdx_perSwath = []; % (double) index(ices) of frequency(ies) of interest under the swath(s)

    % -------------- First go over favoriate channels/frequencies --------------
    % loop through each swath in L1C file
    for i_sw = 1:num_swath
        % Note: name(s) of frequency(ies) is(are) hidden behind the attributes of the Tc variable
        AllCh_name_bad = file_info.Groups(i_sw).Datasets(4).Attributes(6).Value; % several lines with white spaces and extra charaters
        AllCh_name_good = char(join(erase(splitlines(erase(strtrim(AllCh_name_bad),'and')), ' '))); % one single line
        % loop through each frequency in control.favCh
        for i_ch = 1:length(control.favCh)
            if contains(AllCh_name_good, control.favCh{i_ch})
               % Locate the characters [name(s) of frequency(ies)]
               oneCh_name_idx_start = strfind(AllCh_name_good,control.favCh{i_ch});
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

    % -------------- Then go over similar channels/frequencies -------------- 
    if length(ChName_perSwath) == 2 % if the condition is met, it means that the two frequencies of interest are already found 
        % Sanity Check
        if length(Swath_used) ~= length(ChName_perSwath) || length(Swath_used) ~= length(ChIdx_perSwath)
            disp('Error getting used swathes and channels!');
        end
 
        return % Frequencies of interest are found; exit the function  
    else
        % loop through each swath in L1C file
        for i_sw = 1:num_swath
            AllCh_name_bad = file_info.Groups(i_sw).Datasets(4).Attributes(6).Value;
            AllCh_name_good = char(join(erase(splitlines(erase(strtrim(AllCh_name_bad),'and')), ' ')));
            % loop through each frequency in control.favCh_sup
            for i_ch = 1:length(control.favCh_sup)
                if contains(AllCh_name_good, control.favCh_sup{i_ch})
                   oneCh_name_idx_start = strfind(AllCh_name_good, control.favCh_sup{i_ch});
                   idx_one_Ch = cellfun(@str2num, extractBetween(AllCh_name_good,oneCh_name_idx_start-2,oneCh_name_idx_start-2));
                   % Save
				   Swath_used = [Swath_used,string(file_info.Groups(i_sw).Name)];
				   ChName_perSwath = [ChName_perSwath, string(control.favCh_sup{i_ch})];           
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
    
    %disp('Frequencies of interest are:');
    %for item = 1:length(ChName_perSwath)
    %    disp(["        " + ChName_perSwath(item)]);
    %end
	

end

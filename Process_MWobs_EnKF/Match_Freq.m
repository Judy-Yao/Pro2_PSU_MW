% ==========================================================================================
% Generally, this function identifies the relevant attributes to the frequencies of interest
% ==========================================================================================
% HDF5:
% This function reads in available frequencies from a L1C MW obs file in HDF5 format &
% loops the frequenies of interest listed in control.favF &
% identifies the frequency(ies) in the L1C file that matches the frequency(ies) of interest to the study & 
% returns the location of the frequency(ies) of interest (under which swath and with which index under that swath)
%
% NC:
% Instead of reading in the swath name for a specific frequency, this function manually assigns "lores" or "hires" for the frequenies of interest listed in control.favF 
%
% Note:
% - In the current study, we are only interested in a low frequency around 18 GHz and a high frequency around 183 GHz (if 183 GHz doesn' t exist, we will use 89 GHz instead)
% - Each frequency is treated as an object/item as well as its swath and index under that swath 
% ------------------------------------------------------------------------------------------------------------------

function [Swath_used, ChIdx_perSwath, ChName_perSwath] = Match_Freq(isensor, Tb_file, control)
   
	% Determine if the file should be read by ncread or h5read
    [~,~,filext] = fileparts(Tb_file); 

    Swath_used= []; % (strings) name(s) of swaths which include the frequency(ies) of interest
    ChName_perSwath = []; % (strings) name(s) of frequency(ies) of interest under the swath(s)
    ChIdx_perSwath = []; % (double) index(ices) of frequency(ies) of interest under the swath(s)
   
	% HDF5 file
	if contains(filext,"HDF5")
		file_info = h5info(Tb_file);
		num_swath = length(file_info.Groups); % total number of swaths per L1C file

		% loop through each swath in L1C file
		for i_sw = 1:num_swath
			% Note: name(s) of frequency(ies) is(are) buried in the attributes of the Tc variable
			AllCh_name_bad = file_info.Groups(i_sw).Datasets(4).Attributes(6).Value; % several lines with white spaces and extra charaters
			AllCh_name_good = char(join(erase(splitlines(erase(strtrim(AllCh_name_bad),'and')), ' '))); % one single line
			% loop through each frequency in control.favCh
			for i_ch = 1:length(control.favFreq{isensor})
				if contains(AllCh_name_good, control.favFreq{isensor}{i_ch})
					% Locate the characters [name(s) of frequency(ies)]
					oneCh_name_idx_start = strfind(AllCh_name_good,control.favFreq{isensor}{i_ch});
					idx_one_Ch = cellfun(@str2num, extractBetween(AllCh_name_good,oneCh_name_idx_start-2,oneCh_name_idx_start-2));   
					% Save
					Swath_used = [Swath_used,string(file_info.Groups(i_sw).Name)]; %file_info.Groups(i_sw).Name is Character. It is converted to strings in this case. E.g., '/S1' is converted to "/S1"
					ChName_perSwath = [ChName_perSwath, string(control.favFreq{isensor}{i_ch})];
					ChIdx_perSwath = [ChIdx_perSwath,idx_one_Ch];
				else
					continue;
				end       
			end
		end

	 % NC file
	 elseif contains(filext,"nc")
		% loop through each frequency in control.favCh
	    for i_ch = 1:length(control.favFreq{isensor})			
			if contains(control.favFreq{isensor}{i_ch}, 'fcdr_tb19v')
				Swath_used = [Swath_used, string('lores')];
				ChName_perSwath = [ChName_perSwath, string('fcdr_tb19v')];
				ChIdx_perSwath = [ChIdx_perSwath,NaN];
			elseif contains(control.favFreq{isensor}{i_ch}, 'fcdr_tb85v')
				Swath_used = [Swath_used, string('hires')];
				ChName_perSwath = [ChName_perSwath, string('fcdr_tb85v')];
				ChIdx_perSwath = [ChIdx_perSwath,NaN];
			end 
		end

	end

    % Sanity Check
    if length(Swath_used) ~= length(ChName_perSwath) || length(Swath_used) ~= length(ChIdx_perSwath)
    	disp('Error getting used swathes and channels!');
    end	

end

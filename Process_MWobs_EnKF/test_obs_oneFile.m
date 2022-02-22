% -------------- Control variables ----------------------
control = struct;
Tb_file = '/expanse/lustre/projects/pen116/zuy121/Pro2_PSU_MW/Obs/MW_used/Irma2ndRI/DAt201709031000_1C.F17.SSMIS.HDF5';

Swath_used = {'/S1','/S3'};
ChIdx_perSwath = [1,4];
ChName_perSwath = {'19.35GHzV-Pol','183.31+/-6.6GHzH-Pol'};
if_swth_good = [1,1];
DAtime_perSwath = {'201709031000','201709031000'};
loc_sotrm_DAtime = {[18.000,-47.500],[18.000,-47.500]};




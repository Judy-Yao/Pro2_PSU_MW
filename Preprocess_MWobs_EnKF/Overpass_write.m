function [] = Overpass_write(iTb,istorm,Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,Tb_overpass,control) 
	% Preallocating memory
	sat_name = cell(size(Tb_overpass));
	op_lat = cell(size(Tb_overpass));
	op_lon = cell(size(Tb_overpass));
	op_Tb = cell(size(Tb_overpass));
	
	op_scan = cell(size(Tb_overpass));
	op_zenith = cell(size(Tb_overpass));
	op_Sat_lat = cell(size(Tb_overpass));
	op_Sat_lon = cell(size(Tb_overpass));
	op_Sat_alt = cell(size(Tb_overpass));
	op_azimuth = cell(size(Tb_overpass));
	op_Fov_crossTrack = cell(size(Tb_overpass));
	op_Fov_alongTrack = cell(size(Tb_overpass));

	op_times = cell(size(Tb_overpass));
	op_chNum = cell(size(Tb_overpass));
	op_ROI_other = cell(size(Tb_overpass));
	op_ROI_hydro = cell(size(Tb_overpass));
	op_ObsErr = cell(size(Tb_overpass));

	% Loop through each overpass Tb file and read variables with different files into cells
	for io = 1:length(Tb_overpass)
        disp(["    File: " + Tb_overpass(io)]);
        disp('      Frequencies of interest are:');
        for item = 1:length(ChName_all{iTb(io)})
            disp(["        " + ChName_all{iTb(io)}(item)]);
        end

		[filepath,filename,filext] = fileparts(Tb_overpass(io));
		[sat_name{io},op_lat{io},op_lon{io},op_Tb{io},op_Sat_lat{io},op_Sat_lon{io},op_Sat_alt{io},op_azimuth{io},op_scan{io},op_zenith{io},op_Fov_crossTrack{io},op_Fov_alongTrack{io},op_times{io},op_chNum{io},op_ROI_hydro{io},op_ROI_other{io},op_ObsErr{io}] = ProduceforEnKF(iTb(io),Swath_used,ChIdx_all,ChName_all,DAtime_all,loc_DAtime_all,Tb_overpass(io),control);
    end

	% Gather variables of different Tb files with the same ROI combination into the same cell
	myLat = cell(size(control.roi_oh));
    myLon = cell(size(control.roi_oh));
    myTb = cell(size(control.roi_oh));

	mySat_name = cell(size(control.roi_oh));
    mySat_lat = cell(size(control.roi_oh));
    mySat_lon = cell(size(control.roi_oh));
    mySat_alt = cell(size(control.roi_oh));
    myazimuth = cell(size(control.roi_oh)); 
    myScan_angle = cell(size(control.roi_oh));
    myZenith_angle = cell(size(control.roi_oh));
    myFov_crossTrack = cell(size(control.roi_oh));
    myFov_alongTrack = cell(size(control.roi_oh));

    myTimes = cell(size(control.roi_oh));
    myChNum = cell(size(control.roi_oh));
    myRoi_hydro = cell(size(control.roi_oh));
    myRoi_otherVars = cell(size(control.roi_oh));
    myObsErr = cell(size(control.roi_oh));	

	for ir = 1:length(control.roi_oh)
		tem_lat = []; tem_lon = []; tem_Tb = []; tem_scan = []; tem_zenith = [];
		tem_Sat_lat = []; tem_Sat_lon = []; tem_Sat_alt = [];
		tem_azimuth = []; tem_crossTrack = []; tem_alongTrack = [];
		tem_times = []; tem_chNum =[ ]; tem_ROI_other = []; tem_ROI_hydro = []; tem_ObsErr = []; tem_Sat_name = [];
		for io = 1:length(Tb_overpass)
			tem_lat = [tem_lat;op_lat{io}{ir}]; tem_lon = [tem_lon;op_lon{io}{ir}]; tem_Tb = [tem_Tb;op_Tb{io}{ir}]; 
			tem_scan = [tem_scan;op_scan{io}{ir}]; tem_zenith = [tem_zenith;op_zenith{io}{ir}];			
			tem_Sat_lat = [tem_Sat_lat;op_Sat_lat{io}{ir}]; tem_Sat_lon = [tem_Sat_lon; op_Sat_lon{io}{ir}];
			tem_Sat_alt = [tem_Sat_alt; op_Sat_alt{io}{ir}]; tem_azimuth = [tem_azimuth;op_azimuth{io}{ir}];
            tem_crossTrack = [tem_crossTrack; op_Fov_crossTrack{io}{ir}];
			tem_alongTrack = [tem_alongTrack;op_Fov_alongTrack{io}{ir}];			
            tem_times = [tem_times; op_times{io}{ir}];
			tem_chNum = [tem_chNum; op_chNum{io}{ir}]; tem_ROI_other = [tem_ROI_other; op_ROI_other{io}{ir}];
			tem_ROI_hydro = [tem_ROI_hydro; op_ROI_hydro{io}{ir}]; tem_ObsErr = [tem_ObsErr;op_ObsErr{io}{ir}];
			Sat_name_copies = repmat(sat_name{io},[length(op_Tb{io}{ir}),1]);
			tem_Sat_name =[tem_Sat_name;Sat_name_copies];
		end
		% randomize
		randOrder = randperm(length(tem_Tb));
        myLat{ir} = tem_lat(randOrder); myLon{ir} = tem_lon(randOrder); myTb{ir} = tem_Tb(randOrder);
		mySat_lat{ir} = tem_Sat_lat(randOrder); mySat_lon{ir} = tem_Sat_lon(randOrder);
		mySat_alt{ir} = tem_Sat_alt(randOrder); myazimuth{ir} = tem_Sat_alt(randOrder);
		myScan_angle{ir} = tem_scan(randOrder); myZenith_angle{ir} = tem_zenith(randOrder);
        myFov_crossTrack{ir} = tem_crossTrack(randOrder); myFov_alongTrack{ir} = tem_alongTrack(randOrder);
		myTimes{ir} = tem_times(randOrder); myChNum{ir} = tem_chNum(randOrder); myRoi_hydro{ir} = tem_ROI_hydro(randOrder);
		myRoi_otherVars{ir} = tem_ROI_other(randOrder); myObsErr{ir} = tem_ObsErr(randOrder);
		mySat_name{ir} = tem_Sat_name(randOrder);
	end

	% ---------------------------------------------------------------------
    % ---- Write MW records to a file
    % ---------------------------------------------------------------------

    % create text file name
    filename = strcat(control.output_dir,control.storm_phase(istorm),'/microwave_d03_',DAtime_all{iTb(1)}(1),'_so'); % same DA time
	disp("    Output processed MW obs file: "+filename);
    formatSpec = '%12s%16s%12i%12.3f%12.3f%12.3f%12i%12i%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n';
    fileID = fopen(filename,'w');
    % reshape values into columns
    out_lat = cat(1, myLat{:}); out_lon = cat(1, myLon{:}); out_Tb = cat(1, myTb{:});
    out_sat_lat = cat(1, mySat_lat{:}); out_sat_lon = cat(1, mySat_lon{:}); out_sat_alt = cat(1, mySat_alt{:});
    out_azimuth = cat(1, myazimuth{:}); out_Scan_angle = cat(1, myScan_angle{:}); out_Zenith_angle = cat(1, myZenith_angle{:});
    out_Zenith_angle = cat(1, myZenith_angle{:}); out_Fov_crossTrack = cat(1, myFov_crossTrack{:}); out_Fov_alongTrack = cat(1, myFov_alongTrack{:});
    out_times = cat(1, myTimes{:}); out_chNum = cat(1, myChNum{:}); out_obsErr = cat(1, myObsErr{:});
    out_ROI_other = cat(1, myRoi_otherVars{:}); out_ROI_hydro = cat(1, myRoi_hydro{:}); out_Sat_name = cat(1, mySat_name{:});
    
    % write values to a file
    for rd = 1:length(out_Tb)
        if (strcmp(out_Sat_name(rd),'gmi_gpm'))
            if (out_chNum(rd) < 9)
                out_Sat_name(rd) = 'gmi_gpm_lf';
            else
                out_Sat_name(rd) = 'gmi_gpm_hf';
            end
        end
        fprintf(fileID, formatSpec, ...
                out_times(rd), out_Sat_name(rd), out_chNum(rd),...
                out_lat(rd), out_lon(rd), out_Tb(rd), ...
                out_ROI_hydro(rd), out_ROI_other(rd), out_obsErr(rd), ...
                out_Fov_crossTrack(rd), out_Fov_alongTrack(rd),...
                out_Scan_angle(rd), out_Zenith_angle(rd), out_azimuth(rd),...
                out_sat_lat(rd), out_sat_lon(rd), out_sat_alt(rd));
    end
  
    fclose(fileID);


end

% This script reads satellite infrared Tbs observations and gets it geographyically plotted.

% control variables
control = struct;

control.obs_dir = '../../raw_Obs/GOESR_IR/';
control.plot_dir = '../../Visual/Obs/GOESR_IR/';
%control.storm_phase = {'Irma_1stRI','Irma_2ndRI','Jose_RI','Maria_RI'};
control.storm_phase = {'MARIA',};


% Operation for each storm phase
for istorm = 1:length(control.storm_phase)
    switch(control.storm_phase{istorm})  
        case 'Irma_1stRI'
            min_lat = 5;
            max_lat = 20;
            min_lon = -35;
            max_lon = -20;
        case 'Irma_2ndRI'
            min_lat = 10;
            max_lat = 30;
            min_lon = -60;
            max_lon = -40;
        case 'Jose_RI'
            min_lat = 5;
            max_lat = 20;
            min_lon = -50;
            max_lon = -35;
        case 'MARIA'
            min_lat = 0;
            max_lat = 10;
            min_lon = -30;
            max_lon = -15; 
    end
    
    % Operation for each observation file
    Tb_dir = [control.obs_dir, control.storm_phase{istorm},'/*'];
    Tb_files = regexprep(ls(Tb_dir),'\n$', ''); 
    Tb_files = regexp(Tb_files,'\n','split');

    for i = 1:length(Tb_files)
        Tb_file = Tb_files{i};
        disp(Tb_file);
        % subroutine to plot and save plots
        IR_vis(istorm, min_lat, max_lat, min_lon, max_lon, Tb_file, control);
    end
end


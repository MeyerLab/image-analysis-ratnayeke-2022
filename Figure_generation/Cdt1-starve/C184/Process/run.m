%% Rename live images
% manualwells = {};
% ORIGINAL_DIR = {'E:\4TB9\Data\C184-live-nikon\Raw\Timelapse_1\',
%     'E:\4TB9\Data\C184-live-nikon\Raw\Timelapse_2\',
%     'E:\4TB9\Data\C184-live-nikon\Raw\Timelapse_3\'};
% DESTINATION_DIR = 'E:\4TB9\Data\C184-live-nikon\Raw_rename';
% CHANNELS = {'CFP','YFP','RFP'};
% workers = 2;
% delete(gcp('nocreate'))
% parpool(workers);
% nikon_rename(ORIGINAL_DIR,DESTINATION_DIR,CHANNELS,manualwells, workers)
% delete(gcp('nocreate'))

%% Rename fixed images
% manualwells = {};
% ORIGINAL_DIR = {'E:\Nikon\HCA_Live_fixed\20190524_215451_884\'};
% DESTINATION_DIR = 'E:\4TB9\Data\C184-IF1-nikon\Raw';
% CHANNELS = {'DAPI','FarRed','YFP'};
% workers = 2;
% delete(gcp('nocreate'))
% parpool(workers);
% nikon_rename(ORIGINAL_DIR,DESTINATION_DIR,CHANNELS,manualwells, workers)
% delete(gcp('nocreate'))
% 
%% Rename fixed empty wells images
% manualwells = {'2_1','3_1','4_1','5_1'};
% ORIGINAL_DIR = {'J:\C185_PCNA_45_fast15_IF\20190609_164254_275\'};
% DESTINATION_DIR = 'E:\4TB9\Data\C184-IF1-nikon\Raw';
% CHANNELS = {'DAPI','YFP','FarRed'};
% workers = 2;
% delete(gcp('nocreate'))
% parpool(workers);
% nikon_rename(ORIGINAL_DIR,DESTINATION_DIR,CHANNELS,manualwells, workers)
% delete(gcp('nocreate'))

%% Calculate Bias live
% s.experimentpath='E:\4TB9\Data\C184-live-nikon\';
% s.imagepath='Raw_rename\';
% s.shadingpath='';
% s.nucradius=12;%12;
% s.names={
%     'CFP';
%     'YFP';
%     'RFP';
%     };
% s.row_mat = [2:7];
% s.col_mat = [2:11];
% s.site_mat = [1:4];
% s.frame_mat=[1 50 100];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 1; %1 to mask nuclei
% s.adaptive = [0 0 0];
% s.dilate = {.5,.5,.5}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 15;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0);

%% Calculate Bias IF
% s.experimentpath='E:\4TB9\Data\C184-IF1-nikon\';
% s.imagepath='Raw\';
% s.shadingpath='';
% s.nucradius=24;%12;
% s.names={
%     'DAPI';
%     'YFP';
%     'FarRed';
%     };
% s.row_mat = [2 3 4 5];
% s.col_mat = [1];
% s.site_mat = [1:16];
% s.frame_mat=[1];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 0; %1 to mask nuclei
% s.adaptive = [0 0 0];
% s.dilate = {.5,.5,.5}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 15;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0);

%% Process
Parallelwell
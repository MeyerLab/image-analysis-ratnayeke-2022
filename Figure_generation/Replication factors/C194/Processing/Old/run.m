%% Calculate Bias live
% clear s
% s.biaspath='E:\Nikon\C194-live\Bias';
% s.imagepath='E:\Nikon\C194-live\Raw';
% s.shadingpath='';
% s.nucradius=12;%12;
% s.names={
%     'CFP';
%     'YFP';
%     'RFP';
%     };
% s.row_mat = [2:7];
% s.col_mat = [2:11];
% s.site_mat = [1:9];
% s.frame_mat=[10];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 1; %1 to mask nuclei
% s.adaptive = [1 1 0];
% s.dilate = {.75,1.25,.5}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 9;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0);

%% Calculate Bias IF
%clear s
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
%Parallelwell
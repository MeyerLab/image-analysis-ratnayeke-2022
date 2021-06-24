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
% s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% s.microscope = 'nikon';

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

% %% Calculate Bias IF1
% clear s
% s.biaspath='E:\Nikon\C194_IF1\Bias';
% s.imagepath='E:\Nikon\C194_IF1\Raw';
% s.shadingpath='';
% s.nucradius=24;%12;
% s.names={
%     'CFP';
%     'YFP';
%     };
% s.formatCode = 'Well%2$s_Point%2$s_%3$04d_*.tiff';
% s.microscope = 'nikon';
% s.row_mat = [2:5];
% s.col_mat = [1];
% s.site_mat = [1:36];
% s.frame_mat=[1];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 0; %1 to mask nuclei
% s.adaptive = [1 1 0];
% s.dilate = {.75,1.25,.5}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 15;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0);
% 
% %% Calculate Bias IF1
% clear s
% s.biaspath='E:\Nikon\C194_IF1\Bias';
% s.imagepath='E:\Nikon\C194_IF1\Raw';
% s.shadingpath='';
% s.nucradius=24;%12;
% s.names={
%     'CFP';
%     'YFP';
%     };
% s.formatCode = 'Well%2$s_Point%2$s_%3$04d_*.tiff';
% s.microscope = 'nikon';
% s.row_mat = [2:5];
% s.col_mat = [1];
% s.site_mat = [1:36];
% s.frame_mat=[1];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 0; %1 to mask nuclei
% s.adaptive = [1 1 0];
% s.dilate = {.75,1.25,.5}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 15;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0);
% 
% %% Calculate Bias IF3
% clear s
% s.biaspath='E:\Nikon\C194_IF3\Bias';
% s.imagepath='E:\Nikon\C194_IF3\Raw';
% s.shadingpath='';
% s.nucradius=24;%12;
% s.names={
%     'DAPI';
%     'YFP';
%     'Cy5';
%     };
% s.formatCode = 'Well%2$s_Point%2$s_%3$04d_*.tiff';
% s.microscope = 'nikon';
% s.row_mat = [2:5];
% s.col_mat = [1];
% s.site_mat = [1:36];
% s.frame_mat=[1];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 0; %1 to mask nuclei
% s.adaptive = [1 1 0];
% s.dilate = {.75,1.25,.5}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 15;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0);

%% Prep processing
%%% Paths
settings.data_path = 'F:\Data\C-Cdt1\C194-live\Data\';
settings.IF_imagesessions = {'G:\Nikon\C194_IF3\','G:\Nikon\C194_IF2\', 'G:\Nikon\C194_IF1\'};
settings.live_imagepath = 'G:\Nikon\C194-live\Raw\';
settings.bgcmospath = '';
settings.crop_save = 'G:\Nikon\C194-IFcrop\';


%%% General parameters
settings.microscope = 'nikon';
settings.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.formatCodeFixed = 'Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.magnification = 20; %10 = 10x or 20 = 20x
settings.binsize = 1; %1 = bin 1 or 2 = bin 2
settings.postbin = .5; %0 for no bin, number for scaling factor;
settings.match10x20x = 1;
settings.numsites = [9];
settings.scaleLive = 1;
settings.signals = {{'DAPI_','YFP_','Cy5_'}, {'CFP_','YFP_'}, {'CFP_','YFP_'}};
settings.nucLive = 'CFP_';
settings.manualLiveFrame = [];

%%% Quantification parameters
settings.bias = {[1 1 1], [1 1], [1 1]};
settings.biasall = {[0 0 0], [0 0], [0 0]};
settings.sigblur = {[0 0 0], [0 0], [0 0]};
settings.signal_foreground = {[0 0 0], [0 0], [0 0]};
settings.bleedthrough = {[0 0 0], [0 0], [0 0]};
settings.bleedthroughslope = {};
settings.bleedthroughoff = {};
settings.ringcalc = {[1 1 1], [0 1], [0 1]};
settings.ringthresh = {[0 0 0 ], [0 0] , [0 0]};
settings.punctacalc = {[0 0 0], [0 0], [0 0]};
settings.punctaThresh = {{[],[],[]},{}};
settings.punctatopsize = 2;
settings.localbg = {[1 1 1], [0 1], [0 1]};
settings.minringsize = 100;
settings.bgcmoscorrection = 1; %1 = correct for shading; 0 = dont correct for shading;
settings.bgsubmethod = {{'global nuclear','global nuclear', 'global nuclear'}, ...
    {'global nuclear', 'none'},{'global nuclear', 'global nuclear'}}; %'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25 25], [25 25], [25 25]};
settings.frameIF=1;
 
%%% Segmentation parameters
settings.segmethod = {'concavity','thresh','thresh'}; %'log' or 'single' or 'double'
settings.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100; %100
settings.boulderarea = 1500; %1500 % less restrictive for C152 for >4N cells
settings.blobthreshold = -0.03;
settings.blurradius = 3;
settings.soliditythresh = 0.80;
settings.compression = 4;
settings.badFrameCheck = .25; %0 for no check

%%% Tracking parameters
settings.distthresh = 3*settings.nucr;
settings.arealowthresh = -.5;
settings.areahighthresh = .4;

%Timelapse_addIF(settings,5,2,1,0)


%% Run timelapse
rows = [2:7];
cols = [2:11];
sites = [1:9];
Parallelwell({@Timelapse_addIF},settings,rows,cols,sites)

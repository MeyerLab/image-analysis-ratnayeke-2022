%% Rename live sequence
% nikon_rename('K:\Nikon\C222-live\Raw','K:\Nikon\C222-live\Raw_renamed','COPY',false)
% nikon_rename('K:\Nikon\C222-IF1\20191210_113601_509','K:\Nikon\C222-IF1\Raw','COPY',false)
% nikon_rename('K:\Nikon\C222-IF2\20191210_174343_144','K:\Nikon\C222-IF2\Raw','COPY',false)
%% Calculate Bias live
clear p
p.biaspath='K:\Nikon\C222-live\Bias';  % Output path for bias
p.imagepath='K:\Nikon\C222-live\Raw';  % Image directory. If nikon, single folder for all sites
p.shadingpath=''; % cmos offset file for IX micro
p.names={
    'CFP';
    'YFP';
    'RFP';
    };
p.row_mat = [2 :7];
p.col_mat = [2:11];
p.site_mat = [1:9];
p.frame_mat=[10];
p.biasAll = 0; % Save averaged bias across sites
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
p.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
p.microscope = 'nikon';

%%% Settings
p.maskforeground = 1;   % 1 to mask nuclei
p.nucradius=12;%12;
p.blur_radius = 5;
p.adaptive = [0 0 0];   % use add channels together to segment foreground
p.dilate = {.75,.75,1};   % multiples of nucradius, default .5, to expand mask
p.foreground_calc = 0;  % 1 for foreground calc (average foreground signal)
p.method = 'block';     % or pixel
p.blocknum= 11;
p.prctilethresh = 50;
p.compress = .25;
p.prctile_thresh=[0 100]; % pixel remove outliers
p.sigma = 25;             % pixel smooth window


% biasCalc(p,0);


%% Prep live processing 
%% Paths
s.experiment_name='C222-live';
s.image_drive = 'K:\Nikon\C222-live\';
s.savepath=fullfile('F:\Data\C-Cdt1\',s.experiment_name,'Data');
s.imagepath=fullfile(s.image_drive,'Raw'); %single directory
s.biaspath=fullfile(s.image_drive,'Bias');
s.maskpath=fullfile(s.image_drive,'Mask');
s.bgcmospath='';
s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
s.microscope = 'nikon';

%%% General parameters
s.startFrame = 1; 
%s.endFrame = 110;
s.endFrame = getMaxFrame(fullfile(s.imagepath,'2_2_1'),'Time(\d+)_');
s.magnification=10; %10=10x or 20=20x
s.binsize=1; %1=bin 1 or 2=bin 2
s.postbin=0;
s.signals={'CFP_','YFP_','RFP_'};
s.maskwrite = 1; %1 = save an image with the nucmask; 0 = dont save an image with the nucmask
s.maskname = 'nucedge_';
s.register = 0;
s.register_exception = [];

%%% Quantification parameters
s.bgcmoscorrection = 1;
s.bias = [1 1 1];
s.signal_foreground = [0 0 0];
s.bgsubmethod = {'global nuclear','global nuclear','global nuclear'}; %Options:'global nuclear','global cyto','tophat','semi-local nuclear', 'none'
s.compression = 4;
s.bgprctile = [25 25 25];
s.sigblur = [3 3 3];
s.localbgmeasure = [0 0 0];
s.ringcalc = [0 0 0];
s.ringthresh = [0 0 0];
s.punctacalc = 0;
s.punctaThresh = [125 150 200];
s.varThresh = [75 100 125];

%%% Segmentation parameters
s.firstsegmethod = 'log contour'; %Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
s.secondsegmethod = 'log contour'; %Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
s.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
s.blurradius  =  7; %10x: 3
s.soliditythresh = 0.9;
s.debrisarea = 75; % 10x = 100 % 150 if no mitosis
s.boulderarea = 1500; % 10x = 1500
s.blobthreshold = -0.03;

%%% Tracking parameters

s.maxjump = s.nucr*3;
s.masschangethreshold = 0.30;
s.areachangethreshold = 0.60;
s.daughtervariance = 0.10;

Timelapse(s,2,2,8,0)
% % 
%% Run timelapse
rows = [2:7];
cols = [2:11];
sites = [1:9];
% Parallelwell({@Timelapse},s,rows,cols,sites)



%% Calculate Bias IF1 (used from C216)
clear s
s.biaspath='K:\Nikon\C222-IF1\Bias';
s.imagepath='K:\Nikon\C222-IF1\Raw';
s.shadingpath='';
s.nucradius=24;%12;
s.names={
    'DAPI';
    'YFP';
    'FarRed';

    };
s.biasAll = 1;
s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
s.microscope = 'nikon';
s.row_mat = [2:5];
s.col_mat = [1];
s.site_mat = [1:36];
s.frame_mat=[1];
%%% Settings
s.blur_radius = 5;
s.method = 'block'; % or pixel
s.maskforeground = 0; %1 to mask nuclei
s.adaptive = [0 0 0];
s.dilate = {.75,1.25,.5}; %nucradius/2
s.foreground_calc = 0; % 1 for foreground calc
s.blocknum= 11;
s.prctilethresh = 50;
s.compress = .25;
s.prctile_thresh=[0 100]; %pixel
s.sigma = 25; % pixel

% biasCalc(s,0);

  %% Prep processing
%%% Paths
settings.data_path = 'F:\Data\C-Cdt1\C222-live\Data\';
settings.IF_imagesessions = {'K:\Nikon\C222-IF1\', 'K:\Nikon\C222-IF2\'};
settings.live_imagepath = 'K:\Nikon\C222-live\Raw';
settings.bgcmospath = '';
settings.crop_save = 'K:\Nikon\C222-IFcrop\';
settings.mask_save = 'K:\Nikon\C222-IFcrop\';           % Output mask, leave empty if no save


%%% General parameters
settings.microscope = 'nikon';
settings.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.magnification = 20; %10 = 10x or 20 = 20x
settings.binsize = 1; %1 = bin 1 or 2 = bin 2
settings.postbin = .5; %0 for no bin, number for scaling factor;
settings.match10x20x = 1;
settings.numsites = 9;
settings.scaleLive = 1;
settings.maxjit = 500;
settings.signals = {{'DAPI_','YFP_','FarRed_'},{'DAPI_','YFP_'}};
settings.nucLive = 'CFP_';
settings.manualLiveFrame = [];

%%% Quantification parameters
settings.bias = {[1 1 1], [1 1]};
settings.biasall = {[0 0 1], [0 0]};
settings.sigblur = {[0 0 0], [0 0]};
settings.signal_foreground = {[0 0 0], [0 0]};
settings.bleedthrough = {[0 0 0], [0 0]};
settings.bleedthroughslope = {};
settings.bleedthroughoff = {};
settings.ringcalc = {[1 1 1], [0 1]};
settings.ringthresh = {[0 0 0 ], [0 0]};
settings.punctacalc = {[0 0 0], [0 0]};
settings.punctaThresh = {{[],[],[]},{}};
settings.punctatopsize = 2;
settings.localbg = {[1 1 1], [1 1]};
settings.minringsize = 100;
settings.bgcmoscorrection = 1; %1 = correct for shading; 0 = dont correct for shading;
settings.bgsubmethod = {{'global nuclear','global nuclear', 'global nuclear'},{'global nuclear', 'global cyto'}}; %'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25 25], [25 25]};
settings.frameIF=1;

%%% Segmentation parameters
settings.segmethod = {'concavity', 'thresh'}; %'log' or 'single' or 'double'
settings.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100; %100
settings.boulderarea = 1500; %1500 % less restrictive for C152 for >4N cells
settings.blobthreshold = -0.03;
settings.blurradius = 3;
settings.soliditythresh = 0.50;
settings.compression = 4;
settings.badFrameCheck = .25; %0 for no check

%%% Tracking parameters
settings.distthresh = 3*settings.nucr;
settings.arealowthresh = -.5;
settings.areahighthresh = .75;

% Timelapse_addIF(settings,2,7,5,0)


%% Run add IF
rows = [2:7];
cols = [2:11];
sites = [1:9];
% Parallelwell({@Timelapse_addIF},settings,rows,cols,sites)

%% Run Fixed
settings.IF_imagesessions = {'G:\8TB6\Data\Nikon\C222-IF1\', 'G:\8TB6\Data\Nikon\C222-IF2\'};
settings.crop_save = '';
settings.mask_save = ''; 
settings.maskname = 'mask';

% Fixed_IF(settings,2,7,5,0)

rows = [2:3];
cols = [2:11];
sites = [1:36];
Parallelwell({@Fixed_IF},settings,rows,cols,sites)


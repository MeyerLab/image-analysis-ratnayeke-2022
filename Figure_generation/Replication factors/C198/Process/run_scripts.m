% IX_Micro_rename('H:\Nalin\201901003-C198-live\10696','H:\8TB4\Data\C198-live\Raw',{'CFP','YFP','RFP'},'COPY',false);
% IX_Micro_rename('H:\Nalin\201901003-C198-IF1\2019-10-05\10699','H:\8TB4\Data\C198-live-IF1\Raw',{'DAPI','YFP','FarRed'});
%IX_Micro_rename('H:\Nalin\201901003-C198-IF2\2019-10-07\10701','H:\8TB4\Data\C198-live-IF2\Raw',{'DAPI','YFP','FarRed'});
%% Calculate Bias live
% clear s
% s.biaspath='H:\8TB4\Data\C198-live\Bias';
% s.imagepath='H:\8TB4\Data\C198-live\Raw';
% s.shadingpath='F:\MATLAB\BGimages\IX_Liu\cmosoffset_bin1.mat';
% s.nucradius=12;%12;
% s.names={
%     'CFP';
%     'YFP';
%     'RFP';
%     };
% s.formatCode = '';
% s.microscope = 'IXM';
% s.row_mat = [2:6];
% s.col_mat = [2:11];
% s.site_mat = [1:4];
% s.frame_mat=[1];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 1; %1 to mask nuclei
% s.adaptive = [0 0 0];
% s.dilate = {.75,1,1}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 15;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0)

%% Calculate Bias IF1
% clear s
% s.biaspath='H:\8TB4\Data\C198-IF1\Bias';
% s.imagepath='H:\8TB4\Data\C198-IF1\Raw';
% s.shadingpath='F:\MATLAB\BGimages\IX_Liu\cmosoffset_bin1.mat';
% s.nucradius=12;%12;
% s.names={
%     'DAPI';
%     'YFP';
%     'FarRed';
%     };
% s.formatCode = '';
% s.microscope = 'IXM';
% s.row_mat = [2:5];
% s.col_mat = [1];
% s.site_mat = [1:16];
% s.frame_mat=[1];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 0; %1 to mask nuclei
% s.adaptive = [0 0 0];
% s.dilate = {.75,1,1}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 15;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0);

%% Prep processing
% clear s
% %%% Paths
% s.experiment_name='C198-live';
% s.image_drive = 'H:\8TB4\Data\C198-live\';
% s.savepath=fullfile('F:\Data\C-Cdt1\',s.experiment_name,'Data');
% s.imagepath={fullfile(s.image_drive,'Raw')}; %single directory
% s.biaspath=fullfile(s.image_drive,'Bias');
% s.maskpath=fullfile(s.image_drive,'Mask');
% s.bgcmospath='C:\Users\Meyerlab\Documents\GitHub\cell-cycle-tracking\BGimages\IX_Liu\cmosoffset_bin1.mat';
% s.formatCode = '';
% s.microscope = 'IXM';
% 
% %%% General parameters
% s.Frames = [1 20];   %Frames to analyze [start_frame, end_first etc, end_frame]
% s.startFrame = 1; 
% s.magnification=10; %10=10x or 20=20x
% s.binsize=1; %1=bin 1 or 2=bin 2
% s.signals={'CFP_','YFP_','RFP_'};
% s.maskwrite = 1; %1 = save an image with the nucmask; 0 = dont save an image with the nucmask
% s.maskname = 'nucedge_';
% s.register = 0;
% s.register_exception = [];
% 
% %%% Quantification parameters
% s.bgcmoscorrection = 1;
% s.bias = [1 1 1];
% s.signal_foreground = [0 0 0];
% s.bgsubmethod = {'global nuclear','global nuclear','global nuclear'}; %Options:'global nuclear','global cyto','tophat','semi-local nuclear', 'none'
% s.compression = 4;
% s.bgprctile = [25 25 25];
% s.sigblur = [3 3 3];
% s.localbgmeasure = [0 0 0];
% s.ringcalc = [0 0 0];
% s.ringthresh = [0 50 50];
% s.punctacalc = 0;
% s.punctaThresh = [125 150 200];
% s.varThresh = [75 100 125];
% 
% %%% Segmentation parameters
% s.firstsegmethod = 'concavity'; %Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
% s.secondsegmethod = 'concavity'; %Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
% s.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
% s.blurradius  =  3; %10x: 3
% s.soliditythresh = 0.8;
% s.debrisarea = 100; % 10x = 100 % 150 if no mitosis
% s.boulderarea = 1500; % 10x = 1500
% s.blobthreshold = -0.02;
% 
% %%% Tracking parameters
% s.maxjump = s.nucr*3;
% s.masschangethreshold = 0.30;
% s.areachangethreshold = 0.60;
% s.daughtervariance = 0.10;
% 
% Timelapse(s,2,2,3,0)
% 
% %% Run timelapse
% rows = [2:7];
% cols = [2:11];
% sites = [1:4];
% Parallelwell({@Timelapse},s,rows,cols,sites)

%% Prep processing
%%% Paths
settings.data_path = 'F:\Data\C-Cdt1\C198-live\Data\';
settings.IF_imagesessions = {'H:\8TB4\Data\IXM\C198-IF1\','H:\8TB4\Data\IXM\C198-IF2\'};
settings.live_imagepath = 'H:\8TB4\Data\IXM\C198-live\Raw\';
settings.bgcmospath = 'C:\Users\nalin\Documents\GitHub\cell-cycle-tracking\BGimages\IX_Liu\cmosoffset_bin1.mat';
settings.crop_save = 'H:\8TB4\Data\IXM\C198-IFcrop\';


%%% General parameters
settings.microscope = 'IXM';
settings.formatCodeLive = '';
settings.formatCodeFixed = '';
settings.magnification = 20; %10 = 10x or 20 = 20x
settings.binsize = 1; %1 = bin 1 or 2 = bin 2
settings.postbin = 0; %0 for no bin, number for scaling factor;
settings.match10x20x = 1;
settings.numsites = [4];
settings.scaleLive = 2;
settings.signals = {{'DAPI_','YFP_','FarRed_'}, {'DAPI_','YFP_','FarRed_'}};
settings.nucLive = 'CFP_';
settings.manualLiveFrame = [];

%%% Quantification parameters
settings.bias = {[1 1 1], [1 1 1]};
settings.biasall = {[0 0 0], [0 0 0]};
settings.sigblur = {[0 0 0], [0 0 0]};
settings.signal_foreground = {[0 0 0], [0 0 0]};
settings.bleedthrough = {[0 0 0], [0 0 0]};
settings.bleedthroughslope = {};
settings.bleedthroughoff = {};
settings.ringcalc = {[1 1 1], [0 1 1]};
settings.ringthresh = {[0 0 0 ], [0 0 0]};
settings.punctacalc = {[0 0 0], [0 0 0]};
settings.punctaThresh = {{[],[],[]},{}};
settings.punctatopsize = 2;
settings.localbg = {[0 0 0], [0 0 0]};
settings.minringsize = 100;
settings.bgcmoscorrection = 1; %1 = correct for shading; 0 = dont correct for shading;
settings.bgsubmethod = {{'global nuclear','global nuclear', 'none'}, ...
    {'global nuclear', 'global nuclear', 'global nuclear'}}; %'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25 25], [25 25 25]};
settings.frameIF=1;
 
%%% Segmentation parameters
settings.segmethod = {'concavity','thresh'}; %'log' or 'single' or 'double'
settings.nucr = 24; %10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 400; %100
settings.boulderarea = 6000; %1500 % less restrictive for C152 for >4N cells
settings.blobthreshold = -0.03;
settings.blurradius = 3;
settings.soliditythresh = 0.80;
settings.compression = 4;
settings.badFrameCheck = .25; %0 for no check

%%% Tracking parameters
settings.distthresh = 3*settings.nucr;
settings.arealowthresh = -.4;
settings.areahighthresh = .5;

Timelapse_addIF(settings,5,2,1,0)


%% Run timelapse
rows = [2:7];
cols = [2:11];
sites = [1:4];
Parallelwell({@Timelapse_addIF},settings,rows,cols,sites)
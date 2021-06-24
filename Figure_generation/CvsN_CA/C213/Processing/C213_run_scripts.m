%% Rename live sequence
% nikon_rename_time('J:\C214_live\Raw2','J:\C214_live\Raw','COPY',false)

%% Calculate Bias live
clear p
p.biaspath='I:\Nikon\C213-live\Bias';  % Output path for bias
p.imagepath='I:\Nikon\C213-live\Raw';  % Image directory. If nikon, single folder for all sites
p.shadingpath=''; % cmos offset file for IX micro
p.names={
    'CFP';
    'YFP';
    'RFP';
    };
p.row_mat = [2:4];
p.col_mat = [2:11];
p.site_mat = [1:2];
p.frame_mat=[10 50];
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


%biasCalc(p,0);


%% Prep live processing 
%% Paths
s.experiment_name='C213-live';
s.image_drive = 'I:\Nikon\C213-live\';
s.savepath=fullfile('F:\Data\C-Cdt1\',s.experiment_name,'Data');
s.imagepath=fullfile(s.image_drive,'Raw'); %single directory
s.biaspath=fullfile(s.image_drive,'Bias');
s.maskpath=fullfile(s.image_drive,'Mask');
s.bgcmospath='';
s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
s.microscope = 'nikon';

%%% General parameters
s.startFrame = 1; 
%s.endFrame = 80;
s.endFrame = getMaxFrame(s.imagepath,'Time(\d+)_');
s.magnification=10; %10=10x or 20=20x
s.binsize=1; %1=bin 1 or 2=bin 2
s.postbin=0;
s.signals={'CFP_','YFP_','RFP_'};
s.maskwrite = 1; %1 = save an image with the nucmask; 0 = dont save an image with the nucmask
s.maskname = 'nucedge_';
s.register = 0;
s.register_exception = [81:85];

%%% Quantification parameters
s.bgcmoscorrection = 1;
s.bias = [1 1 1];
s.signal_foreground = [0 0 0];
s.bgsubmethod = {'global nuclear','global nuclear','global nuclear'}; %Options:'global nuclear','global cyto','tophat','semi-local nuclear', 'none'
s.compression = 4;
s.bgprctile = [25 25 25];
s.sigblur = [3 3 3];
s.localbgmeasure = [0 0 0];
s.ringcalc = [0 0 1];
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

% Timelapse(s,2,2,1,1)
% % 
%% Run timelapse
rows = [2:4];
cols = [2:11];
sites = [1:2];
Parallelwell({@Timelapse},s,rows,cols,sites)

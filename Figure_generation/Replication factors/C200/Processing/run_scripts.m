
%% Calculate Bias live
% clear p
% p.biaspath='H:\8TB4\Data\Nikon\C200-live\Bias';  % Output path for bias
% p.imagepath='H:\8TB4\Data\Nikon\C200-live\Raw';  % Image directory. If nikon, single folder for all sites
% p.shadingpath=''; % cmos offset file for IX micro
% p.names={
%     'CFP';
%     'YFP';
%     'RFP';
%     };
% p.row_mat = [2:6];
% p.col_mat = [2:11];
% p.site_mat = [1:9];
% p.frame_mat=[10 50];
% p.biasAll = 0; % Save averaged bias across sites
% % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
% p.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% p.microscope = 'nikon';
% 
% %%% Settings
% p.maskforeground = 1;   % 1 to mask nuclei
% p.nucradius=12;%12;
% p.blur_radius = 5;
% p.adaptive = [0 0 0];   % use add channels together to segment foreground
% p.dilate = {.5,.75,.75};   % multiples of nucradius, default .5, to expand mask
% p.foreground_calc = 0;  % 1 for foreground calc (average foreground signal)
% p.method = 'block';     % or pixel
% p.blocknum= 11;
% p.prctilethresh = 50;
% p.compress = .25;
% p.prctile_thresh=[0 100]; % pixel remove outliers
% p.sigma = 25;             % pixel smooth window
% 
% 
% biasCalc(p,0);

%% Prep live processing 
%% Paths
% s.experiment_name='C200-live';
% s.image_drive = 'H:\8TB4\Data\Nikon\C200-live\';
% s.savepath=fullfile('F:\Data\C-Cdt1\',s.experiment_name,'Data');
% s.imagepath=fullfile(s.image_drive,'Raw'); %single directory
% s.biaspath=fullfile(s.image_drive,'Bias');
% s.maskpath=fullfile(s.image_drive,'Mask');
% s.bgcmospath='';
% s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% s.microscope = 'nikon';
% 
% %%% General parameters
% s.startFrame = 1; 
% s.endFrame = 86;
% s.magnification=10; %10=10x or 20=20x
% s.binsize=1; %1=bin 1 or 2=bin 2
% s.postbin=0;
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
% s.ringcalc = [0 1 1];
% s.ringthresh = [0 0 0];
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

% Timelapse(s,6,2,3,0)
% % % 
% %% Run timelapse
% rows = [2:7];
% cols = [2:11];
% sites = [1:9];
%  Parallelwell({@Timelapse},s,rows,cols,sites)

%% Calculate Bias IF2
% clear p
% p.biaspath='H:\8TB4\Data\Nikon\C200-IF2\Bias';  % Output path for bias
% p.imagepath='H:\8TB4\Data\Nikon\C200-IF2\Raw';  % Image directory. If nikon, single folder for all sites
% p.shadingpath=''; % cmos offset file for IX micro
% p.names={
%     'DAPI';
%     'YFP';
%     'FarRed';
%     };
% p.row_mat = [2:5];
% p.col_mat = [1];
% p.site_mat = [1:36];
% p.frame_mat=[1];
% p.biasAll = 0; % Save averaged bias across sites
% % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
% p.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% p.microscope = 'nikon';
% 
% %%% Settings
% p.maskforeground = 0;   % 1 to mask nuclei
% p.nucradius=12;%12;
% p.blur_radius = 5;
% p.adaptive = [0 0 0];   % use add channels together to segment foreground
% p.dilate = {.5,.75,.75};   % multiples of nucradius, default .5, to expand mask
% p.foreground_calc = 0;  % 1 for foreground calc (average foreground signal)
% p.method = 'block';     % or pixel
% p.blocknum= 11;
% p.prctilethresh = 50;
% p.compress = .25;
% p.prctile_thresh=[0 100]; % pixel remove outliers
% p.sigma = 25;             % pixel smooth window


%biasCalc(p,0);

%% Prep IF processing
sIF.data_path = 'F:\Data\C-Cdt1\C200-live\Data\';      % Output location of mat files
sIF.IF_imagesessions = {'H:\8TB4\Data\Nikon\C200-IF1\','H:\8TB4\Data\Nikon\C200-IF2\','H:\8TB4\Data\Nikon\C200-IF3\'};    
                                                            % Cell array containing all imaging sessions to be matched
sIF.live_imagepath = 'H:\8TB4\Data\Nikon\C200-live\Raw\';    % Live imaging raw image path
sIF.bgcmospath = '';
sIF.crop_save = 'H:\8TB4\Data\Nikon\C200-IFcrop\';           % Output cropped images, leave empty if no save
sIF.mask_save = 'H:\8TB4\Data\Nikon\C200-IFcrop\';           % Output mask, leave empty if no save

%%% General parameters
sIF.microscope = 'nikon';
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
sIF.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
sIF.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
sIF.magnification = 20;     % 10 = 10x or 20 = 20x
sIF.binsize = 1;            % 1 = bin 1 or 2 = bin 2
sIF.postbin = 0;           % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
sIF.match10x20x = 1;
sIF.numsites = 9;           % Number of sites imaged for 20x to 10x (double check order is correct)
sIF.scaleLive = 2;          % Scale live to have same pixel size as final IF image
sIF.maxjit = 500;      % maximum jitter between IF images before throwing warning
sIF.signals = {{'DAPI_','YFP_','FarRed_'}, {'DAPI_','YFP_','FarRed_'}, {'DAPI_','FarRed_'}};  % Each imaging session in a cell inside cell array
sIF.nucLive = 'CFP_';       % Label for IX micro live nucleus (assume first tiff in stack for nikon)
sIF.manualLiveFrame = [];   % Manually set last frame for live imaging, otherwise use tracedata if empty

%%% Quantification parameters
sIF.bias = {[1 1 1], [1 1 1], [1 1]};
sIF.biasall = {[0 0 0], [0 0 0], [0 0]};
sIF.sigblur = {[0 0 0], [0 0 0], [0 0]};
sIF.signal_foreground = {[0 0 0], [0 0 0], [0 0]};
sIF.bleedthrough = {[0 0 0], [0 0 0], [0 0]};   % Calculate bleedthrough for channel
sIF.bleedthroughslope = {};              % Cell array of cell arrays
sIF.bleedthroughoff = {};                % Cell array of cell arrays
sIF.ringcalc = {[1 1 1], [0 1 1], [0 1]};
sIF.ringthresh = {[0 0 0 ], [0 0 0], [0 0]};    % Threshold for foreground ring pixels
sIF.punctacalc = {[0 0 0], [0 0 0], [0 0]};     % Analyze puncta in nucleus
sIF.punctaThresh = {{[],[],[]},{}};      % Threshold for tophat filter for puncta
sIF.punctatopsize = 2;                   % Top hat filter size
sIF.localbg = {[0 0 0], [0 0 0], [0 0]};        % Measure local background in box around cell
sIF.minringsize = 100;
sIF.bgcmoscorrection = 1;                % 1 = correct for shading; 0 = dont correct for shading;
sIF.compression = 4;
sIF.bgsubmethod = {{'global nuclear','global nuclear', 'none'}, ...
    {'global nuclear', 'global nuclear', 'global nuclear'}, ...
    {'global nuclear', 'global nuclear'}}; % 'global nuclear','global cyto','tophat','semi-local nuclear'
sIF.bgperctile = {[25 25 25], [25 25 25], [25 25]};  % Set bg as perctile for each channel
sIF.frameIF=1;
 
%%% Segmentation parameters
sIF.segmethod = {'concavity','thresh','thresh'};  % 'log' or 'single' or 'double'
sIF.nucr = 24;                           % 10x bin1 = 12 and 20x bin2 = 12
sIF.debrisarea = 400;                    % 100
sIF.boulderarea = 6000;                  % 1500 
sIF.blobthreshold = -0.03;               % For blobdector 'log' segmentation
sIF.blurradius = 3;                      % Blur for blobdetector
sIF.soliditythresh = 0.80;               % Minimum solidity of nucleus allowed
sIF.badFrameCheck = .25;                 % Min fraction of cells untracked before logging error, 0 for no check

%%% Tracking parameters
sIF.distthresh = 3*sIF.nucr;
sIF.arealowthresh = -.4;
sIF.areahighthresh = .5;

%Timelapse_addIF(sIF,6,2,3,0)
% % % 
% %% Run timelapse
rows = [2:7];
cols = [2:11];
sites = [1:9];
 Parallelwell({@Timelapse_addIF},sIF,rows,cols,sites)

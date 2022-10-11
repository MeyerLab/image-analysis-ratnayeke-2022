cd('F:\Data\C-Cdt1\C229-IF-p27')
%% Rename live sequence
% nikon_rename({'O:\C229_p27_FP\20211109_141751_582'},...
%      'O:\C229_p27_FP\Raw','COPY',false)
nikon_rename({'O:\C229_p27_IF1\20211110_120716_487'},...
     'O:\C229_p27_IF1\Raw','COPY',false)
 nikon_rename({'O:\C229_p27_IF2\20211111_111831_190'},...
     'O:\C229_p27_IF2\Raw','COPY',false)
%% Calculate Bias FP
% clear p
% p.biaspath='L:\D126-live\Bias';  % Output path for bias
% p.imagepath='L:\D126-live\Raw';  % Image directory. If nikon, single folder for all sites
% p.shadingpath=''; % cmos offset file for IX micro
% p.names={
%     'CFP';
%     'YFP';
%     'RFP';
%     };
% p.row_mat = [1:6];
% p.col_mat = [1:12];
% p.site_mat = [1];
% p.frame_mat=[2 50];
% p.biasAll = 1; % Save averaged bias across sites
% % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
% % 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff' if multipoint
% p.formatCode = 'Time%1$05d_Well%2$s*.tiff';
% p.microscope = 'nikon';
% p.nuc_channel = 1;
% 
% %%% Settings
% p.maskforeground = 1;   % 1 to mask nuclei
% p.nucradius=12;%12;
% p.blur_radius = 5;
% p.adaptive = [0 0 0 0];   % use add channels together to segment foreground
% p.dilate = {.75,.75,.75 .75};   % multiples of nucradius, default .5, to expand mask
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
% 


%% Setup IF images
clear settings
settings.data_path = 'F:\Data\C-Cdt1\C229-IF-p27\Data\';      % Output location of mat files
settings.IF_imagesessions = {'O:\C229_p27_IF1\','O:\C229_p27_IF2\','O:\C229_p27_FP\'};    
                                                            % Cell array containing all imaging sessions to be matched
settings.bgcmospath = '';
settings.crop_save = 'O:\C229_p27_IFcrop';           % Output cropped images, leave empty if no save
settings.mask_save = 'O:\C229_p27_IFcrop';           % Output mask, leave empty if no save
settings.maskname = 'mask';
settings.primaryMaskRound = 1;
settings.maskIndex = [1 1 4];

%%% General parameters
settings.microscope = 'nikon';
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
settings.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'; %'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'
settings.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.magnification = 20;     % 10 = 10x or 20 = 20x
settings.binsize = 1;            % 1 = bin 1 or 2 = bin 2
settings.postbin = [.5];           % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
settings.match10x20x = 0;
settings.numsites = 1;           % Number of sites imaged for 20x to 10x (double check order is correct)
settings.scaleLive = 1;          % Scale live to have same pixel size as final IF image
settings.signals = {{'DAPI','YFP','Cy5'},{'DAPI','YFP','Cy5'},{'CFP','YFP','RFP','iRFP'}};  % Each imaging session in a cell inside cell array
settings.nucLive = 'CFP_';       % Label for IX micro live nucleus (assume first tiff in stack for nikon)
settings.liveChannel = 4;         % If nikon use this tiff in stack for nuclear live
settings.manualLiveFrame = [];   % Manually set last frame for live imaging, otherwise use tracedata if empty

%%% Segmentation parameters
settings.segmethod = {'concavity','concavity','log'};  % 'log' or 'single' or 'double'
settings.nucr = 12;                           % 10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100;                    % 100
settings.boulderarea = 1500;                  % 1500 
settings.blobthreshold = -0.03;               % For blobdector 'log' segmentation
settings.blurradius = 3;                      % Blur for blobdetector
settings.soliditythresh = 0.80;               % Minimum solidity of nucleus allowed
settings.badFrameCheck = .25;                 % Min fraction of cells untracked before logging error, 0 for no check
settings.split_mult = 1;
settings.maxjit=100;

%%% Quantification parameters
settings.bias = {[1 1 1 ],[1 1 1 ],[1 1 1 1] };
settings.biasall = {[1 1 1 ],[1 1 1 ],[1 1 1 1] };
settings.sigblur = {[0 0 0],[0 0 0], [0 0 0 0]  };
settings.signal_foreground = { [0 0 0],[0 0 0], [0 0 0 0]   };
settings.bleedthrough = {[0 0 0],[0 0 0], [0 0 0 0]   };   % Calculate bleedthrough for channel
settings.bleedthroughslope = {};              % Cell array of cell arrays
settings.bleedthroughoff = {};                % Cell array of cell arrays
settings.ringcalc = {[0 1 0],[0 1 0], [1 1 1 0]   };
settings.ringthresh = { [0 0 0],[0 0 0], [0 0 0 0]};    % Threshold for foreground ring pixels
settings.punctacalc = { [0 0 0],[0 0 0], [0 0 0 0]};     % Analyze puncta in nucleus
settings.punctaThresh = {{[],[],[]},{}};      % Threshold for tophat filter for puncta
settings.punctatopsize = 2;                   % Top hat filter size
settings.cytopuncta = { [0 0 0],[0 0 0], [0 0 0 0]};
settings.thickenradius = 2*settings.nucr; 

settings.localbg = { [0 0 0],[0 0 0], [0 0 0 0]};        % Measure local background in box around cell
settings.minringsize = 100;
settings.bgcmoscorrection = 1;                % 1 = correct for shading; 0 = dont correct for shading;
settings.compression = 4;
settings.bgsubmethod = {{'global nuclear','global cyto', 'global nuclear'},...
    {'global nuclear','global cyto', 'global nuclear'},...
    {'global nuclear','global cyto','global cyto', 'global nuclear'}}; % 'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25 25],[25 25 25],[25 25 25 25]};  % Set bg as perctile for each channel
settings.frameIF=1;
 
Fixed_IF( settings, 2,5,2,1)


%% Run timelapse
rows = [2:7];
cols = [2:11];
sites = [1:32];
Parallelwell({@Fixed_IF},settings,rows,cols,sites)

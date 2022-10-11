cd('F:\Data\F-G1S\F024-live\Processing')
%% Rename IF 
nikon_rename('N:\F024-live-fixed\20211017_154825_780',...
     'N:\F024-live-fixed\Raw','COPY',false,'SINGLE_SITE',false)
nikon_rename('N:\F024-live-fixed\20211017_160310_752',...
     'N:\F024-live-fixed\Raw','COPY',false,'SINGLE_SITE',false)

%% Calculate Bias IF
clear p
p.biaspath='N:\F024-live-fixed\Bias';  % Output path for bias
p.imagepath='N:\F024-live-fixed\Raw';  % Image directory. If nikon, single folder for all sites
p.shadingpath=''; % cmos offset file for IX micro
p.names={
    'DAPI';
    'YFP';
    'FarRed';
    };
p.row_mat = [2:4];
p.col_mat = [1];
p.site_mat = [1:24];
p.frame_mat=[1];
p.biasAll = 1; % Save averaged bias across sites
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
% 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff' if multipoint
p.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
p.microscope = 'nikon';
p.nuc_channel = 1;

%%% Settings
p.maskforeground = 0;   % 1 to mask nuclei
p.nucradius=12;%12;
p.blur_radius = 5;
p.adaptive = [0 0 0 0];   % use add channels together to segment foreground
p.dilate = {1,.75,.75 .75};   % multiples of nucradius, default .5, to expand mask
p.foreground_calc = 0;  % 1 for foreground calc (average foreground signal)
p.method = 'block';     % or pixel
p.blocknum= 11;
p.prctilethresh = 50;
p.compress = .25;
p.prctile_thresh=[0 100]; % pixel remove outliers
p.sigma = 25;             % pixel smooth window


biasCalc(p,0);

%% Setup IF images
clear settings
settings.data_path = 'F:\Data\F-G1S\F024-live\Data\';      % Output location of mat files
settings.IF_imagesessions = {'N:\F024-live-fixed\'};    
                                                            % Cell array containing all imaging sessions to be matched
settings.bgcmospath = '';
settings.crop_save = '';           % Output cropped images, leave empty if no save
settings.mask_save = 'N:\F024-live-fixed\Raw';           % Output mask, leave empty if no save
settings.maskname = 'mask';
settings.primaryMaskRound = 1;
settings.maskIndex = [1];

%%% General parameters
settings.microscope = 'nikon';
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
% settings.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'; %'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'
settings.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.magnification = 20;     % 10 = 10x or 20 = 20x
settings.binsize = 1;            % 1 = bin 1 or 2 = bin 2
settings.postbin = [.5 ];           % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
% settings.match10x20x = 1;
% settings.numsites = 9;           % Number of sites imaged for 20x to 10x (double check order is correct)
% settings.scaleLive = 1;          % Scale live to have same pixel size as final IF image
settings.signals = {{'DAPI','YFP','FarRed'}};  % Each imaging session in a cell inside cell array
% settings.nucLive = 'CFP_';       % Label for IX micro live nucleus (assume first tiff in stack for nikon)
% settings.liveChannel = 4;         % If nikon use this tiff in stack for nuclear live
% settings.manualLiveFrame = [];   % Manually set last frame for live imaging, otherwise use tracedata if empty

%%% Segmentation parameters
settings.segmethod = {'concavity'};  % 'log' or 'single' or 'double'
settings.nucr = 12;                           % 10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100;                    % 100
settings.boulderarea = 1500;                  % 1500 
settings.blobthreshold = -0.03;               % For blobdector 'log' segmentation
settings.blurradius = 3;                      % Blur for blobdetector
settings.soliditythresh = 0.80;               % Minimum solidity of nucleus allowed
settings.badFrameCheck = .25;                 % Min fraction of cells untracked before logging error, 0 for no check
settings.split_mult = 1;

%%% Quantification parameters
settings.bias = {[1 1 1] };
settings.biasall = {[1 1 1 ]};
settings.sigblur = {[0 0 0 ] };
settings.signal_foreground = { [0 0 0 ] };
settings.bleedthrough = {[0 0 0 ] };   % Calculate bleedthrough for channel
settings.bleedthroughslope = {};              % Cell array of cell arrays
settings.bleedthroughoff = {};                % Cell array of cell arrays
settings.ringcalc = {[0 1 1 ] };
settings.ringthresh = { [0 0 0 ] };    % Threshold for foreground ring pixels
settings.punctacalc = { [0 0 0 ]};     % Analyze puncta in nucleus
settings.punctaThresh = {{[],[],[]},{}};      % Threshold for tophat filter for puncta
settings.punctatopsize = 2;                   % Top hat filter size
settings.cytopuncta = { [0 0 0 ]};
settings.thickenradius = 2*settings.nucr; 

settings.localbg = { [0 0 0 0 0]};        % Measure local background in box around cell
settings.minringsize = 100;
settings.bgcmoscorrection = 1;                % 1 = correct for shading; 0 = dont correct for shading;
settings.compression = 4;
settings.bgsubmethod = {{'global nuclear','none', 'global nuclear'}}; % 'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25 25 ]};  % Set bg as perctile for each channel
settings.frameIF=1;
 
Fixed_IF( settings, 2,5,1,1)
%% Run timelapse
rows = [2:7];
cols = [2:11];
sites = [1:24];
Parallelwell({@Fixed_IF},settings,rows,cols,sites)


%% Rename IF 
nikon_rename('N:\F024-live-fixed-MCM\20211018_120457_834',...
     'N:\F024-live-fixed-MCM\Raw','COPY',false,'SINGLE_SITE',false)


%% Setup IF images
clear settings
settings.data_path = 'F:\Data\F-G1S\F024-live\Data-MCM\';      % Output location of mat files
settings.IF_imagesessions = {'N:\F024-live-fixed-MCM\'};    
                                                            % Cell array containing all imaging sessions to be matched
settings.bgcmospath = '';
settings.crop_save = '';           % Output cropped images, leave empty if no save
settings.mask_save = 'N:\F024-live-fixed-MCM\Raw';           % Output mask, leave empty if no save
settings.maskname = 'mask';
settings.primaryMaskRound = 1;
settings.maskIndex = [1];

%%% General parameters
settings.microscope = 'nikon';
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
% settings.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'; %'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'
settings.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.magnification = 20;     % 10 = 10x or 20 = 20x
settings.binsize = 1;            % 1 = bin 1 or 2 = bin 2
settings.postbin = [.5 ];           % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
% settings.match10x20x = 1;
% settings.numsites = 9;           % Number of sites imaged for 20x to 10x (double check order is correct)
% settings.scaleLive = 1;          % Scale live to have same pixel size as final IF image
settings.signals = {{'DAPI','YFP','FarRed'}};  % Each imaging session in a cell inside cell array
% settings.nucLive = 'CFP_';       % Label for IX micro live nucleus (assume first tiff in stack for nikon)
% settings.liveChannel = 4;         % If nikon use this tiff in stack for nuclear live
% settings.manualLiveFrame = [];   % Manually set last frame for live imaging, otherwise use tracedata if empty

%%% Segmentation parameters
settings.segmethod = {'concavity'};  % 'log' or 'single' or 'double'
settings.nucr = 12;                           % 10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100;                    % 100
settings.boulderarea = 1500;                  % 1500 
settings.blobthreshold = -0.03;               % For blobdector 'log' segmentation
settings.blurradius = 3;                      % Blur for blobdetector
settings.soliditythresh = 0.80;               % Minimum solidity of nucleus allowed
settings.badFrameCheck = .25;                 % Min fraction of cells untracked before logging error, 0 for no check
settings.split_mult = 1;

%%% Quantification parameters
settings.bias = {[1 1 1] };
settings.biasall = {[1 1 1 ]};
settings.sigblur = {[0 0 0 ] };
settings.signal_foreground = { [0 0 0 ] };
settings.bleedthrough = {[0 0 0 ] };   % Calculate bleedthrough for channel
settings.bleedthroughslope = {};              % Cell array of cell arrays
settings.bleedthroughoff = {};                % Cell array of cell arrays
settings.ringcalc = {[0 1 1 ] };
settings.ringthresh = { [0 0 0 ] };    % Threshold for foreground ring pixels
settings.punctacalc = { [0 0 0 ]};     % Analyze puncta in nucleus
settings.punctaThresh = {{[],[],[]},{}};      % Threshold for tophat filter for puncta
settings.punctatopsize = 2;                   % Top hat filter size
settings.cytopuncta = { [0 0 0 ]};
settings.thickenradius = 2*settings.nucr; 

settings.localbg = { [0 0 0 0 0]};        % Measure local background in box around cell
settings.minringsize = 100;
settings.bgcmoscorrection = 1;                % 1 = correct for shading; 0 = dont correct for shading;
settings.compression = 4;
settings.bgsubmethod = {{'global nuclear','none', 'global nuclear'}}; % 'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25 25 ]};  % Set bg as perctile for each channel
settings.frameIF=1;
 
Fixed_IF( settings, 5,5,1,0)
%% Run timelapse
rows = [5:7];
cols = [2:11];
sites = [1:24];
Parallelwell({@Fixed_IF},settings,rows,cols,sites)


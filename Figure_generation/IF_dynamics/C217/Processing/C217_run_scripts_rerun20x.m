%% Prep match IF to live processing
%%% Paths
clear s;
s.data_path = 'F:\Data\C-Cdt1\C217-live\Data_rerun\';      % Output location of mat files
s.IF_imagesessions = {'H:\8TB6\Data\IXM\C217-IF1\',...
    'H:\8TB6\Data\IXM\C217-IF2\',...
    'H:\8TB6\Data\IXM\C217-IF3\',...
    'H:\8TB6\Data\IXM\C217-IF4\',...
    'H:\8TB6\Data\IXM\C217-IF5\'};    
                                                            % Cell array containing all imaging sessions to be matched
s.live_imagepath = 'H:\8TB6\Data\IXM\C217-live\Raw\';    % Live imaging raw image path
s.bgcmospath = 'C:\Users\Meyerlab\Documents\GitHub\cell-cycle-tracking\BGimages\IX_Liu\cmosoffset_bin1.mat';
s.crop_save = 'H:\8TB6\Data\IXM\C217-IFcrop-rerun\';           % Output cropped images, leave empty if no save
s.mask_save = 'H:\8TB6\Data\IXM\C217-IFcrop-rerun\';           % Output mask, leave empty if no save

%%% General parameters
s.microscope = 'IXM';
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
s.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
s.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
s.magnification = 20;     % 10 = 10x or 20 = 20x
s.binsize = 1;            % 1 = bin 1 or 2 = bin 2
s.postbin = 0;           % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
s.match10x20x = 1;
s.numsites = 4;         % Number of sites imaged for 20x to 10x (double check order is correct)
s.maxjit = 500;
s.scaleLive = 2;          % Scale live to have same pixel size as final IF image
s.signals = {{'DAPI_','YFP_','FarRed_'}, {'DAPI_','YFP_'}, {'DAPI_','YFP_'}, {'DAPI_','YFP_'}, {'DAPI_','FITC_'}};  % Each imaging session in a cell inside cell array
s.nucLive = 'CFP_';       % Label for IX micro live nucleus (assume first tiff in stack for nikon)
s.manualLiveFrame = [];   % Manually set last frame for live imaging, otherwise use tracedata if empty

%%% Quantification parameters
s.bias = {[1 1 1], [1 1], [1 1], [1 1], [1 1]};
s.biasall = {[0 0 0], [0 0], [0 0], [0 0], [0 0]};
s.sigblur = {[0 0 0], [0 0], [0 0], [0 0], [0 0]};
s.signal_foreground = {[0 0 0], [0 0], [0 0], [0 0], [0 0]};
s.bleedthrough = {[0 0 0], [0 0], [0 0], [0 0], [0 0]};   % Calculate bleedthrough for channel
s.bleedthroughslope = {};              % Cell array of cell arrays
s.bleedthroughoff = {};                % Cell array of cell arrays
s.ringcalc = {[1 1 1], [0 1], [0 1], [0 1], [0 1]};
s.ringthresh = {[0 0 0], [0 0], [0 0], [0 0], [0 0]};    % Threshold for foreground ring pixels
s.punctacalc = {[0 0 0], [0 0], [0 0], [0 0], [0 0]};     % Analyze puncta in nucleus
s.punctaThresh = {{[],[],[]},{}};      % Threshold for tophat filter for puncta
s.punctatopsize = 2;                   % Top hat filter size
s.localbg = {[0 0 0], [0 0], [0 0], [0 0], [0 0]};       % Measure local background in box around cell
s.minringsize = 100;
s.bgcmoscorrection = 1;                % 1 = correct for shading; 0 = dont correct for shading;
s.compression = 4;
s.bgsubmethod = {{'global nuclear','global nuclear', 'global nuclear'}, ...
    {'global nuclear', 'global nuclear'},...
    {'global nuclear', 'global nuclear'},...
    {'global nuclear', 'global nuclear'},...
    {'global nuclear', 'global nuclear'}}; % 'global nuclear','global cyto','tophat','semi-local nuclear'
s.bgperctile = {[25 25 25], [25 25], [25 25], [25 25], [25 25]};  % Set bg as perctile for each channel
s.frameIF=1;
 
%%% Segmentation parameters
s.segmethod = {'concavity','thresh','thresh','thresh','thresh'};  % 'log' or 'single' or 'double'
s.nucr = 24;                           % 10x bin1 = 12 and 20x bin2 = 12
s.debrisarea = 400;                    % 100
s.boulderarea = 6000;                  % 1500 
s.blobthreshold = -0.03;               % For blobdector 'log' segmentation
s.blurradius = 3;                      % Blur for blobdetector
s.soliditythresh = 0.80;               % Minimum solidity of nucleus allowed
s.badFrameCheck = .25;                 % Min fraction of cells untracked before logging error, 0 for no check

%%% Tracking parameters
s.distthresh = 3*s.nucr;
s.arealowthresh = -.75;
s.areahighthresh = 1;

%  Timelapse_addIF(s,2,7,2,0)


%% Run timelapse
rows = [2:3];
cols = [2:11];
sites = [1:4];
 Parallelwell({@Timelapse_addIF},s,rows,cols,sites)
 
 
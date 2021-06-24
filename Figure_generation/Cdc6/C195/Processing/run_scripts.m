% IX_Micro_rename('H:\Nalin\201901001-C195-IF1\2019-10-01\10686','H:\8TB4\Data\C195-total-IF1\Raw',{'DAPI','CFP','YFP','RFP'});
% IX_Micro_rename('H:\Nalin\201901001-C195-IF2\2019-10-02\10689','H:\8TB4\Data\C195-total-IF2\Raw',{'DAPI','FarRed'});

%% Calculate Bias IF round 1
% clear s
% s.biaspath='H:\8TB4\Data\C195-total-IF1\Bias';
% s.imagepath='H:\8TB4\Data\C195-total-IF1\Raw';
% s.shadingpath='F:\MATLAB\BGimages\IX_Liu\cmosoffset_bin1.mat';
% s.nucradius=24;%12;
% s.names={
%     'DAPI';
%     'CFP';
%     'YFP';
%     'RFP';
%     };
% s.formatCode = '';
% s.microscope = 'IXM';
% s.row_mat = [2:5];
% s.col_mat = [1];
% s.site_mat = [1:37];
% s.frame_mat=[1];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 0; %1 to mask nuclei
% s.adaptive = [0 0 0];
% s.dilate = {.75,1.25,.5}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 15;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0)

%% Calculate Bias IF round 2
% clear s
% s.biaspath='H:\8TB4\Data\C195-total-IF2\Bias';
% s.imagepath='H:\8TB4\Data\C195-total-IF2\Raw';
% s.shadingpath='F:\MATLAB\BGimages\IX_Liu\cmosoffset_bin1.mat';
% s.nucradius=24;%12;
% s.names={
%     'DAPI';
%     'FarRed';
%     };
% s.formatCode = '';
% s.microscope = 'IXM';
% s.row_mat = [2:5];
% s.col_mat = [1];
% s.site_mat = [1:37];
% s.frame_mat=[1];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 0; %1 to mask nuclei
% s.adaptive = [0 0 0];
% s.dilate = {.75,1.25,.5}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 15;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0)


%% Setup Process IF
%%% Paths
settings.data_path = 'F:\Data\C-Cdt1\C195-IF-total\Data\';
settings.IF_imagesessions = {'H:\8TB4\Data\C195-total-IF1\','H:\8TB4\Data\C195-total-IF2\'};
settings.formatCode = '';
settings.bgcmospath = 'F:\MATLAB\BGimages\IX_Liu\cmosoffset_bin1.mat';
settings.maskpath = 'H:\8TB4\Data\C195-total-IFcrop\';
settings.crop_save = 'H:\8TB4\Data\C195-total-IFcrop\';
settings.microscope = 'IXM';

%%% File Options
settings.maskname = 'nucedge_';

%%% General parameters
settings.magnification = 20; %10 = 10x or 20 = 20x
settings.binsize = 1; %1 = bin 1 or 2 = bin 2
settings.postbin = .5;
settings.signals = {{'DAPI_','CFP_','YFP_','RFP_'}, {'DAPI_','FarRed_'}};

%%% Quantification parameters
settings.bias = {[1 1 1 1],[1 1]};
settings.biasall = {[0 0 0 0], [0 0]};
settings.sigblur = {[0 0 0 0], [ 0 0]};
settings.signal_foreground = { [0 0 0 0 ], [ 0 0]};
settings.ringcalc = {[0 0 0 0],[ 0 1]};
settings.ringthresh = {[0 0 0 0], [0 0]};
settings.punctacalc = {[0 0 0 0], [ 0 0]};
settings.punctaThresh = {{[],[]},{[],[],[]}};
settings.punctatopsize = 2;
settings.localbg = {[0 0 0 0],[0 0]};
settings.minringsize = 100;
settings.bleedthrough = {[0 0 0 0], [ 0 0]};
settings.bleedthroughpth = {};
settings.bgcmoscorrection = 1; %1 = correct for shading; 0 = dont correct for shading;
settings.bgsubmethod = {{'global nuclear', 'global nuclear', 'global nuclear', 'global nuclear',},{'global nuclear', 'global nuclear'}};
%'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25 25 25],[25 25]};

%%% Segmentation parameters
settings.segmethod = {'concavity', 'thresh'}; %'log' or 'single' or 'double'
settings.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100; %100
settings.boulderarea = 1500; %1500 % less restrictive for C152 for >4N cells
settings.blobthreshold = -0.03;
settings.blurradius = 3;
settings.soliditythresh = 0.80;
settings.compression = 4;

%Fixed_IF(settings,2,2,1,0)

%% Process
rows = [1:8];
cols = [2:12];
sites = [1:37];
Parallelwell({@Fixed_IF},settings,rows,cols,sites,0);

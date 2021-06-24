function IF_addToTimelapse10x20x_E1099( row,col,site )
% row=3;
% col=3;
% site=1;

%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Paths
settings.data_path = 'F:\Data\C-Cdt1\E1099-live\';
settings.IF_imagepath = 'G:\Data\E1099-IF-2\Raw\';
settings.live_imagepath = 'G:\Data\E1099-live\Raw\';
settings.bgcmospath = 'C:\Users\Meyerlab\OneDrive - Leland Stanford Junior University\Meyer\MATLAB\BGimages\cmosoffset_bin2.mat';
settings.biaspath = 'G:\Data\E1099-IF-2\Bias_all\';
settings.maskpath = 'G:\Data\E1099-IF-2\Mask\';

%%% File Options
settings.maskwrite_option = 0; %1 = save an image with the nucmask; 0 = dont save an image with the nucmask

%%% Experiment parameters
settings.magnification = 20; %10 = 10x or 20 = 20x
settings.binsize = 2; %1 = bin 1 or 2 = bin 2
settings.signals = {'DAPI_','RFP_','FarRed_'};
settings.nucLive = 'CFP_';
settings.postbin = 0;

%%% Options
settings.ringcalc = [0 0 0 0];
settings.localbg = [1 1 1 1];
settings.ringthresh = [100 100 100 100];
settings.bias = [1 1 1 1];
settings.bleedthrough = [0 0 0 0];
settings.bleedthroughpth = {'' '' '' ''};
settings.bgcmoscorrection = 1; %1 = correct for shading; 0 = dont correct for shading;
settings.firstsegmethod = 'double marker'; %'log' or 'single' or 'double'
settings.bgsubmethod = {'global nuclear','global nuclear', 'global nuclear','global nuclear'}; %'global nuclear','global cyto','tophat','semi-local nuclear'
settings.frameIF=1;

%%% Segmentation parameters
settings.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100;
settings.boulderarea = 1500;
settings.blobthreshold = -0.03;
settings.blurradius = 3;
settings.remove_smears_option = 0;
settings.soliditythresh = 0.85;
settings.compression = 4;

%% Process Image
shot10x=[num2str(row),'_',num2str(col),'_',num2str(site)];
if exist([settings.data_path,'\Data\','tracedata_',shot10x,'.mat']) & ~exist([settings.data_path,'\Data\','IF_',shot10x,'.mat'])
    
    timetotal=tic;
    
    numsites=4;
    [upperleft,upperright,lowerleft,lowerright]=siteChange(numsites,site);
    
    IFdata_upperleft = processAllWells(row, col, site, upperleft, settings);
    IFdata_upperright = processAllWells(row, col, site, upperright, settings);
    IFdata_lowerleft = processAllWells(row, col, site, lowerleft, settings);
    IFdata_lowerright = processAllWells(row, col, site, lowerright, settings);
    
    IFdata = ones(size(IFdata_upperleft))*NaN;
    IFdata = addGoodRows(IFdata,IFdata_upperleft);
    IFdata = addGoodRows(IFdata,IFdata_upperright);
    IFdata = addGoodRows(IFdata,IFdata_lowerleft);
    IFdata = addGoodRows(IFdata,IFdata_lowerright);
    
    save([settings.data_path,'\Data\','IF_',shot10x,'.mat'],'IFdata');
    header = output_names_IF(settings.signals,settings.localbg, settings.ringcalc);
    save([settings.data_path,'\Data\settings_IF.mat'],'settings','header'); %saves all the settings to the Data folder. will overwrite with the most recent
    
    toc(timetotal);
end
end


%%
function [IFdata] = processAllWells(row, col, site, subsite, settings)

shot10x=[num2str(row),'_',num2str(col),'_',num2str(site)];
shot20x=[num2str(row),'_',num2str(col),'_',num2str(subsite)];
dataDir=[settings.data_path,'\Data\'];


raw20xDir = [settings.IF_imagepath,shot20x,'\',shot20x,'_'];
rawLiveDir = [settings.live_imagepath,shot10x,'\',shot10x,'_'];
maskDir = [settings.maskpath,shot20x];

if ~exist(maskDir,'dir') && settings.maskwrite_option
    mkdir(maskDir);
end


%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names = settings.signals;
nucname = names{1};

%%% Segmentation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr = settings.nucr;
debrisarea = settings.debrisarea;
boulderarea = settings.boulderarea;
blobthreshold = settings.blobthreshold;
solidity = settings.soliditythresh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% load timelapse data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([dataDir,'tracedata_',shot10x,'.mat'],'tracedata','genealogy','jitters');
[totalcells,totalframes,totalsignals]=size(tracedata);
%%% get previous mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawprev=single(imread([rawLiveDir,settings.nucLive,num2str(totalframes),'.tif']));
[prevheight,prevwidth]=size(rawprev);
[nuc_mask_prev,~]=blobdetector_foreground_2(log(rawprev),nucr,blobthreshold,debrisarea);
%%% bgcmos correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([settings.bgcmospath],'cmosoffset');
bgcmos = cmosoffset;
if any(settings.ringcalc)
    if (settings.magnification == 10 & settings.binsize == 1) | ...
            (settings.magnification == 20 & (settings.binsize == 2 | settings.postbin))
        innerrad = 1; outerrad = 5;
    else
        innerrad = 2; outerrad = 10;
    end
end

for i = 1:length(names)
    if settings.bias(i)
        load([settings.biaspath,names{i},num2str(site),'.mat']); bias_cell{i} = bias;
    end
end

%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(names)
    raw{j} = single(imread([raw20xDir,names{j},num2str(settings.frameIF),'.tif']));
    if settings.bgcmoscorrection
        raw{j} = (raw{j}-bgcmos);
        raw{j}(raw{j}<0) = 0;
    end
    
    if settings.bias(j)
        raw{j} = raw{j}./bias_cell{j};
    end
    
    if settings.postbin
        raw{j} = imresize(raw{j},.5);
    end
    normraw{j}=(raw{j}-min(raw{j}(:)))./(max(raw{j}(:))-min(raw{j}(:)));
    
    if settings.localbg(j)
        foreground{j} = threshmask(raw{j},settings.blurradius);
        %foreground{j}=nuc_mask;
    end
end

%NEfile=[maskdir,nucedgename,'stain.tif'];


%%% remove smears %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if settings.remove_smears_option
%         foregroundthresh = 1000;
%         areathresh = 20000;
%         bgprctile = 20;
%         raw1 = removesmears_1(raw1,foregroundthresh,areathresh,bgprctile);
%         raw2 = removesmears_1(raw2,foregroundthresh,areathresh,bgprctile);
%         if settings.signal3_option;raw3 = removesmears_1(raw3,foregroundthresh,areathresh,bgprctile);end
%         if settings.signal4_option;raw4 = removesmears_1(raw4,foregroundthresh,areathresh,bgprctile);end
%         if settings.signal5_option;raw5 = removesmears_1(raw5,foregroundthresh,areathresh,bgprctile);end
%     end
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
segmethod = settings.firstsegmethod;
blurradius = settings.blurradius;
switch segmethod
    case 'log'
        nuc_mask = blobdetector_4(log(raw{1}),nucr,blobthreshold,debrisarea);
    case 'single'
        nuc_mask = threshmask(raw{1},blurradius);
        nuc_mask = markershed(nuc_mask,round(nucr*2/3));
    case 'double'
        %nuc_mask = threshmask(raw{1},blurradius);
        nuc_mask = imbinarize(raw{1},550);
        nuc_mask = imclose(nuc_mask,strel('disk',3));
        nuc_mask = imfill(nuc_mask,'holes');
        %imagesc(nuc_mask);
        nuc_mask = markershed(nuc_mask,round(nucr*2/3));
        nuc_mask = secondthresh(raw{1},blurradius,nuc_mask,boulderarea*2);
    case 'double marker'
        nuc_mask = threshmask(raw{1},blurradius);
        nuc_mask = secondthresh_all(raw{1},blurradius,nuc_mask);
        nuc_mask = markershed(nuc_mask,round(nucr*.4));
    case 'multithresh'
        nuc_mask = threshmask_multi(raw{1},blurradius,3, 1e-04);
        nuc_mask = markershed_filter(nuc_mask,round(nucr*.4),4);
end
whole_mask = nuc_mask;
nuc_mask = bwareaopen(nuc_mask,debrisarea);

%%% Deflection-Bridging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask = segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
nuc_mask = excludelargeandwarped_3(nuc_mask,boulderarea,settings.soliditythresh);


%{
%%% debugging: view images %%%%%%%%%%
extractmask = bwmorph(nuc_mask,'remove');
tempframe = imadjust(mat2gray(raw{1}));
tempframe(:,:,2) = extractmask;
tempframe(:,:,3) = 0;
figure,imshow(tempframe);

nuc_info = struct2cell(regionprops(nuc_mask,'Area')');
nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
hist(nuc_area,100);

%}
%%% calculate background: semi-local %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compression = settings.compression;

for j = 1:length(names)
    if settings.localbg(j)
        mask_total = whole_mask | foreground{j};
        nanmask{j} = imdilate(mask_total,strel('disk',nucr/2));
        nanmaskcyto{j} = imdilate(mask_total,strel('disk',nucr*2));
    else
        nanmask{j} = imdilate(whole_mask,strel('disk',nucr/2));
        nanmaskcyto{j} = imdilate(whole_mask,strel('disk',nucr*2));
    end
    blur = imfilter(raw{j},fspecial('disk',3),'symmetric');
    %blur = raw{j};
    switch settings.bgsubmethod{j}
        case 'global nuclear'
            real{j} = bgsubmasked_global_NR(blur,nanmask{j},1,compression,25);
        case 'global cyto'
            real{j} = bgsubmasked_global_2(blur,nanmaskcyto{j},1,compression,25);
        case 'tophat'
            real{j} = imtophat(blur,strel('disk',nucr,0));
        case 'semi-local nuclear'
            real{j} = bgsubmasked_global_2(blur,nanmask{j},11,compression,25);
        case 'none'
            real{j} = blur;
    end
end

%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(names)
    if settings.bleedthrough(j)
        load([dataDir,settings.bleedthroughpth{j}],'bleedthroughrate');
        real1bleedthrough = real{j-1}*bleedthroughrate(2)+bleedthroughrate(1);
        real{j} = real{j}-real1bleedthrough;
    end
end

%%% pad raw images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cropheight,cropwidth] = size(raw{1});
tophalf = floor((prevheight-cropheight)/2);
if mod(cropheight,2) == 1
    bottomhalf = tophalf+1;
else
    bottomhalf = tophalf;
end
lefthalf = floor((prevwidth-cropwidth)/2);
if mod(cropwidth,2) == 1
    righthalf = lefthalf + 1;
else
    righthalf = lefthalf;
end

whole_mask_padded = padarray(whole_mask,[tophalf lefthalf],'pre');
whole_mask_padded = padarray( whole_mask_padded,[bottomhalf righthalf],'post');
nuc_mask_padded = padarray(nuc_mask,[tophalf lefthalf],'pre');
nuc_mask_padded = padarray(nuc_mask_padded,[bottomhalf righthalf],'post');

for j = 1:length(names)
    real_padded{j} = padarray(real{j},[tophalf lefthalf],'pre');
    real_padded{j} = padarray( real_padded{j},[bottomhalf righthalf],'post');
    nanmask_padded{j} = padarray(nanmask{j},[tophalf lefthalf],'pre');
    nanmask_padded{j} = padarray( nanmask{j},[bottomhalf righthalf],'post');
end

%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_label = bwlabel(nuc_mask_padded);
nuc_info = struct2cell(regionprops(nuc_mask_padded,real_padded{1},'Area','Centroid','MeanIntensity')');
nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center = squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density = squeeze(cell2mat(nuc_info(3,1,:)));

%%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mass = nuc_density.*nuc_area;
%%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[reljitx, reljity] = registerimages(nuc_mask_prev,nuc_mask_padded);
[padheight, padwidth] = size(nuc_mask_padded);
jitcoors = getcropcoors([padheight,padwidth],reljitx,reljity);
prev_mask_jit = nuc_mask_prev(jitcoors(1,1):jitcoors(1,2),jitcoors(1,3):jitcoors(1,4));
nuc_mask_jit = nuc_mask_padded(jitcoors(2,1):jitcoors(2,2),jitcoors(2,3):jitcoors(2,4));
for j = 1:length(names)
    real_padded_jit{j} = real_padded{j}(jitcoors(2,1):jitcoors(2,2),jitcoors(2,3):jitcoors(2,4));
end
reljitter = [reljitx, reljity];
prevjitter = jitters(totalframes,:);
IFjitter = prevjitter + reljitter;
nuc_center(:,1) = nuc_center(:,1) + IFjitter(1);
nuc_center(:,2) = nuc_center(:,2) + IFjitter(2);
%%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initdata = [nuc_center(:,1), nuc_center(:,2), nuc_area, nuc_mass];
%%% track & correct merges (update centers, masses and labels) %%%%%%%%%%%%
debugpackage = {rawprev, nuc_mask_prev, nuc_mask_padded, prevjitter, reljitter};
distthresh = 2*nucr; arealowthresh = -.4; areahighthresh = .2;
[tracked,nuc_label_track] = adaptivetrack_IF_1(tracedata(:,totalframes,1:4), initdata, nuc_label, nucr,...
    distthresh, arealowthresh, areahighthresh, debugpackage);
approx_perc_lost_cells = (sum(~isnan(tracedata(:,totalframes,1)))/4 - sum(~isnan(tracked(:,1))))/...
    sum(~isnan(tracedata(:,totalframes,1)));
tracked_cells = length(unique(nuc_label));
nuc_label_track_jit = nuc_label_track(jitcoors(2,1):jitcoors(2,2),jitcoors(2,3):jitcoors(2,4));
%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.maskwrite_option
    extractmask = bwmorph(nuc_label_track,'remove');
    imwrite(uint16(extractmask),NEfile);
end
%%% account for all objects, including untracked %%%%%%%%%%%%%%%%%%%%%%%%%%
numtracked = numel(tracked(:,1));
cellid = find(~isnan(tracked(:,1)));
numlivecells = numel(cellid);
trackedobjects = nuc_label_track>0;
untrackedobjects = nuc_mask_padded-trackedobjects;
untracked_label = bwlabel(untrackedobjects);
untracked_label = untracked_label+numtracked;
untracked_label(untracked_label == numtracked) = 0;
allobjects_label = nuc_label + untracked_label;

%%% Measure features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells = numel(tracked(:,1));
nuc_center = tracked(:,[1 2]);
nuc_area = tracked(:,3);
nuc_mass = tracked(:,4);
cellid = find(~isnan(tracked(:,1)));
numlivecells = numel(cellid);
nuc_info = regionprops(nuc_label_track, 'PixelIdxList', 'Area', 'Centroid');

%Initialize storage variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nanvec = ones(numcells,1)*NaN;
sigringmedian{1} = [];
localbg{1} = [];
IF_position  = [nanvec nanvec];
for j = 1:length(names)
    sigmean{j} = nanvec;
    sigmedian{j} = nanvec;
    %sigring75th{j} = nanvec;
    if settings.ringcalc(j)
        sigringmedian{j} = nanvec;
    end
    if settings.localbg(j)
        localbg{j} = nanvec;
    end
end

%Measure features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(names)
    real_temp = real_padded{j};
    real_temp_masked = real_temp;
    real_temp_masked(nanmask_padded{j}) = NaN;
    
    for n = 1:numlivecells
        cc=cellid(n);
        sigmean{j}(cc) = mean(real_temp(nuc_info(cc).PixelIdxList));
        sigmedian{j}(cc) = median(real_temp(nuc_info(cc).PixelIdxList));
        IF_position(cc,:) = nuc_info(cc).Centroid - [tophalf lefthalf];
        if settings.localbg(j)
            localbg{j}(cc) = getlocalbg(real_temp_masked, nuc_info(cc).Centroid,100, 25, .15);
        end
        if settings.ringcalc(j)
            if cc>numel(ring_info)
                break;
            end
            ringall = real_temp(ring_info(cc).PixelIdxList);
            %ringall(ringall>prctile(ringall,98)) = [];
            %sigring75th{j}(n) = prctile(ringall,75);
            ringforeground = ringall(ringall>settings.ringthresh(j));
            if numel(ringforeground)<100
                ringforeground = ringall;
            end
            if numel(ringall)>100
                sigringmedian{j}(cc) = nanmedian(ringforeground);
            end
        end
    end
end
IFdata = [IF_position(:,1),IF_position(:,2),nuc_area,...
    cell2mat(sigmean),cell2mat(sigmedian),cell2mat(localbg),...
    cell2mat(sigringmedian)];
end

%{
%%% debugging: view images %%%%%%%%%%
%Check jitter
extractmask=bwmorph(nuc_mask_jit,'remove');
[tempheight,tempwidth]=size(extractmask);
tempframe=zeros(tempheight,tempwidth,3);
tempframe(:,:,1)=prev_mask_jit;
tempframe(:,:,2)=nuc_mask_jit;
figure,imshow(tempframe);
hold on

%Check tracking
extractmask1=bwmorph(nuc_mask_jit,'remove');
extractmask2=bwmorph(nuc_label_track_jit,'remove');
[tempheight,tempwidth]=size(extractmask1);
tempframe=zeros(tempheight,tempwidth,3);
tempframe(:,:,1)=extractmask1;
tempframe(:,:,2)=extractmask2;
tempframe(:,:,3)=imadjust(mat2gray(real_padded_jit{1}));
figure,imshow(tempframe);
hold on,
scatter(tracedata(:,totalframes,1) - jitters(end,1), tracedata(:,totalframes,2)- jitters(end,2))

%Check IF positions
figure, imshow(real{1},[]), hold on
scatter(IF_position(:,1),IF_position(:,2))

%%% view images %%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_label,'remove');
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray(raw{1}));
tempframe(:,:,2)=extractmask;
%tempframe(:,:,3)=marker_mask;
figure,imshow(tempframe);



%}
function [upperleft,upperright,lowerleft,lowerright]=siteChange(numsites,site)
if numsites==4
    switch site
        case 1
            upperleft=1;
            upperright=2;
            lowerleft=5;
            lowerright=6;
        case 2
            upperleft=3;
            upperright=4;
            lowerleft=7;
            lowerright=8;
        case 3
            upperleft=9;
            upperright=10;
            lowerleft=13;
            lowerright=14;
        case 4
            upperleft=11;
            upperright=12;
            lowerleft=15;
            lowerright=16;
    end
elseif numsites==6
    switch site
        case 1
            upperleft=1;
            upperright=2;
            lowerleft=7;
            lowerright=8;
        case 2
            upperleft=3;
            upperright=4;
            lowerleft=9;
            lowerright=10;
        case 3
            upperleft=5;
            upperright=6;
            lowerleft=11;
            lowerright=12;
        case 4
            upperleft=13;
            upperright=14;
            lowerleft=19;
            lowerright=20;
        case 5
            upperleft=15;
            upperright=16;
            lowerleft=21;
            lowerright=22;
        case 6
            upperleft=17;
            upperright=18;
            lowerleft=23;
            lowerright=24;
    end
end
end

function updateddata=addGoodRows(orgdata,newdata)
goodrows=~isnan(newdata(:,1));
updateddata=orgdata;
updateddata(goodrows,:)=newdata(goodrows,:);
end
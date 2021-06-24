function IF_AddToTimelapse_C103(row,col,site)
row=2;
col=2;
site=1;
debug_mode = 0;


%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Paths
settings.data_path = 'F:\Data\C103-live\';
settings.IF_imagepath = 'D:\Data\C103-IF\Raw\';
settings.live_imagepath = 'H:\Data\C013-live\Raw\';
settings.bgcmospath = 'C:\Users\Meyerlab\OneDrive - Leland Stanford Junior University\Meyer\MATLAB\BGimages\cmosoffset_bin1.mat';
settings.biaspath = 'D:\Data\C013-IF\Bias_all\';
settings.maskpath = 'D:\Data\C013-IF\Mask\';

%%% File Options
settings.separatedirectories_option = 1; %1 = Folder for each well,col,site; 0 = all files in 1 folder
settings.maskwrite_option = 0; %1 = save an image with the nucmask; 0 = dont save an image with the nucmask

%%% Experiment parameters
settings.magnification = 10; %10 = 10x or 20 = 20x
settings.binsize = 1; %1 = bin 1 or 2 = bin 2
settings.signals = {'DAPI_','Cy3_','FarRed_'};
settings.postbin = 0;

%%% Options
settings.ringcalc = [0 0 0];
settings.localbg = [1 1 1];
settings.ringthresh = [100 100 100];
settings.bias = [1 1 1];
settings.bleedthrough = [0 0 0];
settings.bleedthroughpth = {'' '' '' ''};
settings.bgcmoscorrection = 1; %1 = correct for shading; 0 = dont correct for shading;
settings.firstsegmethod = 'multithresh'; %'log' or 'single' or 'double'
settings.bgsubmethod = {'global nuclear','global nuclear', 'global nuclear'}; %'global nuclear','global cyto','tophat','semi-local nuclear'
settings.Frame_of_IF=1;

%%% Segmentation parameters
settings.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100;
settings.boulderarea = 1500;
settings.blobthreshold = -0.03;
settings.blurradius = 3;
settings.remove_smears_option = 0;
settings.soliditythresh = 0.7;
settings.compression = 4;



%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF_imagepath=settings.IF_imagepath;
Frame_of_IF=settings.Frame_of_IF;
biasdir=settings.biaspath;
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
bgcmospath = settings.bgcmospath;
datadir=[settings.data_path,'\Data\'];

if exist([datadir,'tracedata_',shot,'.mat']) & ~exist([datadir,'IF_',shot,'.mat'])
    separatedirectories=1;
    if separatedirectories == 1
        rawdir = [settings.IF_imagepath,shot,'\',shot,'_'];
        rawdir2 = [settings.live_imagepath,shot,'\',shot,'_'];
        maskdir = [settings.maskpath,shot];
    else
        rawdir = [settings.IF_imagepath,shot,'_'];
        maskdir = [settings.maskpath,'\',shot,'_'];
    end
    if ~exist(maskdir,'dir') && settings.maskwrite_option
        mkdir(maskdir);
    end
    maskwrite=0;
    
    
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
    timetotal=tic;
    %%% load timelapse data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
    [totalcells,totalframes,totalsignals]=size(tracedata);
    %%% get previous mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rawprev=single(imread([rawdir2,'CFP_',num2str(totalframes),'.tif']));
    [height,width]=size(rawprev);
    [nuc_mask_prev,~]=blobdetector_foreground_2(log(rawprev),nucr,blobthreshold,debrisarea);
    %%% bgcmos correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([bgcmospath],'cmosoffset');
    bgcmos = cmosoffset;
    if any(settings.ringcalc)
        if (settings.magnification == 10 & settings.binsize == 1) | (settings.magnification == 20 & settings.binsize == 2)
            innerrad = 1; outerrad = 5;
        else
            innerrad = 2; outerrad = 10;
        end
    end
    
    for i = 1:length(names)
        if settings.bias(i)
            load([biasdir,names{i},num2str(site),'.mat']); bias_cell{i} = bias;
        end
    end
    
    %%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:length(names)
        raw{j} = single(imread([rawdir,names{j},num2str(1),'.tif']));
        if settings.bgcmoscorrection
            raw{j} = (raw{j}-bgcmos);
            raw{j}(raw{j}<0) = 0;
        end
        if settings.postbin
            raw{j} = imresize(raw{j},.5);
        end
        if settings.bias(j)
            raw{j} = raw{j}./bias_cell{j};
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
    nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
    % nuc_mask=excludelargeandwarped_3(nuc_mask,boulderarea/2,0.95);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Only for DeBug Mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if debug_mode
        keyboard;
    end
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
            nanmask{j}=imdilate(mask_total,strel('disk',nucr/2));
            nanmaskcyto{j}=imdilate(mask_total,strel('disk',nucr*2));
        else
            nanmask{j}=imdilate(whole_mask,strel('disk',nucr/2));
            nanmaskcyto{j}=imdilate(whole_mask,strel('disk',nucr*2));
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
            load([datadir,settings.bleedthroughpth{j}],'bleedthroughrate');
            real1bleedthrough = real{j-1}*bleedthroughrate(2)+bleedthroughrate(1);
            real{j} = real{j}-real1bleedthrough;
        end
    end
    
    
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_label=bwlabel(nuc_mask);
    nuc_info=struct2cell(regionprops(nuc_mask,real{1},'Area','Centroid','MeanIntensity')');
    nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
    nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
    nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
    
    %%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mass=nuc_density.*nuc_area;
    %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [reljitx,reljity]=registerimages(nuc_mask_prev,nuc_mask);
    reljitter=[reljitx,reljity];
    prevjitter=jitters(totalframes,:);
    IFjitter=prevjitter+reljitter;
    nuc_center(:,1)=nuc_center(:,1)+IFjitter(1);
    nuc_center(:,2)=nuc_center(:,2)+IFjitter(2);
    %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    initdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
    %%% track & correct merges (update centers, masses and labels) %%%%%%%%%%%%
    debugpackage={rawprev,nuc_mask_prev,nuc_mask,prevjitter,reljitter};
    %[tracked,nuc_label]=adaptivetrack_IF(tracedata(:,totalframes,1:4),initdata,nuc_label,nucr,debugpackage);
    distthresh=2*nucr; arealowthresh=-.4; areahighthresh=.5;
    [tracked,nuc_label]=adaptivetrack_IF_1(tracedata(:,totalframes,1:4),initdata,nuc_label,nucr,distthresh,arealowthresh,areahighthresh,debugpackage);
    lost_cells = sum(~isnan(tracedata(:,totalframes,1))) - sum(~isnan(tracked(:,1)));
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if maskwrite
        extractmask=bwmorph(nuc_label,'remove');
        imwrite(uint16(extractmask),NEfile);
    end
    
    %%% Measure features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numcells=numel(tracked(:,1));
    nuc_center=tracked(:,[1 2]);
    nuc_area=tracked(:,3);
    nuc_mass=tracked(:,4);
    cellid=find(~isnan(tracked(:,1)));
    numlivecells=numel(cellid);
    nuc_info=regionprops(nuc_label,'PixelIdxList','Centroid');
    
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
        real_temp = real{j};
        real_temp_masked = real_temp;
        real_temp_masked(nanmask{j}) = NaN;
        
        for n = 1:numlivecells
            cc=cellid(n);
            sigmean{j}(cc) = mean(real_temp(nuc_info(cc).PixelIdxList));
            sigmedian{j}(cc) = median(real_temp(nuc_info(cc).PixelIdxList));
            IF_position(cc,:) = nuc_info(cc).Centroid;
            if settings.localbg(j)
                localbg{j}(cc) = getlocalbg(real_temp_masked, nuc_info(cc).Centroid,100, 25, .15);
            end
            if settings.ringcalc(j)
                if cc>numel(ring_info)
                    break;
                end
                ringall = real_temp(ring_info(cc).PixelIdxList);
                ringall(ringall>prctile(ringall,98)) = [];
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
    
    
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IFdata = [IF_position(:,1),IF_position(:,2),nuc_area,...
        cell2mat(sigmean),cell2mat(sigmedian),cell2mat(localbg),...
        cell2mat(sigringmedian)];
    save([datadir,'IF_',shot,'.mat'],'IFdata');
    header = output_names_IF(settings.signals,settings.localbg, settings.ringcalc);
    save([settings.data_path,'\Data\settings_IF.mat'],'settings','header'); %saves all the settings to the Data folder. will overwrite with the most recent
    
    toc(timetotal);
end
%{
%%%%%%%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%%%
%%% view images %%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_label,'remove');
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray(raw{1}));
tempframe(:,:,2)=extractmask;
%tempframe(:,:,3)=marker_mask;
figure,imshow(tempframe);

%%% overlay nuclear images %%%%%%%%%%%%%%%%%
[rawprevcrop,raw1crop]=cropboth(rawprev,raw{1},reljitx,reljity);
tempframe=imadjust(mat2gray(rawprevcrop));
tempframe(:,:,2)=imadjust(mat2gray(raw1crop));
tempframe(:,:,3)=0;
figure,imshow(tempframe)

%%% view images %%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_label,'remove');
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray(raw{1}));
tempframe(:,:,2)=extractmask;
%tempframe(:,:,3)=marker_mask;
figure,imshow(tempframe);


%}
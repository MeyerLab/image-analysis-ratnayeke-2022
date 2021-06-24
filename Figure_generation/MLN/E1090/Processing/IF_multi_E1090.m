function IF_multi_E1090(row,col,site)
% row = 1;
% col = 1;
% site = 1;
debug_mode = 0;

%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Paths
settings.imagingSessions = {'G:\Data\E1090-IF25\', 'G:\Data\E1090-preimage25\'};
settings.dataPath = 'G:\Data\E1090-full\Data\';
settings.bgcmospath = 'C:\Users\Meyerlab\OneDrive - Leland Stanford Junior University\Meyer\MATLAB\BGimages\BG_bin1.mat';
settings.maskpath = 'G:\Data\E1090-full\Mask\';

%%% File Options
settings.separatedirectories_option = 1; %1 = Folder for each well,col,site; 0 = all files in 1 folder
settings.maskwrite_option = 0; %1 = save an image with the nucmask; 0 = dont save an image with the nucmask

%%% Experiment parameters
settings.magnification = 20; %10 = 10x or 20 = 20x
settings.binsize = 1; %1 = bin 1 or 2 = bin 2
settings.postbin = 1;
settings.signals = {{'DAPI_','Cy3_','FarRed_'}, {'FarRed_', 'YFP_', 'RFP_', 'CFP_'}};

%%% Options
settings.ringcalc = {[0 0 0 0], [0 1 1 0]};
settings.localbg = {[1 1 1], [0 0 0 1]};
settings.ringthresh = {[0 0 0 0], [0 100 100 0]};
settings.bias = {[0 0 0], [0 0 0 0]};
settings.biasall = {[0 0 0], [0 0 0 0]};
settings.bleedthrough =  {[0 0 0], [0 0 0 0]};
settings.bleedthroughpth = {{'','',''}, {'', '', '',''}};
settings.bgcmoscorrection = 1; %1 = correct for shading; 0 = dont correct for shading;
settings.segmethod = {'double marker','thresh'}; %'log' or 'single' or 'double'
settings.bgsubmethod = {{'global nuclear','global nuclear','global nuclear'}, 
    {'global nuclear','global cyto','global cyto', 'global nuclear'}}; %'global nuclear','global cyto','tophat','semi-local nuclear'

%%% Segmentation parameters
settings.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100; %100
settings.boulderarea = 2000; %
settings.blobthreshold = -0.03;
settings.blurradius = 3;
settings.remove_smears_option = 0;
settings.soliditythresh = 0.7;
settings.compression = 4;
settings.exclude = {};


%%% Designate which Well is being analyzied %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datadir = settings.dataPath;

if ~exist(datadir,'dir')
    mkdir(datadir);
end
if any(ismember(settings.exclude,shot))
    IFdata = ones(1,16)*NaN;
    save([datadir,'IFdata_',shot,'.mat'],'IFdata');
    return
end
if ~exist([datadir,'IFdata_',shot,'.mat']);
    %% Initialize variables
    separatedirectories = settings.separatedirectories_option;
    if separatedirectories
        maskdir = [settings.maskpath,shot];
    else
        maskdir = [settings.maskpath,'\',shot,'_'];
    end
    if ~exist(maskdir,'dir') && settings.maskwrite_option
        mkdir(maskdir);
    end
    
    %%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    names = settings.signals;
    nucname = names{1}{1};
    
    %%% Segmentation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucr = settings.nucr;
    debrisarea = settings.debrisarea;
    boulderarea = settings.boulderarea;
    blobthreshold = settings.blobthreshold;
    solidity = settings.soliditythresh;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timetotal = tic;
    %%% bgcmos correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([settings.bgcmospath],'cmosoffset');
    bgcmos = cmosoffset;
    [height,width]=size(bgcmos);
    if settings.postbin
        height = height/2;
        width = width/2;
    end
    for i = 1:length(names)
        for j = 1:length(names{i})
            if settings.biasall{i}(j)
                biasdir = [settings.imagingSessions{i} 'Bias_all\'];
            else
                biasdir = [settings.imagingSessions{i} 'Bias\'];
            end
            if settings.bias{i}(j)
                load([biasdir,names{i}{j},num2str(site),'.mat']);
                bias_cell{i}{j} = bias;
            end
        end
    end
    
    %% load images and segment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:length(names)
        if separatedirectories
            rawdir = [settings.imagingSessions{i},'Raw\',shot,'\',shot,'_'];
        else
            rawdir = [settings.imagingSessions{i},'Raw\',shot,'_'];
        end
        
        for j = 1:length(names{i})
            raw{i}{j} = single(imread([rawdir,names{i}{j},num2str(1),'.tif']));
            if settings.bgcmoscorrection
                raw{i}{j} = (raw{i}{j}-bgcmos);
                raw{i}{j}(raw{i}{j}<0) = 0;
            end
            if settings.postbin
                raw{i}{j} = imresize(raw{i}{j},.5);
            end
            if settings.bias{i}(j)
                raw{i}{j} = raw{i}{j}./bias_cell{i}{j};
            end
            normraw{i}{j} = (raw{i}{j}-min(raw{i}{j}(:)))./(max(raw{i}{j}(:))-min(raw{i}{j}(:)));
        end
        
        %%%remove smears %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if settings.remove_smears_option
            foregroundthresh = 1000;
            areathresh = 20000;
            bgprctile = 20;
            raw1 = removesmears_1(raw1,foregroundthresh,areathresh,bgprctile);
            raw2 = removesmears_1(raw2,foregroundthresh,areathresh,bgprctile);
            if settings.signal3_option;raw3 = removesmears_1(raw3,foregroundthresh,areathresh,bgprctile);end
            if settings.signal4_option;raw4 = removesmears_1(raw4,foregroundthresh,areathresh,bgprctile);end
            if settings.signal5_option;raw5 = removesmears_1(raw5,foregroundthresh,areathresh,bgprctile);end
        end
        %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        segmethod = settings.segmethod{i};
        blurradius = settings.blurradius;
        switch segmethod
            case 'log'
                nuc_mask{i} = blobdetector_4(log(raw{i}{1}),nucr,blobthreshold,debrisarea);
            case 'thresh'
                nuc_mask{i} = threshmask(raw{i}{1},blurradius);
            case 'single'
                nuc_mask{i} = threshmask(raw{i}{1},blurradius);
                nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*2/3),4);
            case 'double'
                nuc_mask{i} = threshmask(raw{i}{1},blurradius);
                nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*2/3),4);
                nuc_mask{i} = secondthresh(raw{i}{1},blurradius,nuc_mask{i},boulderarea*2);
            case 'double marker'
                nuc_mask{i} = threshmask(raw{i}{1},blurradius);
                nuc_mask{i} = secondthresh_all(raw{i}{1},blurradius,nuc_mask{i});
                nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*.4),4);
            case 'multithresh'
                nuc_mask{i} = threshmask_multi(raw{i}{1},blurradius,4, 1e-04);
                nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*.4),4);
                %nuc_mask{i} = secondthresh(raw{i}{1},blurradius,nuc_mask{i},boulderarea);
        end
        %%% generate overlay mask foreground (nuc+channel)
        nuc_mask{i} = bwareaopen(nuc_mask{i},debrisarea);
        %%% Deflection-Bridging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i==1
            nuc_mask{i} = segmentdeflections_bwboundaries(nuc_mask{i},nucr,debrisarea);
            nuc_mask{i} = excludelargeandwarped_3(nuc_mask{i},boulderarea,solidity);
        end
        %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask{i} = imclearborder(nuc_mask{i});
    end
    %% align images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jitmatx = [];
    jitmaty = [];
    for i = 1:length(names)
        [reljitx{i},reljity{i}] = registerimages(nuc_mask{1},nuc_mask{i});
        jitmatx = [jitmatx; reljitx{i}];
        jitmaty = [jitmaty; reljity{i}];
    end
    cropcoors = getcropcoors([height width], jitmatx, jitmaty);
    for i = 1:length(names)
        for j = 1:length(names{i})
            raw{i}{j} = raw{i}{j}(cropcoors(i+1,1):cropcoors(i+1,2),cropcoors(i+1,3):cropcoors(i+1,4));
            if settings.localbg{i}(j)
                overlayraw = normraw{i}{1} + normraw{i}{j};
                overlayraw = overlayraw(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
                foreground_aligned{i}{j} = threshmask_adapt(overlayraw,blurradius);
            end
        end

    end
    %clear normraw;
    %foreground_aligned{i} = threshmask_adapt(raw{i}{1},blurradius);
    nuc_mask_aligned = nuc_mask{1}(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
    
    
    %% calculate background: semi-local %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compression = settings.compression;
    for i=1:length(names)
        for j = 1:length(names{i})
            if settings.localbg{i}(j)
                nanmask_aligned{i}{j} = imdilate(foreground_aligned{i}{j},strel('disk',nucr/2));
                nanmaskcyto_aligned{i}{j} = imdilate(foreground_aligned{i}{j},strel('disk',nucr*2));
            else
                nanmask_aligned{i}{j} = imdilate(nuc_mask_aligned,strel('disk',nucr/2));
                nanmaskcyto_aligned{i}{j} = imdilate(nuc_mask_aligned,strel('disk',nucr*2));
            end
            %blur = imfilter(raw{j},fspecial('disk',3),'symmetric');
            blur = raw{i}{j};
            switch settings.bgsubmethod{i}{j}
                case 'global nuclear'
                    real{i}{j} = bgsubmasked_global_NR(blur,nanmask_aligned{i}{j},1,compression,25);
                case 'global cyto'
                    real{i}{j} = bgsubmasked_global_2(blur,nanmaskcyto_aligned{i}{j},1,compression,25);
                case 'tophat'
                    real{i}{j} = imtophat(blur,strel('disk',nucr,0));
                case 'semi-local nuclear'
                    real{i}{j} = bgsubmasked_global_2(blur,nanmask_aligned{i}{j},11,compression,25);
                case 'none'
                    real{i}{j} = blur;
            end
        end
        %clear foreground_aligned;
        %%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:length(names{i})
            if settings.bleedthrough{i}(j)
                load([datadir,settings.bleedthroughpth{i}{j}],'bleedthroughrate');
                real1bleedthrough = real{i}{j-1}*bleedthroughrate(2)+bleedthroughrate(1);
                real{i}{j} = real{i}{j}-real1bleedthrough;
            end
        end
    end
    
    %% debugging: view images %%%%%%%%%%
    %{
        extractmask = bwmorph(foreground_aligned{1}{4},'remove');
        tempframe = imadjust(mat2gray(raw{1}{4}));
        tempframe(:,:,2) = extractmask;
        tempframe(:,:,3) = 0;
        figure,imshow(tempframe);
    
        extractmask = bwmorph(nuc_mask{1},'remove');
        tempframe = imadjust(mat2gray(raw{1}{1}));
        tempframe(:,:,2) = extractmask;
        tempframe(:,:,3) = 0;
        figure,imshow(tempframe);
        
        rawnucmask = raw{1};
        rawnucmask(~nuc_mask) = NaN;
        figure, imagesc(rawnucmask);
        
        nuc_info = struct2cell(regionprops(nuc_mask,'Area')');
        nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
        hist(nuc_area,100);
        
        anti_mask = bwareaopen(nuc_mask,debrisarea);
        temp_mask = nuc_mask-anti_mask;
        extractmask = bwmorph(temp_mask,'remove');
        
        anti_mask = bwareaopen(smear_mask,10000);
        extractmask = bwmorph(anti_mask,'remove'); 
        
        %% check masks of each session
        for session = 1:length(names)
            extractmask = bwmorph(nuc_mask_aligned,'remove');
            tempframe = imadjust(mat2gray(real{session}{1}));
            tempframe(:,:,2) = extractmask;
            tempframe(:,:,3) = 0;
            figure,imshow(tempframe);
        end
        
        %%% overlay nuclear images %%%%%%%%%%%%%%%%%
        tempframe=imadjust(mat2gray(real{1}{1}));
        tempframe(:,:,2)=imadjust(mat2gray(real{2}{1}));
        tempframe(:,:,3)=0;
        figure,imshow(tempframe)
    %}
    %%% Only for DeBug Mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if debug_mode
        keyboard;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_info = struct2cell(regionprops(nuc_mask_aligned,'Area','Centroid')');
    nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
    nuc_center = squeeze(cell2mat(nuc_info(2,1,:)))';
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask = bwmorph(nuc_mask_aligned,'remove');
    if settings.maskwrite_option
        imwrite(uint16(extractmask),NEfile);
    end
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_label = bwlabel(nuc_mask_aligned);
    numcells = numel(nuc_area);
    nuc_info = regionprops(nuc_label,'PixelIdxList','Centroid');
    
    %initialize ring variables
    if any(cell2mat(settings.ringcalc))
        if (settings.magnification == 10 & settings.binsize == 1) | (settings.magnification == 20 & settings.binsize == 2)
            innerrad = 1; outerrad = 5;
        else
            innerrad = 2; outerrad = 10;
        end
        ring_label = getcytoring_thicken(nuc_label,innerrad,outerrad,real{2});
        ring_info = regionprops(ring_label,'PixelIdxList');
    end
    
    %Initialize storage variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nanvec = ones(numcells,1)*NaN;
    sigringmedian{1} = [];
    localbg{1} = [];
    for i = 1:length(names)
        for j = 1:length(names{i})
            sigmean{i}{j} = nanvec;
            sigmedian{i}{j} = nanvec;
            %sigring75th{j} = nanvec;
            if settings.ringcalc{i}(j)
                sigringmedian{i}{j} = nanvec;
            end
            if settings.localbg{i}(j)
                localbg{i}{j} = nanvec;
            end
            
            
            
            %Measure features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            real_temp = real{i}{j};
            real_temp_masked = real_temp;
            real_temp_masked(nanmask_aligned{i}{j}) = NaN;
            for cc = 1:numcells
                sigmean{i}{j}(cc) = mean(real_temp(nuc_info(cc).PixelIdxList));
                sigmedian{i}{j}(cc) = median(real_temp(nuc_info(cc).PixelIdxList));
                if settings.localbg{i}(j)
                    localbg{i}{j}(cc) = getlocalbg(real_temp_masked, nuc_info(cc).Centroid,100, 50,.15);
                end
                if settings.ringcalc{i}(j)
                    if cc>numel(ring_info)
                        break;
                    end
                    ringall = real_temp(ring_info(cc).PixelIdxList);
                    ringall(ringall>prctile(ringall,98)) = [];
                    %sigring75th{j}(n) = prctile(ringall,75);
                    ringforeground = ringall(ringall>settings.ringthresh{i}(j));
                    if numel(ringforeground)<100
                        ringforeground = ringall;
                    end
                    if numel(ringall)>100
                        sigringmedian{i}{j}(cc) = nanmedian(ringforeground);
                    end
                end
            end
        end
    end
    
    
    %% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IFdata = [nuc_center(:,1),nuc_center(:,2),nuc_area,...
        cell2mat([sigmean{:}]),cell2mat([sigmedian{:}]),cell2mat([localbg{:}]),...
        cell2mat([sigringmedian{:}])];
    save([datadir,'IFdata_',shot,'.mat'],'IFdata');
    header = output_names_multiIF(settings.signals,settings.localbg, settings.ringcalc);
    save([datadir,'settings.mat'],'settings','header'); %saves all the settings to the Data folder. will overwrite with the most recent
    toc(timetotal);
end
close all;
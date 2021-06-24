function Timelapse_C150(row,col,site,debug_mode)
% row = 2;
% col = 2;
% site = 1;
% debug_mode = 1;
timetotal = tic;

%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Paths
experiment_name='C150-live';
image_drive = 'I:\4TB3\Data\';
savepath=['F:\Data\C-Cdt1\',experiment_name,'\Data\'];
imagepath=[image_drive,experiment_name,'\Raw\'];
biaspath=[image_drive,experiment_name,'\Bias\'];
maskpath=[image_drive,experiment_name,'\Mask\'];
bgcmospath='F:\MATLAB\BGimages\IX_XL\cmosoffset_bin1.mat';

%%% General parameters
startFrame=1;   %Frame to start analyzing from
endFrame=89;  %Frame to stop analyzing
magnification=10; %10=10x or 20=20x
binsize=1; %1=bin 1 or 2=bin 2
signals={'CFP_','YFP_','RFP_'};
maskwrite = 1; %1 = save an image with the nucmask; 0 = dont save an image with the nucmask
maskname = 'nucedge_';

%%% Quantification parameters
bgcmoscorrection = 1;
bias = [1 1 1];
signal_foreground = [1 1 1];
bgsubmethod = {'global nuclear','global nuclear','global nuclear'}; %Options:'global nuclear','global cyto','tophat','semi-local nuclear', 'none'
compression = 4;
sigblur = [3 3 3];
localbgmeasure = [0 1 1];
ringcalc = [0 1 0];
ringthresh = [0 50 50];
punctacalc = 0;
punctaThresh = [175 200 225];
varThresh = [75 100 125];

%%% Segmentation parameters
firstsegmethod = 'concavity'; %Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
secondsegmethod = 'concavity'; %Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
blurradius  =  3; %10x: 3
soliditythresh = 0.8;
debrisarea = 100; % 10x = 100
boulderarea = 1500; % 10x = 1500
blobthreshold = -0.02;

%%% Tracking parameters
maxjump = nucr*3;
masschangethreshold = 0.30;
areachangethreshold = 0.60;
daughtervariance = 0.10;

%% TRACKING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
if ~exist([savepath,'tracedata_',shot,'.mat'],'file')
    
    %% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize path variables
    rawdir = [imagepath,shot,'\',shot,'_'];
    maskdir = [maskpath,shot,'\',shot,'_'];
    nucname = signals{1};
    parameternum = 4+2*length(signals)+sum(ringcalc)+sum(localbgmeasure);
    if punctacalc
        parameternum = parameternum + 4 + 2*length(punctaThresh) + 2*length(varThresh);
    end
    
    if ~exist(savepath,'dir')
        mkdir(savepath);
    end
    if ~exist(maskdir,'dir') && maskwrite
        mkdir(maskdir);
    end
    
    %%% Initialize tracking variables
    frames = startFrame:endFrame;
    totalframes = numel(frames);
    badframes = ones(endFrame,1)*NaN;
    if startFrame>1
        badframes(1:startFrame-1) = 0;
    end
    jitters = zeros(endFrame,2);
    blocksize = 10000;
    maxcellnum = blocksize;
    tracedata = ones(maxcellnum,endFrame,parameternum)*NaN;
    tracking = ones(maxcellnum,5)*NaN;
    [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width] = ...
        timelapsesetup_4(rawdir,nucname,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
    regheight = 1:0.5*height; regwidth = 1:0.5*width;
    segparams = struct('nucr', nucr,'blurradius', blurradius,'debrisarea', debrisarea,...
        'boulderarea',boulderarea, 'blobthreshold', blobthreshold);
    trackparams = {nucr,maxjump,debrisarea,masschangethreshold,areachangethreshold,daughtervariance};
    
    %%% Load bgcmos/bias
    load(bgcmospath,'cmosoffset');
    bgcmos  =  cmosoffset;
    for i = 1:length(signals)
        if bias(i)
            load([biaspath,signals{i},num2str(site),'.mat']); 
            bias_cell{i} = bias;
        end
    end
    
    %% Imaging processing for each frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = firstgoodindex:totalframes
        f = frames(i);
        fprintf('frame %0.0f\n',f);
        timeframe = tic;
        
        %%% Load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:length(signals)
            raw{j} = double(imread([rawdir,signals{j},num2str(f),'.tif']));
            if bgcmoscorrection
                raw{j} = (raw{j}-bgcmos);
            end
            if bias(j)
                raw{j} = raw{j}./bias_cell{j};
            end
            if signal_foreground(j)
               %foreground{j} = threshmask_adapt(raw{j},blurradius);
               foreground{j} = threshmask(raw{j},blurradius);
               if sum(foreground{j}(:))/length(foreground{j}(:)) > .9 % account for large foregrounds
                   foreground{j} = zeros(height, width);
               end
            end
        end
        
        %%% Segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i == firstgoodindex
            nuc_mask = segmentImage(raw{1},firstsegmethod, segparams);
        else
            if isequal(secondsegmethod, 'apriori')
                nuc_mask = threshmask(raw{1});
            else
                nuc_mask = segmentImage(raw{1},secondsegmethod, segparams);
            end
            %%% Calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lastgoodframe = find(badframes == 0,1,'last');
            [reljitx,reljity] = registerimages(imfill(extractmask(regheight,regwidth),'holes'),nuc_mask(regheight,regwidth));
            jitters(f,:) = jitters(lastgoodframe,:)+[reljitx,reljity];
            if isequal(secondsegmethod, 'apriori')
                [nuc_mask,marker_mask] = apriori_markermask(nuc_mask,nuc_center,jitters(f,:));
            end
        end
        %nuc_mask = imopen(nuc_mask,strel('disk',3));
        whole_mask = nuc_mask;
      
        %%% Clean up segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask = bwareaopen(nuc_mask,debrisarea);
        nuc_mask = segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
        nuc_mask = excludelargeandwarped_3(nuc_mask,boulderarea,soliditythresh);
        nuc_mask = imclearborder(nuc_mask);
        
        %%% Check segmentation
        if debug_mode
            keyboard;
        end
        %{
        %%% debugging: view images %%%%%%%%%%
        %Check nuclear segmentation
        extractmask = bwmorph(nuc_mask,'remove');
        tempframe = zeros(height,width,3);
        tempframe(:,:,1) = imadjust(mat2gray((raw{1})));
        tempframe(:,:,2) = extractmask;;
        %tempframe(:,:,3) = marker_mask;
        figure,imshow(tempframe);

        %Check foreground segmentation
        extractmask = bwmorph(foreground{2},'remove');
        tempframe = zeros(height,width,3);
        tempframe(:,:,1) = imadjust(mat2gray((raw{1})));
        tempframe(:,:,2) = extractmask;;
        %tempframe(:,:,3) = marker_mask;
        figure,imshow(tempframe);

        nuc_info = struct2cell(regionprops(nuc_mask,'Area')');
        nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
        hist(nuc_area,100);

        anti_mask = bwareaopen(nuc_mask,debrisarea);
        temp_mask = nuc_mask-anti_mask;
        extractmask = bwmorph(temp_mask,'remove');

        anti_mask = bwareaopen(nuc_mask,3500);
        extractmask = bwmorph(anti_mask,'remove');
        %}
        
        %%% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:length(signals)
            if signal_foreground(j)
                mask_total  =  whole_mask | foreground{j};
                nanmask{j} = imdilate(mask_total,strel('disk',nucr));
                nanmaskcyto{j} = imdilate(mask_total,strel('disk',nucr*2));
            else
                nanmask{j} = imdilate(whole_mask,strel('disk',nucr));
                nanmaskcyto{j} = imdilate(whole_mask,strel('disk',nucr*2));
            end
            
            %%% Blur signal
            if sigblur > 0
                blurIm = imfilter(raw{j},fspecial('disk',sigblur(j)),'symmetric');
            else
                blurIm  =  raw{j};
            end
            
            %%% Background subtract signals
            switch bgsubmethod{j}
                case 'global nuclear'
                    real{j} = bgsubmasked_global_NR(blurIm,nanmask{j},1,compression,25);
                    unblurred{j} = bgsubmasked_global_NR(raw{j},nanmask{j},1,compression,25);
                case 'global cyto'
                    real{j} = bgsubmasked_global_2(blurIm,nanmaskcyto{j},1,compression,25);
                    unblurred{j} = bgsubmasked_global_2(raw{j},nanmaskcyto{j},1,compression,25);
                case 'tophat'
                    real{j} = imtophat(blurIm,strel('disk',nucr,0));
                    unblurred{j} = imtophat(raw{j},strel('disk',nucr,0));
                case 'semi-local nuclear'
                    real{j} = bgsubmasked_global_2(blurIm,nanmask{j},11,compression,25);
                    unblurred{j} = bgsubmasked_global_2(raw{j},nanmask{j},11,compression,25);
                case 'none'
                    real{j} = blurIm;
                    unblurred{j} = raw{j};
            end
        end
        
        %%% Extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [nuc_label,numcells] = bwlabel(nuc_mask);
        nuc_info = struct2cell(regionprops(nuc_mask,real{1},'Area','Centroid','MeanIntensity')');
        nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
        
        %%% detect bad frame 
        mednuc = median(nuc_area);
        if i>firstgoodindex+1 && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)
            fprintf('badframe: frame %0.0f\n',f);
            badframes(f) = 1;
            extractmask = bwmorph(nuc_mask,'remove');
            if maskwrite
                imwrite(uint16(extractmask),[maskdir,maskname,num2str(f),'.tif']);
            end
            continue;
        end
        blurthreshhigh = 1.2*mednuc;
        blurthreshlow = 0.8*mednuc;
        numthresh = 0.5*numcells;
        nuc_center = squeeze(cell2mat(nuc_info(2,1,:)))';
        nuc_density = squeeze(cell2mat(nuc_info(3,1,:)));
        
        %%% calculate masses 
        nuc_mass = nuc_density.*nuc_area;
        curdata = [nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        
        %%% Track cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i>firstgoodindex
            nuc_center(:,1) = nuc_center(:,1)+jitters(f,1);
            nuc_center(:,2) = nuc_center(:,2)+jitters(f,2);
            %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            curdata = [nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
            debugpackage = {extractmask,jitters(lastgoodframe,:),[reljitx,reljity]};
            %%% track & correct merges (update centers, masses and labels) %%%%
            [tracedata,curdata,tracking,nuc_label] = ...
                adaptivetrack_9(f,lastgoodframe,f,tracedata,curdata,tracking,real{1},nuc_label,jitters(f,:),trackparams,debugpackage);
            badframes(f) = 0;
            
            %%% Check tracking
            if debug_mode
                keyboard;
            end
            %{
            matchinglabel=prev_label.*nuc_label;
            tempframe=zeros(height,width,3);
            tempframe(:,:,1)=prev_label>0;
            tempframe(:,:,2)=nuc_label>0;
            tempframe(:,:,3)=matchinglabel;
            figure,imshow(tempframe);
            %}
        end
        
        %%% Save mask
        extractmask = bwmorph(nuc_label,'remove');
        if maskwrite
            imwrite(uint16(extractmask),[maskdir,maskname,num2str(f),'.tif']);
        end
        prev_label = nuc_label;
        
        %%% Reformat features 
        cellid = find(~isnan(curdata(:,1)));
        numlivecells = numel(cellid);
        curdata = curdata(cellid,:);
        nuc_center = curdata(:,[1 2]);
        nuc_area = curdata(:,3);
        nuc_mass = curdata(:,4);
        nuc_info = regionprops(nuc_label,'PixelIdxList','Centroid', 'BoundingBox');
        nanvec = ones(numlivecells,1)*NaN;
        
        %% Measure features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% Initialize measurement variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sigringmedian{1} = [];
        localbg{1} = [];
        
        % PCNA variables
        PCNAMean = [];
        filtMean = [];
        varMean = [];
        varStd = [];
        filtMaskedMean{1} = [];
        filtMaskedArea{1} = [];
        varMaskedMean{1} = [];
        varMaskedArea{1} = [];
        
        for j = 1:length(signals)
            sigmean{j} = nanvec;
            sigmedian{j} = nanvec;
            if ringcalc(j)
                sigringmedian{j} = nanvec;
            end
            if localbgmeasure(j)
                localbg{j} = nanvec;
            end
        end
        
        if any(ringcalc)
            if (magnification == 10 && binsize == 1) || (magnification == 20 && binsize == 2)
                innerrad = 1; outerrad = 5;
            else
                innerrad = 2; outerrad = 10;
            end
            ring_label = getcytoring_thicken(nuc_label,innerrad,outerrad,real{2});
            ring_info = regionprops(ring_label,'PixelIdxList');
            
            %%% Check cytoring
            if debug_mode
                keyboard;
            end
            %{
            %%% debugging: view images %%%%%%%%%%
            temp_ring  =  ring_label > 0;
            extractmask = bwmorph(temp_ring,'remove');
            signal  =  real{1};
            intensities  =  signal.*temp_ring;
            tempframe = zeros(height,width,3);
            tempframe(:,:,1) = imadjust(mat2gray(real{1}));
            tempframe(:,:,2) = extractmask;;
            %tempframe(:,:,3) = marker_mask;
            figure,imshow(tempframe);
            %}
        end
        
        %%% Measure for each signal
        for j = 1:length(signals)
            real_temp = real{j};
            real_masked = real_temp;
            real_masked(nanmask{j}) = NaN;
            
            %% Puncta analysis
            if j == punctacalc
               %Preprocess PCNA image
               punctaSig = unblurred{punctacalc};
               punctaBlur = imgaussfilt(punctaSig, 8); %blurred to remove nucleoli
               %punctaSemiBlur = imgaussfilt(punctaSig, 1);
               punctaSemiBlur = punctaSig;
               
               %Masking to remove nucleoli and border
               nuclear_mask = nuc_label > 0;
               punctaDiff = real_temp - punctaBlur;
               erode_nuc = imerode(nuclear_mask, strel('disk',7,0));
               nucleolus_mask = (punctaDiff < -10);
               nucleolus_mask = imdilate(nucleolus_mask, strel('disk',1,0));
               analysis_mask = erode_nuc & ~nucleolus_mask;
               
               %Image transformations quantification
               punctaFilt = imtophat(punctaSig, strel('disk',2,0));
               punctaVar = stdfilt(punctaSemiBlur,ones(3));

               %Mask transformed images
               nanVar = punctaVar;
               nanVar(~analysis_mask) = NaN;
               
               %Threshold transformed images
               for t = 1:length(punctaThresh)
                   filtMask{t} = punctaFilt > punctaThresh(t) & nuclear_mask;
                   punctaFiltMasked{t} = punctaFilt * filtMask{t};
               end
               for t = 1:length(varThresh)
                   varMask{t} = punctaVar > varThresh(t) & analysis_mask;
                   punctaVarMasked{t} = punctaVar * varMask{t};
               end
               
               %Initialize PCNA storage
               PCNAMean = nanvec;
               filtMean = nanvec;
               varMean = nanvec;
               varStd = nanvec;
               for t = 1:length(punctaThresh)
                   filtMaskedMean{t} = nanvec;
                   filtMaskedArea{t} = nanvec;
               end
               for t = 1:length(varThresh)
                   varMaskedMean{t} = nanvec;
                   varMaskedArea{t} = nanvec;
               end
               
               % Check masking/puncta
               if debug_mode
                   keyboard;
               end
               %{
               extractmask = bwmorph(analysis_mask,'remove');
               tempframe=imadjust(mat2gray(punctaVar));
               tempframe(:,:,2)=(extractmask);
               tempframe(:,:,3)=0;
               figure,imshow(tempframe)
               %} 
            end
            
            %% Measurements for each cell
            for n = 1:numlivecells
                cc = cellid(n);
                sigmean{j}(n) = mean(real_temp(nuc_info(cc).PixelIdxList));
                sigmedian{j}(n) = median(real_temp(nuc_info(cc).PixelIdxList));
                if localbgmeasure(j)
                    localbg{j}(n) = getlocalbg(real_masked, nuc_info(cc).Centroid,100, 25, .15);
                end
                if ringcalc(j)
                    if ~(cc > numel(ring_info))
                        ringall = real_temp(ring_info(cc).PixelIdxList);
                        ringall(ringall>prctile(ringall,98)) = [];
                        ringforeground = ringall(ringall>ringthresh(j));
                        if numel(ringforeground)<50
                            ringforeground = ringall;
                        end
                        if numel(ringall)>= 50
                            sigringmedian{j}(n) = nanmedian(ringforeground);
                        end
                    end
                end
                if j == punctacalc
                    PCNAMean(n) = mean(punctaSemiBlur(nuc_info(cc).PixelIdxList));
                    filtMean(n) = mean(punctaFilt(nuc_info(cc).PixelIdxList));
                    varMean(n) = nanmean(nanVar(nuc_info(cc).PixelIdxList));
                    varStd(n) = nanstd(nanVar(nuc_info(cc).PixelIdxList));
                    for t = 1:length(punctaThresh)
                        filtMaskedMean{t}(n) = nanmean(punctaFiltMasked{t}(nuc_info(cc).PixelIdxList));
                        filtMaskedArea{t}(n) = sum(filtMask{t}(nuc_info(cc).PixelIdxList));
                    end
                    for t = 1:length(varThresh)
                        varMaskedMean{t}(n) = nanmean(punctaVarMasked{t}(nuc_info(cc).PixelIdxList));
                        varMaskedArea{t}(n) = sum(varMask{t}(nuc_info(cc).PixelIdxList));
                    end
                end
            end
        end
        
        %% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tracedata(cellid,f,:) = ...
            [nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,...
            cell2mat(sigmean),cell2mat(sigmedian),cell2mat(localbg),...
            cell2mat(sigringmedian), PCNAMean, ...
            filtMean, varMean, varStd, cell2mat(filtMaskedMean), cell2mat(filtMaskedArea), ...
            cell2mat(varMaskedMean), cell2mat(varMaskedArea)];
        if maxcellnum-max(cellid)<blocksize
            tempdata = ones(blocksize,endFrame,parameternum)*NaN;
            temptrack = ones(blocksize,5)*NaN;
            tracedata = [tracedata;tempdata];
            tracking = [tracking;temptrack];
            maxcellnum = maxcellnum+blocksize;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toc(timeframe);
    end
    [tracedata,genealogy,jitters] = postprocessing_nolinking(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
    %%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save([savepath,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
    names = output_names_puncta(signals,localbgmeasure, ringcalc, punctacalc, punctaThresh, varThresh);
    save([savepath,'settings.mat'],'names'); %saves all the settings to the Data folder. will overwrite with the most recent
    
    toc(timetotal);
    clear all;
end
end

function mask = segmentImage(image, segMethod, s)
switch segMethod
    case 'concavity'
        %imblur = imfilter(image,fspecial('gaussian',s.blurradius),'symmetric');
        imblur = imfilter(image,fspecial('gaussian',5),'symmetric');
        mask = ThreshImage(imblur);
        mask = logical(imfill(mask,'holes'));
    case 'log'
        mask = blobdetector_4(log(image),s.nucr,s.blobthreshold,s.debrisarea);
    case 'single adapt'
        mask = threshmask_adapt(image,s.blurradius);
        mask = markershed_filter(mask,round(s.nucr*2/3),6);
    case 'double adapt'
        mask = threshmask_adapt(image,s.blurradius);
        mask = markershed_filter(mask,round(s.nucr*2/3),6);
        mask = secondthresh(image,s.blurradius,mask,s.boulderarea*2);
    case 'single'
        mask = threshmask(image,s.blurradius);
        mask = markershed_filter(mask,round(s.nucr*2/3),6);
    case 'double'
        mask = threshmask(image,s.blurradius);
        mask = markershed_filter(mask,round(s.nucr*2/3),6);
        mask = secondthresh(image,s.blurradius,mask,s.boulderarea*2);        
end
end
function Timelapse_E1099(row,col,site)
% row = 3;
% col = 1;
% site = 1;
debug_mode=0;

%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Paths
settings.experiment_name='E1099-live';
settings.projectpath=['G:\Data\',settings.experiment_name,'\Data\'];
settings.imagepath=['G:\Data\',settings.experiment_name,'\Raw\'];
settings.bgcmospath='C:\Users\Meyerlab\OneDrive - Leland Stanford Junior University\Meyer\MATLAB\BGimages\BG_bin1.mat';
settings.biaspath=['G:\Data\',settings.experiment_name,'\Bias\'];
settings.maskpath=['G:\Data\',settings.experiment_name,'\Mask\'];
%%% Experiment parameters
settings.StartFrame=2;   %Frame to start analyzing from
settings.EndFrame=96;  %Frame to stop analyzing  %must be greater than StartFrame
settings.magnification=10; %10=10x or 20=20x
settings.binsize=1; %1=bin 1 or 2=bin 2
settings.signals={'CFP_','YFP_','RFP_'};

%%% Options
settings.maskwrite_option=1; %1=save an image with the nucmask; 0=dont save an image with the nucmask
settings.ringcalc=[0 0 0];
settings.localbg=[0 1 1];
settings.ringthresh=[0 100 100];
settings.bias=[1 1 1];
settings.bgcmoscorrection=1; %1=correct for shading; 0=dont correct for shading;
settings.firstsegmethod='double'; %'log' or 'single' or 'double'
settings.secondsegmethod='apriori'; %'log' or 'double' or 'apriori'
settings.bgsubmethod={'global nuclear','global nuclear','global nuclear'}; %'global nuclear','global cyto','tophat','semi-local nuclear'

%%% Segmentation parameters
settings.nucr=12; %10x bin1=12 and 20x bin2=12
settings.debrisarea=100;
settings.boulderarea=1500;
settings.blobthreshold=-0.03;

%%% Tracking parameters
settings.maxjump=settings.nucr*3;
settings.masschangethreshold=0.30;
settings.areachangethreshold=0.60;
settings.daughtervariance=0.10;
settings.blurradius=3;
settings.blocksize=10000;
settings.soliditythresh=0.8;
settings.compression=4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Designate which Well is being analyzied %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
namenucedge='nucedge_';
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datadir=settings.projectpath;
imagepath=settings.imagepath;
bgcmospath=settings.bgcmospath;
if ~exist(datadir,'dir')
    mkdir(datadir);
end
if ~exist([datadir,'tracedata_',shot,'.mat']);
    rawdir=[settings.imagepath,shot,'\',shot,'_'];
    biasdir=[settings.biaspath];
    maskdir=[settings.maskpath,shot];
    maskwrite=settings.maskwrite_option;
    if ~exist(maskdir,'dir') && maskwrite
        mkdir(maskdir);
    end
    maskdir=[maskdir,'\',shot,'_'];
    
    %%% General settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SF=settings.StartFrame;EF=settings.EndFrame;
    names=settings.signals;
    nucname=names{1};
    parameternum=4+2*length(settings.signals)+sum(settings.ringcalc)+sum(settings.localbg);
    %%% Segmentation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucr=settings.nucr;
    debrisarea=settings.debrisarea;
    boulderarea=settings.boulderarea;
    blobthreshold=settings.blobthreshold;
    blurradius=settings.blurradius;
    %%% Tracking parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trackparams={nucr,settings.maxjump,debrisarea,settings.masschangethreshold,settings.areachangethreshold,settings.daughtervariance};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    frames=SF:EF;
    totalframes=numel(frames);
    badframes=ones(EF,1)*NaN;
    if SF>1
        badframes(1:SF-1)=0;
    end
    jitters=zeros(EF,2);
    blocksize=settings.blocksize;
    maxcellnum=blocksize;
    tracedata=ones(maxcellnum,EF,parameternum)*NaN;
    tracking=ones(maxcellnum,5)*NaN;
    timetotal=tic;
    [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=...
        timelapsesetup_4(rawdir,nucname,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
    regheight=1:0.5*height; regwidth=1:0.5*width;
    
    if any(settings.ringcalc)
        if (settings.magnification==10 & settings.binsize==1) | (settings.magnification==20 & settings.binsize==2)
            innerrad=1; outerrad=5;
        else
            innerrad=2; outerrad=10;
        end
    end
    %%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([bgcmospath],'cmosoffset');
    bgcmos = cmosoffset;
    for i=1:length(names)
        if settings.bias(i)
            load([biasdir,names{i},num2str(site),'.mat']); bias_cell{i}=bias;
        end
    end
    
    %%% Imaging processing, Frame by Frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=firstgoodindex:totalframes
        f=frames(i); fprintf('frame %0.0f\n',f);
        timeframe=tic;
        %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j=1:length(names)
            raw{j}=double(imread([rawdir,names{j},num2str(f),'.tif']));
            if settings.bgcmoscorrection
                raw{j}=(raw{j}-bgcmos);
            end
            if settings.bias(j)
                raw{j}=raw{j}./bias_cell{j};
            end
            normraw{j}=(raw{j}-min(raw{j}(:)))./(max(raw{j}(:))-min(raw{j}(:)));
            
            if settings.localbg(j)
                overlay = normraw{j} + normraw{1};
                foreground{j} = threshmask_adapt(overlay,blurradius);
                %foreground{j}=nuc_mask;
%                 if sum(foreground{j}(:))/length(foreground{j}(:)) > .9
%                     foreground{j} = zeros(height, width);
%                 end
            end
        end
        
        %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i==firstgoodindex
            switch settings.firstsegmethod
                case 'log'
                    nuc_mask=blobdetector_4(log(raw{1}),nucr,blobthreshold,debrisarea);
                case 'single'
                    nuc_mask=threshmask(raw{1},blurradius);
                    nuc_mask=markershed_filter(nuc_mask,round(nucr*2/3),6);
                case 'double'
                    nuc_mask=threshmask_adapt(raw{1},blurradius);
                    nuc_mask=markershed_filter(nuc_mask,round(nucr*2/3),6);
                    nuc_mask=secondthresh(raw{1},blurradius,nuc_mask,boulderarea*2);
            end
            whole_mask = nuc_mask;
            nuc_mask=bwareaopen(nuc_mask,debrisarea);
            %%% Deflection-Bridging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
            nuc_mask=excludelargeandwarped_3(nuc_mask,boulderarea,settings.soliditythresh);
        else
            %nuc_mask=threshmask_adapt(raw{1},blurradius);
            nuc_mask=threshmask(raw{1},blurradius);
            whole_mask = nuc_mask;
            %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lastgoodframe=find(badframes==0,1,'last');
            [reljitx,reljity]=registerimages(imfill(extractmask(regheight,regwidth),'holes'),nuc_mask(regheight,regwidth));
            jitters(f,:)=jitters(lastgoodframe,:)+[reljitx,reljity];
            switch settings.secondsegmethod
                case 'log'
                    nuc_mask=blobdetector_4(log(raw{1}),nucr,blobthreshold,debrisarea);
                case 'double'
                    nuc_mask=threshmask_adapt(raw{1},blurradius);
                    nuc_mask=markershed_filter(nuc_mask,round(nucr*2/3),6);
                    nuc_mask=secondthresh(raw{1},blurradius,nuc_mask,boulderarea*2);
                case 'apriori'
                    [nuc_mask,marker_mask]=apriori_markermask(nuc_mask,nuc_center,jitters(f,:));
            end
            nuc_mask=bwareaopen(nuc_mask,debrisarea);
        end
        %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask=imclearborder(nuc_mask);
        %%% Only for DeBug Mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if debug_mode
            keyboard;
        end
        %{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray((raw{1})));
tempframe(:,:,2)=extractmask;;
%tempframe(:,:,3)=marker_mask;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
hist(nuc_area,100);

anti_mask=bwareaopen(nuc_mask,debrisarea);
temp_mask=nuc_mask-anti_mask;
extractmask=bwmorph(temp_mask,'remove');

anti_mask=bwareaopen(nuc_mask,1000);
extractmask=bwmorph(anti_mask,'remove');
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        compression=settings.compression;

        for j=1:length(names)
            if settings.localbg(j)
                mask_total = whole_mask | foreground{j};
                nanmask{j}=imdilate(mask_total,strel('disk',nucr/2));
                nanmaskcyto{j}=imdilate(mask_total,strel('disk',nucr*2));
            else
                nanmask{j}=imdilate(whole_mask,strel('disk',nucr/2));
                nanmaskcyto{j}=imdilate(whole_mask,strel('disk',nucr*2));
            end
            
            blur=imfilter(raw{j},fspecial('disk',settings.blurradius),'symmetric');
            %blur=raw{j};
            switch settings.bgsubmethod{j}
                case 'global nuclear'
                    real{j}=bgsubmasked_global_NR(blur,nanmask{j},1,compression,25);
                case 'global cyto'
                    real{j}=bgsubmasked_global_2(blur,nanmaskcyto{j},1,compression,25);
                case 'tophat'
                    real{j}=imtophat(blur,strel('disk',nucr,0));
                case 'semi-local nuclear'
                    real{j}=bgsubmasked_global_2(blur,nanmask{j},11,compression,25);
            end
        end
        
        %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [nuc_label,numcells]=bwlabel(nuc_mask);
        nuc_info=struct2cell(regionprops(nuc_mask,real{1},'Area','Centroid','MeanIntensity')');
        nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
        %%%%%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mednuc=median(nuc_area);
        if i>firstgoodindex+1 && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)
            fprintf('badframe: frame %0.0f\n',f);
            badframes(f)=1;
            extractmask=bwmorph(nuc_mask,'remove');
            if maskwrite
                imwrite(uint16(extractmask),[maskdir,namenucedge,num2str(f),'.tif']);
            end
            continue;
        end
        blurthreshhigh=1.1*mednuc;
        blurthreshlow=0.8*mednuc;
        numthresh=0.5*numcells;
        nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
        nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
        %%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mass=nuc_density.*nuc_area;
        curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i>firstgoodindex
            nuc_center(:,1)=nuc_center(:,1)+jitters(f,1);
            nuc_center(:,2)=nuc_center(:,2)+jitters(f,2);
            %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
            debugpackage={extractmask,jitters(lastgoodframe,:),[reljitx,reljity]};
            %%% track & correct merges (update centers, masses and labels) %%%%
            [tracedata,curdata,tracking,nuc_label]=...
                adaptivetrack_9(f,lastgoodframe,f,tracedata,curdata,tracking,real{1},nuc_label,jitters(f,:),trackparams,debugpackage);
            badframes(f)=0;
        end
        %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        extractmask=bwmorph(nuc_label,'remove');
        if maskwrite
            imwrite(uint16(extractmask),[maskdir,namenucedge,num2str(f),'.tif']);
        end
        %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cellid=find(~isnan(curdata(:,1)));
        numlivecells=numel(cellid);
        curdata=curdata(cellid,:);
        nuc_center=curdata(:,[1 2]);
        nuc_area=curdata(:,3);
        nuc_mass=curdata(:,4);
        nuc_info=regionprops(nuc_label,'PixelIdxList','Centroid');
        nanvec=ones(numlivecells,1)*NaN;
        
        %Initialize storage variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sigringmedian{1}=[];
        localbg{1}=[];
        
        for j=1:length(names)
            sigmean{j}=nanvec;
            sigmedian{j}=nanvec;
            %sigring75th{j}=nanvec;
            if settings.ringcalc(j)
                sigringmedian{j}=nanvec;
            end
            if settings.localbg(j)
                localbg{j}=nanvec;
            end
        end
        
        if any(settings.ringcalc)
            ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,real{2});
            %{
            %%% debugging: view images %%%%%%%%%%
            temp_ring = ring_label > 0;
            extractmask=bwmorph(temp_ring,'remove');
            signal = real{2};
            intensities = signal.*temp_ring;
            tempframe=zeros(height,width,3);
            tempframe(:,:,1)=imadjust(mat2gray(real{2}));
            tempframe(:,:,2)=extractmask;;
            %tempframe(:,:,3)=marker_mask;
            figure,imshow(tempframe);
            %}
            ring_info=regionprops(ring_label,'PixelIdxList');
        end
        
        
        %Measure features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j=1:length(names)
            real_temp=real{j};
            real_masked=real_temp;
            real_masked(nanmask{j})=NaN;
            for n=1:numlivecells
                cc=cellid(n);
                sigmean{j}(n)=mean(real_temp(nuc_info(cc).PixelIdxList));
                sigmedian{j}(n)=median(real_temp(nuc_info(cc).PixelIdxList));
                if settings.localbg(j)
                    localbg{j}(n)=getlocalbg(real_masked, nuc_info(cc).Centroid,100, 25, .15);
                end
                if settings.ringcalc(j)
                    if cc>numel(ring_info)
                        break;
                    end
                    ringall=real_temp(ring_info(cc).PixelIdxList);
                    ringall(ringall>prctile(ringall,98))=[];
                    %sigring75th{j}(n)=prctile(ringall,75);
                    ringforeground=ringall(ringall>settings.ringthresh(j));
                    if numel(ringforeground)<50
                        ringforeground=ringall;
                    end
                    if numel(ringall)>=50
                        sigringmedian{j}(n)=nanmedian(ringforeground);
                    end
                end
             end
        end
        
        
        %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tracedata(cellid,f,:)=...
            [nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,...
            cell2mat(sigmean),cell2mat(sigmedian),cell2mat(localbg),...
            cell2mat(sigringmedian)];
        if maxcellnum-max(cellid)<blocksize
            tempdata=ones(blocksize,EF,parameternum)*NaN;
            temptrack=ones(blocksize,5)*NaN;
            tracedata=[tracedata;tempdata];
            tracking=[tracking;temptrack];
            maxcellnum=maxcellnum+blocksize;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toc(timeframe);
    end
    [tracedata,genealogy,jitters]=postprocessing_nolinking(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
    %%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
    names=output_names(settings.signals,settings.localbg, settings.ringcalc);
    save([settings.projectpath,'settings.mat'],'settings','names'); %saves all the settings to the Data folder. will overwrite with the most recent

    toc(timetotal);
    clear all;
end
end